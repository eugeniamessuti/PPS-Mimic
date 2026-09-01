# Methods: `2_Mimic_analysis_simp.R` (Simplified — no cohort LOS filter, no admin censoring)

> **Revision note:** Per supervisor feedback, this version removes two elements of the original method:
> (1) the cohort-inclusion filter requiring LOS > 48 hours, and (2) administrative censoring at day 60.
> Both changes are applied together, across all analyses (Reference, PPS, CDC). The HAI definition
> (culture-confirmed infection occurring >48h after admission) is **unchanged**. Two decisions that
> follow from removing the day-60 cap and the LOS>48h filter — how "no censoring" is implemented, and
> what upper bound now governs the PPS/CDC survey-date sampling window — have since been resolved (see
> Sections 0 and 5) and are no longer open items.

## 0. Data preparation

Before any modeling, the analysis dataset (`0_my_mimic_simp.csv`) is built from three MIMIC-IV tables — `admissions`, `patients`, and `microbiologyevents` — linked via `subject_id` and `hadm_id` (script `0_my_mimic_simp.R`).

**Cohort assembly and variable derivation:**

- **Cohort inclusion**: ~~only hospitalizations with a length of stay longer than 48 hours are retained~~ — **removed**. All hospitalizations are retained regardless of length of stay, including those shorter than 48 hours.
- **Age at admission**: derived as `anchor_age + (year(admittime) - anchor_year)`, MIMIC-IV's de-identification-preserving convention for reconstructing age at a given event from the de-identified `anchor_age`/`anchor_year` fields. *(Unchanged.)*
- **HAI definition**: a hospitalization is flagged as having a hospital-acquired infection if there is at least one microbiology record with a non-missing organism (`org_name`) and a `charttime` more than 48 hours after `admittime`; the HAI time is the earliest such qualifying `charttime` (`slice_min(charttime, n = 1, with_ties = FALSE)` — ties broken by taking the first row only). *(Unchanged, per instruction.)*
- **LOS and time-to-HAI**: computed as `difftime()` in days between `dischtime`/`hai_time` and `admittime`. *(Unchanged.)*
- **Administrative censoring at day 60**: ~~already applied at this data-preparation stage~~ — **removed**. `los_time` and `hai_time` are carried forward at their actual recorded values, with no cap. **Resolved:** censoring is removed in full — every admission in MIMIC-IV has a recorded `dischtime`, so `cens_los` is set to `1` for every row and `cens_los_time` is the actual (uncapped) LOS; similarly `cens_hai` is `1` whenever an HAI was ever recorded and `cens_hai_time` is the actual time-to-HAI (or the actual LOS, for patients discharged without HAI). No other censoring mechanism applies.

**Data-quality exclusions applied inside the modeling script (`2_Mimic_analysis_simp.R`):**

- Only complete cases on `age` and `sex` are retained before modeling. *(Unchanged.)* **546,028** admissions had complete `age`/`sex` covariates at this stage.
- **New exclusion — non-positive length of stay** (`los_time <= 0`, i.e. `dischtime <= admittime` in the raw MIMIC-IV data): **180** admissions were excluded on this criterion. This is a known MIMIC data-quality artifact rather than a real same-day negative-duration stay — the MIT-LCP `mimic-code` maintainers report that whenever `admittime > dischtime`, both timestamps fall on the same calendar day, consistent with administrative entry errors (e.g. swapped admission/discharge fields, or "12:00" typed for midnight "00:00") rather than an actual clinical event (see `github.com/MIT-LCP/mimic-code`, issue #209). Under the original (non-simplified) script this was silently absorbed by the LOS>48h cohort filter; since that filter has been removed, these admissions surface directly and are excluded here instead, as a data-quality step independent of the removed filter. Without a valid, non-corrupted discharge time there is no usable time-at-risk interval for these patients, so they cannot contribute to the model.
- Timestamp-consistency exclusion: hospitalizations where the recorded HAI time falls *after* the discharge time (`hai_time > los_time`) are dropped, for the same data-quality reason as before (delayed lab-reporting timestamps). *(Unchanged — this is a data-quality exclusion, not the cohort-inclusion filter that was removed.)* Of 95,572 patients recorded as infected at this stage, **88** were excluded on this criterion.
- Discharge times that are observed events get a `+0.001` day nudge (`exit_time`) to avoid exact time-ties with same-day infection or entry rows, since `msm` requires strictly increasing times per subject. *(Unchanged.)*
- An additional `+1e-6 × row_number()` jitter is applied within each subject's rows in the reference and PPS panels, for the same tie-breaking reason. *(Unchanged.)*
- Parallelization (`parallel`, `doParallel`, `foreach`) is a computational detail for speed and does not change the statistical method. *(Unchanged.)*

**Resulting cohort (`mimic_0`):**

```
Admissions with complete age/sex covariates       : 546,028
Excluded (los_time <= 0, dischtime <= admittime)  :    180
Excluded (hai_time > los_time)                    :     88
---------------------------------------------------------
Final analysis cohort                             : 545,760

  HAI (state 2)                    :  95,484  (17.5%)
  Discharged, no HAI (state 3)     : 450,276  (82.5%)
  Admin censored, no HAI           :       0  ( 0.0%)  -- structurally
                                                            empty under
                                                            this scheme
```

## 1. Model

We fit a **continuous-time, time-homogeneous multi-state Markov model with constant transition intensities**, estimated by maximum likelihood using the `msm` R package (Jackson, C.H. (2011), *Multi-State Models for Panel Data: The msm Package for R*, Journal of Statistical Software, 38(8)). *(Unchanged.)*

Each transition intensity is assumed **constant over the follow-up window**.

**Covariate effects** (age, sex) enter as a **log-linear multiplicative model on the transition intensities**: each covariate shifts the constant rate by a multiplicative factor `exp(β)`, which `msm` reports as a **hazard ratio** via `hazard.msm()`. *(Unchanged.)*

## 2. Multi-State structure

A 3-state, irreversible **illness–death model**:

- **State 1** = Admitted (at risk, uninfected)
- **State 2** = HAI (hospital-acquired infection acquired)
- **State 3** = Discharged (absorbing state)

Allowed transitions (from the `qmat` object): 1→2, 1→3, 2→3. No reverse transitions are permitted. *(Unchanged.)*

Three rate parameters are estimated: **λ₁₂**, **λ₁₃**, **λ₂₃**. *(Unchanged.)*

## 3. Study design features

- **Entry into the risk set is now at day 0** (`entry_time <- 0`), not day 2. In the original method, delayed entry at day 2 was tied to excluding admissions shorter than 48 hours; since that cohort-inclusion filter has been removed, there is no longer a design basis for a day-2 landmark, so all patients now enter the risk set at admission (day 0).
- **Right (administrative) censoring at day 60**: **removed**. Patients are followed to their actual recorded discharge time (or actual recorded infection time), with no fixed follow-up horizon.
- The comparison of a full-cohort reference analysis against PPS and CDC survey-based subsampling protocols is otherwise unchanged in structure (Sections 5–7), but is now run on the simplified cohort/censoring scheme described above.

## 4. Reference (full-cohort) analysis

This is the benchmark against which the candidate estimation methods are compared. It uses the complete MIMIC data (now including all hospitalizations, any length of stay) with the exact recorded times of infection and discharge for every patient.

- `crudeinits.msm()` computes non-optimized starting values for the Q-matrix from raw observed transition counts.
- `msm(state ~ time, ..., obstype = obs_type)` fits the model by maximum likelihood.
- `obstype = 2` is used throughout the reference panel, as before.

**Model without covariates** and **model with covariates** (`~ age + sex`) are fit as before, with λ's, delta-method CIs, and hazard ratios extracted the same way (`qmatrix.msm()`, `hazard.msm()`). *(Methodologically unchanged — only the underlying data now reflects no LOS filter and no admin censoring.)*

**Delta-method CI for the reference λ₂₃/λ₁₃ ratio**: unchanged procedure.

## 5. Simulating a cross-sectional survey (100 Monte Carlo runs)

Both the PPS and CDC sections repeat, 100 times, the design: draw a random survey date per patient and keep only patients still hospitalized on that date (`los_time > PPS_date`).

**Resolved:** the survey date was previously drawn as `PPS_date ~ Uniform(pps_date_min, admin_cens)`, using `admin_cens = 60` as the upper bound of the sampling window. With `admin_cens` removed, the upper bound is now data-driven: the maximum observed LOS in the cleaned cohort, `pps_date_max <- ceiling(max(mimic_0$los_time))`. In the current cohort this evaluates to **516**, so `PPS_date ~ Uniform(0, 516)`. This same bound is reused as the horizon for the Stacked Probability Plot (Figure B) time grid.

**PPS_date clamp (CDC protocol only).** In the CDC panel construction, a non-prevalent patient's "known uninfected as of the survey" row is timestamped at `PPS_date`. Since `entry_time` is now 0, the clamp condition (`PPS_date <= entry_time`) becomes `PPS_date <= 0`, which — given `PPS_date` is drawn from a continuous distribution starting above 0 — will essentially never trigger in practice. The clamp logic itself does not need to change, but it becomes a near-vacuous edge case rather than a meaningfully active correction.

Within each protocol, two weighting variants are fit per run:
- **Unweighted** (`fit_uw`): every subject contributes equally.
- **Weighted** (`fit_w`): `msm(..., subject.weights = weight)`, `weight <- 1 / los_time`. *(Unchanged rationale — length-biased sampling correction.)*

## 6. The PPS protocol

**Length-sampling bias, but no interval censoring.** Unchanged in concept: once a patient is sampled, the panel is built from the full, exact recorded infection and discharge times, and every row is coded `obstype = 2`.

**Patient-panel structure.** With entry at day 0 instead of day 2, and no administrative cap, the illustrative panel becomes:

**Resolved:** the original patient types 3 and 4 ("LOS censored") were defined as patients not yet discharged by day 60. Since administrative censoring has been removed in full, and admissions with a corrupted (non-positive) length of stay are excluded upstream (Section 0), every remaining patient's discharge is genuinely, fully observed — there is no "censored" case left to represent. Only the two patient types below exist under this simplified scheme.

| ID | Time | State | Obstype |
|---|---|---|---|
| 1 | 0  | 1 | 2 |
| 1 | Time to infection | 2 | 2 |
| 1 | LOS + 0.001 | 3 | 2 |
| 2 | 0 | 1 | 2 |
| 2 | LOS + 0.001 | 3 | 2 |

Every row is `obstype = 2`: infection and discharge times are treated as exactly known throughout.

**Two model fits.** Unchanged: one **weighted** (`fit_w`) and one **unweighted** (`fit_uw`) per run.

## 7. CDC

**Definition.** Unchanged in concept — same survey-day sampling design as PPS, plus interval censoring: only patients already infected at the time of the PPS are considered infected; HAI occurring after the PPS date is not gathered by the protocol.

**Patient-panel structure.** With entry at day 0:

**Resolved:** as in Section 6, the original "ICU length of stay was censored" patient types (using state `99`/`obstype = 1` at administrative censoring) no longer have a design basis — administrative censoring has been removed in full, so every remaining patient's eventual discharge is fully observed. The `censor = 99` argument and the state-`99` placeholder have been removed from the CDC panel-construction code accordingly. Only the two patient types below exist under this simplified scheme.

| ID | Time | State | Observation type |
|---|---|---|---|
| 1 | 0 | 1 | 2 |
| 1 | Time to infection | 2 | 2 |
| 1 | LOS + 0.001 | 3 | 2 |
| 2 | 0 | 1 | 2 |
| 2 | PPS date | 1 | 2 |
| 2 | LOS + 0.001 | 3 | 3 |

- **Patient 1** (infected as of PPS): exact infection time known, coded `obstype = 2`, followed by exact discharge time, `obstype = 2`.
- **Patient 2** (uninfected as of PPS): a row at the PPS date records "uninfected as of the survey" (`obstype = 2`); the discharge row is `obstype = 3` (exact discharge time known, but whether the patient passed through HAI between survey and discharge is not) — this is the interval-censoring correction, unchanged in concept.

**Four models.** Unchanged: Unweighted + IC-naive, Weighted + IC-naive, Unweighted + IC-corrected, Weighted + IC-corrected (correct method).

**IC-naive panel.** Same two patient types, coded without the interval-censoring correction (every row `obstype = 2`):

| ID | Time | State | Observation type |
|---|---|---|---|
| 1 | 0 | 1 | 2 |
| 1 | Time to infection | 2 | 2 |
| 1 | LOS + 0.001 | 3 | 2 |
| 2 | 0 | 1 | 2 |
| 2 | PPS date | 1 | 2 |
| 2 | LOS + 0.001 | 3 | 2 |

## 8. Output

Figures A–D and the excess-LOS table are unchanged in definition (hazard-rate boxplots, stacked probability plot, covariate HR boxplots, λ₂₃/λ₁₃ ratio boxplot, excess LOS formula). All are computed on the simplified cohort/censoring scheme described above.

```
Excess LOS = (1 / (λ₁₂ + λ₁₃)) × (λ₁₃ / λ₂₃ − 1)
```
