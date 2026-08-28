# Methods: `2_Mimic_msm_analysis.R`

## 0. Data preparation

Before any modeling, the analysis dataset (`my_mimic.csv`) is built from three MIMIC-IV tables — `admissions`, `patients`, and `microbiologyevents` — linked via `subject_id` and `hadm_id` (script `0_my_mimic.R`).

**Cohort assembly and variable derivation:**

- **Cohort inclusion**: only hospitalizations with a length of stay longer than 48 hours (`dischtime > admittime + 48 hours`) are retained at this stage.
- **Age at admission**: derived as `anchor_age + (year(admittime) - anchor_year)`, MIMIC-IV's de-identification-preserving convention for reconstructing age at a given event from the de-identified `anchor_age`/`anchor_year` fields.
- **HAI definition**: a hospitalization is flagged as having a hospital-acquired infection if there is at least one microbiology record with a non-missing organism (`org_name`) and a `charttime` more than 48 hours after `admittime`; the HAI time is the earliest such qualifying `charttime` (`slice_min(charttime, n = 1, with_ties = FALSE)` — ties broken by taking the first row only).

- **LOS and time-to-HAI**: computed as `difftime()` in days between `dischtime`/`hai_time` and `admittime`.
- **Administrative censoring at day 60** is already applied at this data-preparation stage, before the modeling script: `cens_los`/`cens_los_time` cap LOS at 60 days, and `cens_hai`/`cens_hai_time` cap time-to-HAI at 60 days in the same way — these are the columns read directly into `mimic_0_all` in the modeling script.

**Further preparation inside the modeling script (`2_Mimic_msm_analysis.R`):**

- Only complete cases on `age` and `sex` are retained before modeling.
- Timestamp-consistency exclusion: hospitalizations where the recorded HAI time falls *after* the (censored) discharge time (`hai_time > los_time`) are dropped. This arises in MIMIC because microbiology result timestamps are not constrained to precede `dischtime` — a delayed lab-reporting artifact — which makes the ordering of infection and discharge unreliable for these patients. The number excluded is logged at run time. Exact same-day ties (`hai_time == los_time`) are not excluded here; they are instead resolved by the `+0.001` exit-time nudge below.
- Discharge times that are observed events get a `+0.001` day nudge (`exit_time`) to avoid exact time-ties with same-day infection or entry rows, since `msm` requires strictly increasing times per subject.
- An additional `+1e-6 × row_number()` jitter is applied within each subject's rows in the reference and PPS panels, for the same tie-breaking reason.
- Parallelization (`parallel`, `doParallel`, `foreach`) is a computational detail for speed across the 100 runs and does not change the statistical method.

## 1. Model

We fit a **continuous-time, time-homogeneous multi-state Markov model with constant transition intensities**, estimated by maximum likelihood using the `msm` R package (Jackson, C.H. (2011), *Multi-State Models for Panel Data: The msm Package for R*, Journal of Statistical Software, 38(8)). This package fits continuous-time Markov (or hidden Markov) multi-state models by maximum likelihood, where observations of the process can be made at arbitrary times or the exact times of transition can be known, and covariates can be fitted to the transition intensities.

Each transition intensity is assumed **constant over the follow-up window**.

**Covariate effects** (age, sex) enter as a **log-linear multiplicative model on the transition intensities**: each covariate shifts the constant rate by a multiplicative factor `exp(β)`, which `msm` reports as a **hazard ratio** via `hazard.msm()`. The **baseline hazard** is parametric and homogeneous (piecewise-constant): the transition rates qᵣₛ are held constant within the observation window rather than estimated nonparametrically, as they would be under a Cox-type model.

## 2. Multi-State structure

A 3-state, irreversible **illness–death model**:

- **State 1** = Admitted (at risk, uninfected)
- **State 2** = HAI (hospital-acquired infection acquired)
- **State 3** = Discharged (absorbing state)

Allowed transitions (from the `qmat` object): 1→2, 1→3, 2→3. No reverse transitions are permitted.

Three rate parameters are estimated:
- **λ₁₂** — instantaneous rate of acquiring HAI
- **λ₁₃** — instantaneous rate of discharge without ever acquiring HAI
- **λ₂₃** — instantaneous rate of discharge after acquiring HAI

## 3. Study design features

- **Left truncation / delayed entry at day 2** (`entry_time <- 2`): patients enter the risk set at day 2, consistent with the ≥48h clinical convention for classifying an infection as hospital-*acquired*.
- **Right (administrative) censoring at day 60** (`admin_cens <- 60`): patients not yet discharged or infected by day 60 are censored there. This is consistent with the CDC methodology proposed by Magill et al. (Magill SS, Edwards JR, Bamberg W, et al. Multistate point-prevalence survey of health care–associated infections. N Engl J Med. 2014;370(13):1198–1208).
- We compare a full-cohort reference analysis against two survey-based subsampling protocols (PPS and CDC). Both protocols share the same underlying sampling design (Section 5) and are each run under two weighting schemes (unweighted and weighted). The CDC protocol additionally introduces interval censoring, so it is run under two panel-coding schemes (IC-naive and IC-corrected) crossed with the two weighting schemes, yielding four CDC models in total, versus two for PPS.

## 4.Reference (full-cohort) analysis

This is the benchmark against which the three candidate estimation methods are compared. It uses the complete MIMIC data with the exact recorded times of infection and discharge for every patient (no survey subsampling).

- `crudeinits.msm()` computes non-optimized starting values for the Q-matrix from raw observed transition counts, passed to `msm()` as starting values for the likelihood optimizer.
- `msm(state ~ time, ..., obstype = obs_type)` fits the model by maximum likelihood.
- `obstype = 2` is used throughout the reference panel. In `msm`, obstype 2 means an exact transition time, with the state at the previous observation retained until the current observation — i.e., each patient's exact transition times are treated as known, which matches the reference data.

**Model without covariates — the λ's and their 95% CIs.** From this model we obtain the three fitted constant transition rates λ₁₂, λ₁₃, λ₂₃ via `qmatrix.msm(fit, ci = "none")`. Internally, `msm` estimates the transition intensities on the **log scale**, and a covariance matrix for these log-scale estimates is obtained from the Hessian of the maximized log-likelihood. The 95% confidence intervals reported alongside the point estimates (`qmatrix.msm(fit, ci = "delta")`) are obtained by the **delta method**: the log-scale covariance matrix is propagated through the exponential transformation back to the natural (rate) scale, assuming asymptotic normality on the log scale, and the resulting interval is then back-transformed to give the 95% CI on the rate scale.

**Model with covariates — the hazard ratios.** A second model adds `covariates = ~ age + sex`. Here, each covariate enters as a linear effect on the log-transition-intensity scale (i.e. `log(q_rs) = log(q_rs,0) + β·covariate`), so that `exp(β)` is the multiplicative effect of the covariate on the rate — the **hazard ratio**. `hazard.msm()` extracts these exponentiated coefficients per transition, together with their 95% confidence intervals, using the **same delta-method approach** described above: the Hessian-based covariance of the log-scale coefficients is propagated through the exponential transformation to obtain approximate standard errors and confidence limits on the hazard-ratio scale.

**Delta-method CI for the reference λ₂₃/λ₁₃ ratio.** For the reference analysis, a 95% CI for the ratio λ₂₃/λ₁₃ is derived by:
1. Taking `msm`'s delta-method CIs for λ₂₃ and λ₁₃ separately.
2. Back-calculating an implied log-scale standard error for each rate from its CI width: `(log(U) − log(L)) / (2 × 1.96)`.
3. Combining the two log-SEs as `sqrt(SE_log(λ23)² + SE_log(λ13)²)` to get the SE of `log(ratio)`, then exponentiating back for the final CI.


## 5. Simulating a cross-sectional survey (100 Monte Carlo runs)

Both the PPS and CDC sections repeat, 100 times, the following design: draw a random survey date per patient, `PPS_date ~ Uniform(pps_date_min, admin_cens)`, and keep only patients still hospitalized on that date (`los_time > PPS_date`). This simulates a one-day, cross-sectional survey ("point prevalence survey") as the source of infection-status information rather than continuous retrospective surveillance. Each of the 100 runs re-draws survey timing at random, producing a Monte Carlo sampling distribution of the estimates.

**`pps_date_min` is decoupled from `entry_time`.** `PPS_date` is drawn from `Uniform(0, admin_cens)` — i.e. `pps_date_min <- 0`, not `entry_time` (day 2). This means a sampled `PPS_date` can legitimately fall *before* a patient's day-2 entry into the risk set.

**PPS_date clamp (CDC protocol only).** In the CDC panel construction (`build_cdc_panel()` and `build_cdc_panel_no_ic()`), a non-prevalent patient's "known uninfected as of the survey" row is timestamped at `PPS_date`. If `PPS_date <= entry_time`, this row would fall at or before the patient's day-2 entry row, which `msm` rejects as a non-increasing (backward) transition time. Since there is no information about a patient's state before `entry_time` in any case — they are not yet in the risk set — such a `PPS_date` carries no more information than "uninfected at day 2." The code therefore clamps this row's time to `entry_time + 1e-4` whenever `PPS_date <= entry_time`. This clamp is applied identically in both the IC-corrected and IC-naive panels. No equivalent clamp is needed for prevalent (already-infected) patients, since `hai_time` is by construction always ≥ `entry_time` under the ≥48h HAI convention.

Within each protocol, two weighting variants are fit per run:
- **Unweighted** (`fit_uw`): every subject contributes equally to the likelihood.
- **Weighted** (`fit_w`): `msm(..., subject.weights = weight)`, where `weight <- 1 / los_time`, an **inverse-length-of-stay weight**. Rationale: in a one-day cross-sectional survey, longer-stay patients have proportionally higher odds of being captured on any given survey day than short-stay patients — a length-biased sampling problem. Weighting each sampled patient by 1/(length of stay) is a design-based correction for that over-representation.

## 6 the PPS protocol
 
 
**Length-sampling bias, but no interval censoring.** This sampling mechanism introduces a **length-biased sampling** problem: a patient with a longer hospital stay has a higher probability of being "caught" by any given survey day than a patient with a short stay, so the pool of sampled patients systematically over-represents long-stay patients relative to the full cohort — this is the bias the 1/LOS weighting (Section 5) is designed to correct. However, once a patient is sampled, it does **not** introduce interval censoring: the panel is built from the full, exact recorded infection and discharge times for every sampled patient , and every row is coded `obstype = 2`. No censoring code is used. All infections are fully gathered/observed by the protocol — the only distortion introduced by this design is in *which* patients are selected, not in what is known about them once selected.
 
**Patient-panel structure.** The table below shows how four illustrative patients are represented as rows for `msm`:

4 types of patients were defined:
1.	Patient that acquired a HAI during their ICU stay and whose ICU length of stay was not censored
2.	Patient that didn´t acquired a HAI during their ICU stay and whose ICU length of stay was not censored
3.	Patient that acquired a HAI during their ICU stay and whose ICU length of stay was censored
4.	Patient that didn´t acquired a HAI during their ICU stay and whose ICU length of stay was censored

 
| ID | Time | State | Obstype |
|---|---|---|---|
| 1 | 2  | 1 | 2 |
| 1 | Time to infection | 2 | 2 |
| 1 | ICU LOS + 0.001 | 3 | 2 |
| 2 | 2 | 1 | 2 |
| 2 | ICU LOS + 0.001 | 3 | 2 |
| 3 | 2 | 1 | 2 |
| 3 | Time to infection | 2 | 2 |
| 3 | 60 | 2 | 2 |
| 4 | 2 | 1 | 2 |
| 4 | 60 | 1 | 2 |
 
 
Every row across all four patients is `obstype = 2`: infection and discharge times are treated as exactly known throughout, consistent with there being no interval censoring in this protocol.
 
**Two model fits.** For the PPS protocol, two models are fit per run: one **weighted** (`fit_w`, `subject.weights = 1/los_time`) and one **unweighted** (`fit_uw`).

## 7. CDC

**Definition.** The CDC method refers to the way in which data is collected — the same sampling design as for the PPS protocol (a random survey day, `PPS_date`, drawn per patient, retaining only patients still hospitalized on that day). Additionally, the CDC protocol introduces **interval censoring**: only patients who already had a HAI at the time of the PPS are considered infected. Any HAI that occurs *after* the PPS date is not gathered by the protocol — from the survey's point of view, such a patient looks uninfected as of the survey, and what happens to them afterward (whether or when they become infected before discharge) is not observed directly.

![See 1 Tornado Plot CDC ][def]

**Patient-panel structure.** Four patient types are defined:

1. Patient with HAI at time PPS whose ICU length of stay was not censored
2. Patient without HAI at time PPS whose ICU length of stay was not censored
3. Patient with HAI at time PPS whose ICU length of stay was censored
4. Patient without HAI at time PPS whose ICU length of stay was censored

| ID | Time | State | Observation type |
|---|---|---|---|
| 1 | 2 | 1 | 2 |
| 1 | Time to infection | 2 | 2 |
| 1 | LOS + 0.001 | 3 | 2 |
| 2 | 2 | 1 | 2 |
| 2 | PPS date | 1 | 2 |
| 2 | LOS + 0.001 | 3 | 3 |
| 3 | 2 | 1 | 2 |
| 3 | Time to infection | 2 | 2 |
| 3 | 60 | 2 | 2 |
| 4 | 2 | 1 | 2 |
| 4 | PPS date | 1 | 2 |
| 4 | 60 | 99 | 1 |

- **Patient 1** (infected as of PPS, not censored): the exact infection time is known (already infected at the survey), so it is coded `obstype = 2` as usual, followed by an exact discharge time, also `obstype = 2`.
- **Patient 2** (uninfected as of PPS, not censored): a row is added *at the PPS date itself*, coded `obstype = 2`, recording "uninfected as of the survey." The subsequent discharge row is coded `obstype = 3` — the exact discharge time is known, but whether the patient passed through HAI (state 2) at some point between the survey and discharge is not: this is the interval-censoring correction, since the protocol only tells us the patient's state at the survey and at the moment of discharge, not in between.
- **Patient 3** (infected as of PPS, censored): the exact infection time is known, and the patient is still known to be infected (state 2) at administrative censoring, so this final row stays `obstype = 2`.
- **Patient 4** (uninfected as of PPS, censored): a row at the PPS date records "uninfected as of the survey" (`obstype = 2`), but by administrative censoring we no longer know the patient's state — they may have become infected in the interim — so the state is recoded to `99` (an internal placeholder) with `obstype = 1`, `msm`'s coding for "state unknown, observed only at an arbitrary time."

**Four models.** For the CDC protocol we fit **four models**, in order to observe how the result would be biased relative to the reference line if we chose to ignore the length-sampling bias and/or the interval censoring:

- Unweighted + IC-naive
- Weighted + IC-naive
- Unweighted + IC-corrected
- Weighted + IC-corrected — the correct method

**IC-naive panel.** For the models that do not address interval censoring, patients are instead defined according to the following table:

| ID | Time | State | Observation type |
|---|---|---|---|
| 1 | 2 | 1 | 2 |
| 1 | Time to infection | 2 | 2 |
| 1 | LOS + 0.001 | 3 | 2 |
| 2 | 2 | 1 | 2 |
| 2 | PPS date | 1 | 2 |
| 2 | LOS + 0.001 | 3 | 2 |
| 3 | 2 | 1 | 2 |
| 3 | Time to infection | 2 | 2 |
| 3 | 60 | 2 | 2 |
| 4 | 2 | 1 | 2 |
| 4 | PPS date | 1 | 2 |
| 4 | 60 | 1 | 2 |

Here, every row is coded `obstype = 2`: Patient 2's discharge is treated as an exact, fully known transition rather than the ambiguous `obstype = 3` used in the IC-corrected panel, and Patient 4's final observation is treated as an exact, known state (`1`, uninfected) at censoring rather than the unknown/censored placeholder (`99`, `obstype = 1`). This is the "naive" coding: it ignores the fact that, in a real survey design, a patient's infection status between the PPS date and their eventual discharge or censoring is not actually observed.

The `msm` function was then fitted the same way as in the Reference and PPS analyses, using these patient definitions.

## 8. Output

**Figure A** compares the constant hazard rates across all methods. These are boxplots, across the 100 Monte Carlo runs, of the three fitted constant hazard rates λ₁₂, λ₁₃, λ₂₃, with the full-cohort reference estimate and its delta-method 95% CI overlaid. The baseline against which all methods are compared is the full-cohort Reference analysis described in Section 2.

**Figure B** is the **Stacked Probability Plot (SPP)** — the exact analytic solution of the 4×4 transition-intensity matrix (states: Admitted / HAI / Discharged-no-HAI / Discharged-post-HAI) via matrix exponentiation (`expm::expm(Q4 * t)`), evaluated at the mean of the 100 runs' hazard estimates. This is the same underlying calculation as the model's cumulative risk: state-occupation probabilities at any follow-up time *t* are obtained directly from the fitted model as `P(t) = exp(tQ)`. The reference panel is the one using the λ's from the Reference analysis (Section 2).

**Figure C** shows boxplots of the fitted covariate hazard ratios for age and sex, per transition, on a log scale, again with the Reference analysis (Section 2) hazard ratios and their delta-method 95% CI overlaid as the comparison baseline.

**Figure D** shows a boxplot of the ratio λ₂₃ / λ₁₃ across runs, again against the Reference analysis (Section 2) estimate of the same ratio.

Finally, a table is attached with the excess length of stay attributable to HAI for each method of the CDC protocol, calculated following this formula:

```
Excess LOS = (1 / (λ₁₂ + λ₁₃)) × (λ₁₃ / λ₂₃ − 1)
```


[def]: pps_timing_diagram.png