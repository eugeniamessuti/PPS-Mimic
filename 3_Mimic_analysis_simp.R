############################################################
# MIMIC HAI ANALYSIS — PPS & CDC PROTOCOLS, 100 RUNS EACH
# SIMPLIFIED VERSION per supervisor request:
#   - Cohort inclusion filter (LOS > 48h) REMOVED upstream in
#     0_my_mimic_simp.R -- all hospitalizations are included here.
#   - Administrative censoring at day 60 REMOVED: entry_time is now 0
#     (day of admission), and there is no follow-up cap -- every patient
#     is followed to their actual recorded discharge time.
#   - HAI definition (>48h post-admission) is UNCHANGED.
#   - NEW data-quality exclusion: admissions with los_time <= 0
#     (dischtime <= admittime -- a MIMIC timestamp-entry artifact, see
#     MIT-LCP mimic-code issue #209) are excluded in Section 1. These were
#     previously masked by the removed LOS>48h cohort filter; without that
#     filter they surface directly and must be excluded on their own merits.
#   - Because there is no day-60 cap, the CDC protocol's "administratively
#     censored, unknown state" branch (state 99 / obstype 1 / censor = 99)
#     is REMOVED: every patient's eventual discharge is fully observed,
#     so that branch never applied and is not needed.
#   - The PPS/CDC survey-date sampling window's upper bound, previously
#     admin_cens = 60, is now data-driven: the maximum observed LOS in
#     the cohort (pps_date_max). This same bound is reused as the
#     horizon for the Stacked Probability Plot time grid.
#   - Since no patient can be censored under this scheme, the Reference
#     and PPS panel builders (build_msm_long_ref, run_one_pps) no longer
#     branch on a "cens" outcome, and the los_status column has been
#     removed entirely (it was always 1 and played no further role).
#
# Illness-death model: 1 = Admitted (at risk), 2 = HAI, 3 = Discharged
#
# FIGURES PRODUCED (per protocol, PPS and CDC):
#   Fig A : Boxplots of constant hazard rates (lambda12, lambda13, lambda23)
#   Fig B : Stacked probability plot from mean hazards (matrix exponential)
#   Fig C : Boxplots of covariate HRs (age, sex) per transition
#   Fig D : Boxplot of ratio lambda23 / lambda13
#
# CDC protocol uses FOUR methods:
#   - Unweighted + IC-naive
#   - Unweighted + IC-corrected
#   - Weighted   + IC-naive
#   - Weighted   + IC-corrected
#
# REQUIRES: my_mimic data.frame already loaded in the environment with
#   columns: hadm_id, gender, age_at_admission, hai, hai_time, los,
#            cens_los, cens_los_time, cens_hai, cens_hai_time
#   (produced by 0_my_mimic_simp.R -- cens_los is always 1, cens_los_time
#   and cens_hai_time are uncapped, under this simplified scheme)
#
# REQUIRES PACKAGES: msm, dplyr, tibble, tidyr, ggplot2, expm,
#                     parallel, doParallel, foreach, patchwork
############################################################

library(msm)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(expm)
library(parallel)
library(doParallel)
library(foreach)
library(patchwork)

# 1. Access my_mimic data --------------------------------------------
setwd("C:/Users/messuti/Documents/MIMIC-IV")
my_mimic <- read.csv("0_my_mimic_simp.csv")
variable.names(my_mimic)


############################################################
# SECTION 0: CONFIGURATION
############################################################

entry_time   <- 0    # SIMPLIFIED: patients enter the risk set at admission
                      # (day 0). Previously day 2, tied to the now-removed
                      # LOS>48h cohort filter -- that design basis no longer
                      # applies, so entry reverts to day 0 for everyone.
pps_date_min <- 0
n_runs       <- 100

# NOTE: admin_cens (previously 60) is removed. There is no follow-up cap;
# patients are followed to their actual recorded discharge time.
# pps_date_max (the upper bound for sampling a survey day, and the SPP
# plot horizon) is computed below, once the analysis cohort is built --
# it is data-driven (max observed LOS), not a fixed design constant.

cat("============================================================\n")
cat("CONFIGURATION (SIMPLIFIED: no LOS filter, no admin censoring)\n")
cat(sprintf("  entry_time   : %d  <-- all patients enter at admission\n", entry_time))
cat("  admin_cens   : none -- patients followed to actual discharge\n")
cat(sprintf("  pps_date_min : %d\n", pps_date_min))
cat("============================================================\n\n")

############################################################
# SECTION 1: BUILD THE ANALYSIS COHORT
############################################################

# No patient can be censored under this scheme: every admission's
# discharge is a fully observed event (cens_los is always 1 in
# 0_my_mimic_simp.R), so a los_status column is not needed here.
mimic_0_all <- my_mimic %>%
  transmute(
    id            = hadm_id,
    age           = age_at_admission,
    sex           = ifelse(gender == "F", 1L, 0L),   # 1 = female, 0 = male (ref)
    hai_status    = as.integer(cens_hai),             # 1 = infected (fully observed)
    hai_time      = ifelse(cens_hai == 1, cens_hai_time, NA_real_),
    los_time      = cens_los_time                     # actual, uncapped discharge time
  )

# ---- Keep only complete cases on the covariates ----
n_complete <- sum(complete.cases(mimic_0_all[, c("age", "sex")]))
mimic_0 <- mimic_0_all[complete.cases(mimic_0_all[, c("age", "sex")]), ]

# ---- Exclude admissions with non-positive length of stay
# (los_time <= 0, i.e. dischtime <= admittime in the raw MIMIC-IV data).
# This is a known MIMIC data-quality issue, not a real same-day negative
# stay: the MIT-LCP mimic-code maintainers report that whenever
# admittime > dischtime, both timestamps fall on the same calendar day,
# consistent with administrative entry errors (e.g. swapped admit/discharge
# fields, or "12:00" typed for midnight "00:00") rather than an actual
# clinical event (see github.com/MIT-LCP/mimic-code, issue #209). Under the
# previous (non-simplified) script this was silently excluded by the
# LOS>48h cohort filter; since that filter has been removed, these
# admissions now surface directly and must be excluded here instead, as a
# data-quality step independent of that removed filter. Without a valid,
# non-corrupted discharge time there is no usable time-at-risk interval for
# these patients, so they cannot contribute to the model. ----
n_before_los_excl <- nrow(mimic_0)
bad_los_ids <- mimic_0$id[mimic_0$los_time <= 0]

cat("--- Data-consistency check: los_time <= 0 (dischtime <= admittime) ---\n")
cat(sprintf("  Admissions before exclusion : %d\n", n_before_los_excl))
cat(sprintf("  Excluded (los_time <= 0)    : %d\n", length(bad_los_ids)))
cat("  (Timestamp-entry artifact; see MIT-LCP mimic-code issue #209.)\n\n")

mimic_0 <- mimic_0[!(mimic_0$id %in% bad_los_ids), ]

# ---- Exclude patients whose recorded HAI time falls after their
# discharge time. This happens in MIMIC because microbiology result
# timestamps aren't constrained to fall before dischtime (delayed lab
# reporting / timestamp inconsistency), making the infection and discharge
# events impossible to order reliably on the same timeline. This is a
# data-quality exclusion, distinct from the (now-removed) LOS>48h cohort
# filter -- it is retained unchanged. ----
n_infected_before <- sum(mimic_0$hai_status == 1)
bad_ids <- mimic_0$id[mimic_0$hai_status == 1 & mimic_0$hai_time > mimic_0$los_time]

cat("--- Data-consistency check: hai_time vs los_time ---\n")
cat(sprintf("  Infected patients before exclusion : %d\n", n_infected_before))
cat(sprintf("  Excluded (hai_time > los_time)     : %d\n", length(bad_ids)))
cat("  (HAI culture resulted after discharge in the source data.)\n\n")

mimic_0 <- mimic_0[!(mimic_0$id %in% bad_ids), ]

# ---- Adjusted exit/discharge time: nudge by +0.001 so a discharge that
# lands on the exact same day as an infection (or as the day-0 entry)
# never ties in the msm long format. Every discharge is an observed event
# under this scheme (no censoring exists), so this is unconditional. ----
mimic_0 <- mimic_0 %>%
  mutate(
    exit_time = los_time + 0.001
  )

# ---- obs.time.1 / obs.cause.1 : first transition out of state 1.
# No patient can be censored under this scheme (every discharge is
# observed), so there are only two possible outcomes: infected first, or
# discharged without ever being infected. ----
mimic_0 <- mimic_0 %>%
  mutate(
    obs.time.1  = ifelse(hai_status == 1, hai_time, exit_time),
    obs.cause.1 = ifelse(hai_status == 1, "2", "3"),
    obs.time.2  = ifelse(hai_status == 1, exit_time, NA_real_),
    obs.cause.2 = ifelse(hai_status == 1, "3", NA_character_)
  )

cat("--- Cohort summary ---\n")
cat(sprintf("  N patients   : %d\n", nrow(mimic_0)))
cat(sprintf("  HAI (state 2): %d (%.1f%%)\n",
            sum(mimic_0$obs.cause.1 == "2"),
            mean(mimic_0$obs.cause.1 == "2") * 100))
cat(sprintf("  Discharged, no HAI (state 3): %d (%.1f%%)\n",
            sum(mimic_0$obs.cause.1 == "3"),
            mean(mimic_0$obs.cause.1 == "3") * 100))
# NOTE: no "admin censored" category is printed here -- no patient can be
# censored under this scheme, so obs.cause.1 only ever takes the values
# "2" (infected) or "3" (discharged without HAI); the two counts above
# always sum to nrow(mimic_0).
cat("============================================================\n\n")

# ---- Data-driven upper bound for the PPS/CDC survey-date sampling
# window, replacing the removed admin_cens = 60 constant: the maximum
# observed LOS in the cohort. Also reused as the SPP plot horizon
# (see Section on shared helpers, below). ----
pps_date_max <- ceiling(max(mimic_0$los_time, na.rm = TRUE))

cat(sprintf("--- Data-driven PPS survey window: Uniform(%d, %.0f) ---\n\n",
            pps_date_min, pps_date_max))

############################################################
# SECTION 2: REFERENCE ANALYSIS (FULL COHORT)
############################################################

cat("============================================================\n")
cat("SECTION 2 — Reference Analysis (Full cohort)\n")
cat("============================================================\n\n")

qmat <- matrix(
  c(0, 1, 1,
    0, 0, 1,
    0, 0, 0),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("1","2","3"), c("1","2","3"))
)

build_msm_long_ref <- function(data_sub) {

  msm_0 <- data_sub
  msm_0$time     <- entry_time
  msm_0$state    <- 1L
  msm_0$obs_type <- 2L

  # No patient can be censored under this scheme, so obs.cause.1 only ever
  # takes the values "2" (infected) or "3" (discharged without HAI).
  msm_1 <- data_sub
  msm_1$time     <- data_sub$obs.time.1
  msm_1$state    <- ifelse(data_sub$obs.cause.1 == "2", 2L, 3L)
  msm_1$obs_type <- 2L

  # Discharge after HAI is always fully observed under this scheme, so
  # obs.cause.2 only ever takes the value "3" for these (already-infected)
  # patients.
  msm_2 <- data_sub %>% filter(obs.cause.1 == "2")
  msm_2$time     <- msm_2$obs.time.2
  msm_2$state    <- 3L
  msm_2$obs_type <- 2L

  bind_rows(msm_0, msm_1, msm_2) %>%
    arrange(id, time) %>%
    group_by(id) %>%
    mutate(time = time + (row_number() - 1) * 1e-6) %>%
    ungroup()
}

msm_long_ref <- build_msm_long_ref(mimic_0)

# ---- Fit reference MSM (no covariates) --------------------------------
tra_init_ref <- crudeinits.msm(
  state ~ time, subject = id,
  data = msm_long_ref, qmatrix = qmat
)

ref_msm_nocov <- msm(
  state ~ time, subject = id,
  data = msm_long_ref, qmatrix = tra_init_ref,
  obstype = obs_type
)

ref_Q    <- qmatrix.msm(ref_msm_nocov, ci = "none")
ref_Q_ci <- qmatrix.msm(ref_msm_nocov, ci = "delta")

cat("--- Reference MSM constant hazard rates ---\n")
cat(sprintf("  lambda12 (HAI)                 : %.6f\n", ref_Q["1","2"]))
cat(sprintf("  lambda13 (discharge no HAI)    : %.6f\n", ref_Q["1","3"]))
cat(sprintf("  lambda23 (discharge post HAI)  : %.6f\n", ref_Q["2","3"]))
cat("\n")

# ---- Fit reference MSM (with covariates: age + sex) -------------------
ref_msm_cov <- msm(
  state ~ time, subject = id,
  data    = msm_long_ref,
  qmatrix = tra_init_ref,
  obstype = obs_type,
  covariates = ~ age + sex
)

ref_hr <- hazard.msm(ref_msm_cov)

cat("--- Reference MSM covariate HRs ---\n")
print(ref_hr)
cat("\n")

# ---- Helper: safe CI extraction ---------------------------------------
safe_q_ci <- function(fit, from, to, bound) {
  tryCatch({
    q_ci <- qmatrix.msm(fit, ci = "delta")
    as.numeric(q_ci[[bound]][from, to])
  }, error = function(e) NA_real_)
}

# ---- Lambda ratio: lambda23 / lambda13 --------------------------------
extract_lam_delta <- function(fit, from, to) {
  q_ci   <- qmatrix.msm(fit, ci = "delta")
  lam    <- as.numeric(q_ci$estimates[from, to])
  lam_L  <- as.numeric(q_ci$L[from, to])
  lam_U  <- as.numeric(q_ci$U[from, to])
  log_se <- (log(lam_U) - log(lam_L)) / (2 * 1.96)
  list(lam = lam, log_se = log_se)
}

lam23_ref <- extract_lam_delta(ref_msm_nocov, "2", "3")
lam13_ref <- extract_lam_delta(ref_msm_nocov, "1", "3")

ratio_ref    <- lam23_ref$lam / lam13_ref$lam
ratio_log_se <- sqrt(lam23_ref$log_se^2 + lam13_ref$log_se^2)
ratio_ref_L  <- exp(log(ratio_ref) - 1.96 * ratio_log_se)
ratio_ref_U  <- exp(log(ratio_ref) + 1.96 * ratio_log_se)

cat("--- Lambda ratio: lambda23 / lambda13 ---\n")
cat(sprintf("  Reference MSM : %.4f  (95%% CI: %.4f - %.4f)\n",
            ratio_ref, ratio_ref_L, ratio_ref_U))
cat("\n")

ref_results <- list(
  msm_nocov  = ref_msm_nocov,
  msm_cov    = ref_msm_cov,
  Q          = ref_Q,
  Q_ci       = ref_Q_ci,
  ratio_est  = ratio_ref,
  ratio_CI_L = ratio_ref_L,
  ratio_CI_U = ratio_ref_U,
  hr_ref     = ref_hr
)

cat("SECTION 2 COMPLETE\n")
cat("============================================================\n\n")

############################################################
# SHARED HELPERS FOR SECTIONS 3 & 4
############################################################

# ---- Exact SPP via 4-state matrix exponential --------------------------
# t_since_entry is measured from the day-`entry_time` landmark (t = 0 at
# landmark). Since entry_time = 0 under this simplified scheme, this is
# equivalent to time since admission.
compute_spp_exact <- function(h12, h13, h23, t_since_entry_grid) {

  spp_at_t <- function(t) {
    Q4 <- matrix(
      c(-(h12 + h13),  h12,   h13,   0,
              0,      -h23,    0,   h23,
              0,        0,     0,    0,
              0,        0,     0,    0),
      nrow = 4, byrow = TRUE
    )
    P <- expm(Q4 * t)
    c(P[1,1], P[1,2], P[1,3], P[1,4])
  }

  probs <- t(sapply(t_since_entry_grid, spp_at_t))

  data.frame(
    time                = t_since_entry_grid + entry_time,
    Admitted            = probs[, 1],
    HAI                 = probs[, 2],
    Discharged_no_HAI   = probs[, 3],
    Discharged_post_HAI = probs[, 4]
  ) %>%
    pivot_longer(-time, names_to = "State", values_to = "Probability") %>%
    mutate(State = factor(State, levels = c(
      "Discharged_no_HAI", "Admitted",
      "Discharged_post_HAI", "HAI"
    )))
}

state_colors_spp <- c(
  "Admitted"             = "#DB2F73",
  "HAI"                  = "#FFF2A6",
  "Discharged_no_HAI"    = "#7DA9E3",
  "Discharged_post_HAI"  = "#708238"
)
state_labels_spp <- c(
  "Discharged_no_HAI"   = "Discharged w/o HAI",
  "Admitted"            = "Admitted (uninfected)",
  "Discharged_post_HAI" = "Discharged after HAI",
  "HAI"                 = "HAI"
)

plot_spp <- function(combined_df, title) {
  ggplot(combined_df, aes(x = time, y = Probability, fill = State)) +
    geom_area(alpha = 0.85, color = "white", linewidth = 0.2) +
    facet_wrap(~ Panel, nrow = 1) +
    scale_fill_manual(
      name   = "State",
      values = state_colors_spp,
      labels = state_labels_spp,
      guide  = guide_legend(nrow = 1, byrow = TRUE, title.position = "left")
    ) +
    scale_y_continuous(
      expand = c(0, 0), limits = c(0, 1.001),
      breaks = seq(0, 1, 0.25), labels = sprintf("%.2f", seq(0, 1, 0.25))
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = title, x = "Days Since Hospital Admission",
         y = "Probability of Being in State") +
    theme_bw(base_size = 12) +
    theme(
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(face = "bold", size = 11, colour = "#2E6DA4"),
      strip.background  = element_rect(fill = "grey95"),
      axis.text         = element_text(colour = "black", size = 10),
      axis.title        = element_text(colour = "black", size = 12),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.title      = element_text(size = 11),
      legend.text       = element_text(size = 10),
      legend.key.size   = unit(0.7, "lines"),
      plot.title        = element_text(hjust = 0.5, face = "bold", size = 13)
    )
}

transition_labeller <- as_labeller(
  c(
    "12" = "lambda[12]~(Admission~to~HAI)",
    "13" = "lambda[13]~(Admission~to~Discharge)",
    "23" = "lambda[23]~(HAI~to~Discharge)"
  ),
  default = label_parsed
)

theme_boxplot <- function() {
  theme_bw(base_size = 12) +
    theme(
      strip.text         = element_text(face = "bold", size = 11, colour = "#2E6DA4"),
      strip.background   = element_rect(fill = "grey95"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text.x        = element_text(face = "bold", size = 10),
      axis.text.y        = element_text(size = 10),
      axis.title.y       = element_text(size = 11),
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle      = element_text(hjust = 0.5, colour = "grey40", size = 10),
      legend.position    = "bottom",
      legend.text        = element_text(size = 10),
      legend.title       = element_text(size = 10, face = "bold"),
      legend.key.width   = unit(1.8, "lines")
    )
}

fill_colors_pps <- c("PPS Unweighted" = "brown2", "PPS Weighted" = "darkolivegreen3")

# ---- CDC has 4 methods. Colors per spec:
#   Unweighted + IC-naive     -> red
#   Unweighted + IC-corrected -> orange
#   Weighted   + IC-naive     -> yellow
#   Weighted   + IC-corrected -> green
fill_colors_cdc <- c(
  "Unweighted + IC-naive"     = "red",
  "Unweighted + IC-corrected" = "orange",
  "Weighted + IC-naive"       = "yellow",
  "Weighted + IC-corrected"   = "green3"
)

# ---- Fig A: hazard rate boxplot maker ---------------------------------
make_hz_boxplot <- function(runs_long, ref_df, fill_colors, protocol_label, subtitle,
                             x_breaks = c(1, 2), x_labels = c("Unweighted", "Weighted"),
                             x_limits = c(0.4, 2.6)) {
  ggplot(
    runs_long,
    aes(x = x_pos, y = hazard, fill = weighting,
        group = interaction(x_pos, weighting))
  ) +
    geom_rect(
      data = ref_df,
      aes(ymin = ref_ci_L, ymax = ref_ci_U, xmin = -Inf, xmax = Inf),
      fill = "grey80", alpha = 0.55, inherit.aes = FALSE
    ) +
    geom_hline(
      data = ref_df,
      aes(yintercept = ref_est, colour = "Reference Analysis",
          linetype = "Reference Analysis"),
      linewidth = 0.9, inherit.aes = FALSE
    ) +
    geom_boxplot(outlier.size = 1.0, width = 0.55, alpha = 0.90) +
    facet_wrap(~ transition, nrow = 1, scales = "free_y",
               labeller = transition_labeller) +
    scale_fill_manual(name = NULL, values = fill_colors) +
    scale_colour_manual(name = NULL, values = c("Reference Analysis" = "black")) +
    scale_linetype_manual(name = NULL, values = c("Reference Analysis" = "dotted")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = x_limits) +
    labs(title = sprintf("%s - Constant Hazard Rates", protocol_label),
         subtitle = subtitle, x = NULL, y = "Estimated hazard rate") +
    theme_boxplot() +
    guides(
      fill     = guide_legend(order = 1, nrow = 1, override.aes = list(linetype = 0)),
      colour   = guide_legend(order = 2, nrow = 1,
                              override.aes = list(linetype = "dotted", linewidth = 0.9, shape = NA)),
      linetype = "none"
    )
}

# ---- Fig C: HR covariate boxplot maker --------------------------------
hr_labeller_cov <- as_labeller(
  c("age" = "Age~(per~year)", "sex" = "Sex~(female~vs~male)"),
  default = label_parsed
)

make_hr_boxplot <- function(hr_long, ref_hr_df, fill_colors, protocol_label, subtitle,
                             x_breaks = c(1, 2), x_labels = c("Unweighted", "Weighted"),
                             x_limits = c(0.4, 2.6)) {
  ggplot(
    hr_long,
    aes(x = x_pos, y = hr_val, fill = weighting,
        group = interaction(x_pos, weighting))
  ) +
    geom_rect(
      data = ref_hr_df,
      aes(ymin = ref_hr_L, ymax = ref_hr_U, xmin = -Inf, xmax = Inf),
      fill = "grey80", alpha = 0.55, inherit.aes = FALSE
    ) +
    geom_hline(
      data = ref_hr_df,
      aes(yintercept = ref_hr_est, colour = "Reference Analysis",
          linetype = "Reference Analysis"),
      linewidth = 0.9, inherit.aes = FALSE
    ) +
    geom_boxplot(outlier.size = 1.0, width = 0.55, alpha = 0.90) +
    facet_grid(covariate ~ transition, scales = "free_y",
               labeller = labeller(covariate = hr_labeller_cov,
                                   transition = transition_labeller)) +
    scale_fill_manual(name = NULL, values = fill_colors) +
    scale_colour_manual(name = NULL, values = c("Reference Analysis" = "black")) +
    scale_linetype_manual(name = NULL, values = c("Reference Analysis" = "dotted")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = x_limits) +
    scale_y_log10() +
    labs(title = sprintf("%s - Covariate Hazard Ratios", protocol_label),
         subtitle = subtitle, x = NULL, y = "Hazard Ratio (log scale)") +
    theme_boxplot() +
    guides(
      fill     = guide_legend(order = 1, nrow = 1, override.aes = list(linetype = 0)),
      colour   = guide_legend(order = 2, nrow = 1,
                              override.aes = list(linetype = "dotted", linewidth = 0.9, shape = NA)),
      linetype = "none"
    )
}

# ---- Fig D: ratio lambda23/lambda13 boxplot maker ---------------------
make_ratio_boxplot <- function(ratio_long, ref_ratio_est, ref_ratio_L, ref_ratio_U,
                               fill_colors, protocol_label, subtitle,
                               x_breaks = c(1, 2), x_labels = c("Unweighted", "Weighted"),
                               x_limits = c(0.4, 2.6)) {
  ref_df_r <- data.frame(
    dummy = "ratio", ref_est = ref_ratio_est,
    ref_ci_L = ref_ratio_L, ref_ci_U = ref_ratio_U
  )

  ggplot(
    ratio_long,
    aes(x = x_pos, y = ratio_val, fill = weighting,
        group = interaction(x_pos, weighting))
  ) +
    geom_rect(
      data = ref_df_r,
      aes(ymin = ref_ci_L, ymax = ref_ci_U, xmin = -Inf, xmax = Inf),
      fill = "grey80", alpha = 0.55, inherit.aes = FALSE
    ) +
    geom_hline(
      data = ref_df_r,
      aes(yintercept = ref_est, colour = "Reference Analysis",
          linetype = "Reference Analysis"),
      linewidth = 0.9, inherit.aes = FALSE
    ) +
    geom_boxplot(outlier.size = 1.0, width = 0.55, alpha = 0.90) +
    scale_fill_manual(name = NULL, values = fill_colors) +
    scale_colour_manual(name = NULL, values = c("Reference Analysis" = "black")) +
    scale_linetype_manual(name = NULL, values = c("Reference Analysis" = "dotted")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = x_limits) +
    labs(title = sprintf("%s - Ratio \u03bb23 / \u03bb13", protocol_label),
         subtitle = subtitle, x = NULL, y = expression(lambda[23] / lambda[13])) +
    theme_boxplot() +
    guides(
      fill     = guide_legend(order = 1, nrow = 1, override.aes = list(linetype = 0)),
      colour   = guide_legend(order = 2, nrow = 1,
                              override.aes = list(linetype = "dotted", linewidth = 0.9, shape = NA)),
      linetype = "none"
    )
}

# ---- Reference data frames for Fig A -----------------------------------
ref_hz_df <- data.frame(
  transition = factor(c("12", "13", "23"), levels = c("12","13","23")),
  ref_est    = c(as.numeric(ref_Q["1","2"]), as.numeric(ref_Q["1","3"]), as.numeric(ref_Q["2","3"])),
  ref_ci_L   = c(safe_q_ci(ref_msm_nocov, "1", "2", "L"),
                 safe_q_ci(ref_msm_nocov, "1", "3", "L"),
                 safe_q_ci(ref_msm_nocov, "2", "3", "L")),
  ref_ci_U   = c(safe_q_ci(ref_msm_nocov, "1", "2", "U"),
                 safe_q_ci(ref_msm_nocov, "1", "3", "U"),
                 safe_q_ci(ref_msm_nocov, "2", "3", "U"))
)

# ---- Reference data frame for Fig C (HR covariates) --------------------
extract_ref_hr_df <- function(hr_list) {
  bind_rows(lapply(c("age", "sex"), function(cov) {
    bind_rows(lapply(c("12", "13", "23"), function(tr) {

      tr_key <- switch(tr,
        "12" = "State 1 - State 2",
        "13" = "State 1 - State 3",
        "23" = "State 2 - State 3"
      )

      hr_mat <- tryCatch(hr_list[[cov]], error = function(e) NULL)
      if (is.null(hr_mat)) {
        return(data.frame(covariate = cov, transition = tr,
                          ref_hr_est = NA, ref_hr_L = NA, ref_hr_U = NA))
      }

      rn  <- rownames(hr_mat)
      idx <- grep(paste0(gsub(" ", ".", tr_key), "|",
                         gsub("-", ".", tr_key), "|",
                         paste0(strsplit(tr,"")[[1]], collapse = ".*")),
                  rn)
      if (length(idx) == 0) idx <- seq_len(nrow(hr_mat))[
        grepl(sub("State ","",strsplit(tr_key," - ")[[1]][1]), rn) &
        grepl(sub("State ","",strsplit(tr_key," - ")[[1]][2]), rn)
      ]

      if (length(idx) == 0) {
        return(data.frame(covariate = cov, transition = tr,
                          ref_hr_est = NA, ref_hr_L = NA, ref_hr_U = NA))
      }

      data.frame(
        covariate  = cov, transition = tr,
        ref_hr_est = hr_mat[idx[1], 1],
        ref_hr_L   = hr_mat[idx[1], 2],
        ref_hr_U   = hr_mat[idx[1], 3]
      )
    }))
  })) %>%
    mutate(
      covariate  = factor(covariate,  levels = c("age", "sex")),
      transition = factor(transition, levels = c("12","13","23"))
    )
}

ref_hr_df <- extract_ref_hr_df(ref_hr)

# ---- Shared safe extractors --------------------------------------------
safe_q <- function(Q, from, to) {
  tryCatch(as.numeric(Q[from, to]), error = function(e) NA_real_)
}

extract_hr_val <- function(hr_list, cov, tr_from, tr_to) {
  if (is.null(hr_list)) return(NA_real_)
  tryCatch({
    mat <- hr_list[[cov]]
    rn  <- rownames(mat)
    idx <- grep(paste0(tr_from, ".*", tr_to), rn)
    if (length(idx) == 0) return(NA_real_)
    as.numeric(mat[idx[1], 1])
  }, error = function(e) NA_real_)
}

# SIMPLIFIED: grid runs from entry_time (0) to the data-driven horizon
# pps_date_max (max observed LOS), replacing the old 0:(admin_cens - entry_time).
time_grid_since_entry <- 0:(pps_date_max - entry_time)
n_cores <- max(1L, detectCores() - 1L)

############################################################
# SECTION 3: PPS PROTOCOL — 100 RUNS
############################################################

cat("============================================================\n")
cat("SECTION 3 — PPS 100 runs\n")
cat("============================================================\n\n")

run_one_pps <- function(data, run_id, entry_time, pps_date_min, pps_date_max, qmat) {

  data$PPS_date <- runif(nrow(data), min = pps_date_min, max = pps_date_max)  # SIMPLIFIED: data-driven max
  PPS_0 <- data[data$los_time > data$PPS_date, ]

  PPS_0$prevalent_infection <- ifelse(
    PPS_0$hai_status == 1 & PPS_0$hai_time < PPS_0$PPS_date, 1L, 0L
  )

  # Full retrospective panel for the subsampled cohort (same structure as
  # the full-cohort reference panel — only the subject subset differs)
  msm_0 <- PPS_0
  msm_0$time     <- entry_time
  msm_0$state    <- 1L
  msm_0$obs_type <- 2L

  # No patient can be censored under this scheme (see build_msm_long_ref()).
  msm_1 <- PPS_0
  msm_1$time     <- PPS_0$obs.time.1
  msm_1$state    <- ifelse(PPS_0$obs.cause.1 == "2", 2L, 3L)
  msm_1$obs_type <- 2L

  msm_2 <- PPS_0 %>% filter(obs.cause.1 == "2")
  msm_2$time     <- msm_2$obs.time.2
  msm_2$state    <- 3L
  msm_2$obs_type <- 2L

  msm_long <- bind_rows(msm_0, msm_1, msm_2) %>%
    arrange(id, time) %>%
    group_by(id) %>%
    mutate(time = time + (row_number() - 1) * 1e-6) %>%
    ungroup() %>%
    mutate(weight = 1 / los_time)

  tra <- crudeinits.msm(state ~ time, subject = id, data = msm_long, qmatrix = qmat)

  fit_uw <- msm(state ~ time, subject = id, data = msm_long,
                qmatrix = tra, obstype = obs_type)
  fit_w  <- msm(state ~ time, subject = id, data = msm_long,
                qmatrix = tra, obstype = obs_type, subject.weights = weight)

  fit_uw_cov <- tryCatch(
    msm(state ~ time, subject = id, data = msm_long, qmatrix = tra,
        obstype = obs_type, covariates = ~ age + sex),
    error = function(e) NULL
  )
  fit_w_cov <- tryCatch(
    msm(state ~ time, subject = id, data = msm_long, qmatrix = tra,
        obstype = obs_type, subject.weights = weight, covariates = ~ age + sex),
    error = function(e) NULL
  )

  get_hr <- function(fit) {
    if (is.null(fit)) return(NULL)
    tryCatch(hazard.msm(fit), error = function(e) NULL)
  }

  list(
    run   = run_id,
    Q_uw  = qmatrix.msm(fit_uw, ci = "none"),
    Q_w   = qmatrix.msm(fit_w,  ci = "none"),
    hr_uw = get_hr(fit_uw_cov),
    hr_w  = get_hr(fit_w_cov),
    error = NA_character_
  )
}

cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl,
  varlist = c("mimic_0", "entry_time", "pps_date_min", "pps_date_max", "qmat", "run_one_pps"),
  envir = environment()
)
clusterEvalQ(cl, { library(msm); library(dplyr); library(tibble) })

PPS_results <- foreach(
  i = seq_len(n_runs),
  .packages      = c("msm", "dplyr", "tibble"),
  .errorhandling = "pass"
) %dopar% {
  set.seed(123 + i)
  tryCatch(
    run_one_pps(mimic_0, i, entry_time, pps_date_min, pps_date_max, qmat),
    error = function(e) list(run = i, Q_uw = NULL, Q_w = NULL,
                             hr_uw = NULL, hr_w = NULL,
                             error = conditionMessage(e))
  )
}

stopCluster(cl)
registerDoSEQ()

n_err_pps <- sum(sapply(PPS_results, function(x) !is.na(x$error)))
message(sprintf("PPS complete. Failed runs: %d", n_err_pps))

pps_hz <- bind_rows(lapply(PPS_results, function(res) {
  if (!is.na(res$error)) return(NULL)
  data.frame(
    run      = res$run,
    hz_12_uw = safe_q(res$Q_uw, "1","2"), hz_13_uw = safe_q(res$Q_uw, "1","3"), hz_23_uw = safe_q(res$Q_uw, "2","3"),
    hz_12_w  = safe_q(res$Q_w,  "1","2"), hz_13_w  = safe_q(res$Q_w,  "1","3"), hz_23_w  = safe_q(res$Q_w,  "2","3")
  )
}))

pps_hz_long <- pps_hz %>%
  pivot_longer(cols = hz_12_uw:hz_23_w, names_to = c("transition", "wt"),
               names_pattern = "hz_([0-9]+)_(uw|w)", values_to = "hazard") %>%
  mutate(
    weighting  = ifelse(wt == "uw", "PPS Unweighted", "PPS Weighted"),
    weighting  = factor(weighting, levels = c("PPS Unweighted","PPS Weighted")),
    x_pos      = ifelse(wt == "uw", 1, 2),
    transition = factor(transition, levels = c("12","13","23"))
  )

pps_hr_long <- bind_rows(lapply(PPS_results, function(res) {
  if (!is.na(res$error)) return(NULL)
  bind_rows(lapply(c("age","sex"), function(cov) {
    bind_rows(lapply(list(c("1","2","12"), c("1","3","13"), c("2","3","23")), function(tr) {
      data.frame(
        run = res$run, covariate = cov, transition = tr[3],
        hr_uw = extract_hr_val(res$hr_uw, cov, tr[1], tr[2]),
        hr_w  = extract_hr_val(res$hr_w,  cov, tr[1], tr[2])
      )
    }))
  }))
})) %>%
  pivot_longer(c(hr_uw, hr_w), names_to = "wt", values_to = "hr_val") %>%
  mutate(
    weighting  = ifelse(wt == "hr_uw", "PPS Unweighted", "PPS Weighted"),
    weighting  = factor(weighting, levels = c("PPS Unweighted","PPS Weighted")),
    x_pos      = ifelse(wt == "hr_uw", 1, 2),
    covariate  = factor(covariate,  levels = c("age","sex")),
    transition = factor(transition, levels = c("12","13","23"))
  )

pps_ratio_long <- pps_hz %>%
  mutate(ratio_uw = hz_23_uw / hz_13_uw, ratio_w = hz_23_w / hz_13_w) %>%
  pivot_longer(c(ratio_uw, ratio_w), names_to = "wt", values_to = "ratio_val") %>%
  mutate(
    weighting = ifelse(wt == "ratio_uw", "PPS Unweighted", "PPS Weighted"),
    weighting = factor(weighting, levels = c("PPS Unweighted","PPS Weighted")),
    x_pos     = ifelse(wt == "ratio_uw", 1, 2)
  )

mean_pps_hz_uw <- colMeans(pps_hz[, c("hz_12_uw","hz_13_uw","hz_23_uw")], na.rm = TRUE)
mean_pps_hz_w  <- colMeans(pps_hz[, c("hz_12_w","hz_13_w","hz_23_w")],  na.rm = TRUE)

panel_levels_pps <- c("Full Cohort (Reference)", "PPS Unweighted", "PPS Weighted")

spp_pps <- bind_rows(
  compute_spp_exact(ref_Q["1","2"], ref_Q["1","3"], ref_Q["2","3"],
                    time_grid_since_entry) %>% mutate(Panel = "Full Cohort (Reference)"),
  compute_spp_exact(mean_pps_hz_uw["hz_12_uw"], mean_pps_hz_uw["hz_13_uw"], mean_pps_hz_uw["hz_23_uw"],
                    time_grid_since_entry) %>% mutate(Panel = "PPS Unweighted"),
  compute_spp_exact(mean_pps_hz_w["hz_12_w"], mean_pps_hz_w["hz_13_w"], mean_pps_hz_w["hz_23_w"],
                    time_grid_since_entry) %>% mutate(Panel = "PPS Weighted")
) %>% mutate(Panel = factor(Panel, levels = panel_levels_pps))

fig_pps_A <- make_hz_boxplot(pps_hz_long, ref_hz_df, fill_colors_pps, "PPS",
                              "100 PPS runs — MIMIC cohort (simplified: no LOS filter, no admin censoring)")
print(fig_pps_A)

fig_pps_B <- plot_spp(spp_pps, "PPS — Stacked Probability Plot (mean hazards across 100 runs)")
print(fig_pps_B)

fig_pps_C <- make_hr_boxplot(pps_hr_long, ref_hr_df, fill_colors_pps, "PPS",
                              "100 PPS runs — MIMIC cohort (simplified)")
print(fig_pps_C)

fig_pps_D <- make_ratio_boxplot(pps_ratio_long, ref_results$ratio_est,
                                ref_results$ratio_CI_L, ref_results$ratio_CI_U,
                                fill_colors_pps, "PPS", "100 PPS runs — MIMIC cohort (simplified)")
print(fig_pps_D)

cat("SECTION 3 COMPLETE\n")
cat("============================================================\n\n")

############################################################
# SECTION 4: CDC PROTOCOL — 100 RUNS
#
# FOUR methods:
#   1. Unweighted + IC-naive
#   2. Unweighted + IC-corrected
#   3. Weighted   + IC-naive
#   4. Weighted   + IC-corrected
#
# SIMPLIFIED: no administrative-censoring override (no state 99,
# no censor = 99) -- every patient's discharge is fully observed, so
# that branch never applied and has been removed.
############################################################

cat("============================================================\n")
cat("SECTION 4 — CDC 100 runs\n")
cat("============================================================\n\n")

build_cdc_panel <- function(PPS_sub, entry_time) {

  # Non-prevalent patients get a "known uninfected as of PPS_date" row.
  # PPS_date is sampled independently (Uniform(pps_date_min, pps_date_max))
  # and can in principle fall at or below entry_time; since entry_time is
  # now 0 and PPS_date is drawn from a continuous distribution, this is a
  # near-vacuous edge case, but the clamp is retained defensively: there is
  # no information about a patient's state before entry_time in any case,
  # so a PPS_date <= entry_time carries no more information than
  # "uninfected at entry_time".
  #
  # NOTE: no equivalent clamp is needed on the infection-time branch;
  # hai_time is, by construction, always >= entry_time.
  #
  # SIMPLIFIED: the administrative-censoring override (state 99, obstype 1,
  # for patients not yet discharged/infected by a fixed cutoff) has been
  # removed. Every patient's discharge is fully observed under this scheme,
  # so that branch never applied.

  df <- data.frame(
    id                     = PPS_sub$id,
    age                    = PPS_sub$age,
    sex                    = PPS_sub$sex,
    los_time               = PPS_sub$los_time,
    exit_time              = PPS_sub$exit_time,
    prevalent_infection    = PPS_sub$prevalent_infection,
    PPS_date               = PPS_sub$PPS_date,
    cens_time_to_infection = ifelse(PPS_sub$prevalent_infection == 1, PPS_sub$hai_time, NA_real_)
  )

  df$time     <- entry_time
  df$state    <- 1L
  df$obs_type <- 2L

  df1          <- df
  df1$time     <- ifelse(df1$prevalent_infection == 1,
                          df1$cens_time_to_infection,
                          ifelse(df1$PPS_date > entry_time,
                                 df1$PPS_date, entry_time + 1e-4))
  df1$state    <- ifelse(df1$prevalent_infection == 1, 2L, 1L)
  df1$obs_type <- 2L

  df2          <- df
  df2$time     <- df2$exit_time
  df2$state    <- 3L
  df2$obs_type <- ifelse(df2$prevalent_infection == 1, 2L, 3L)

  cdc_msm <- rbind(df, df1, df2) %>% arrange(id)

  cdc_msm$weight <- 1 / cdc_msm$los_time
  cdc_msm
}

build_cdc_panel_no_ic <- function(PPS_sub, entry_time) {

  # Same underlying survey-limited data as build_cdc_panel() -- prevalent_infection
  # marks what the PPS survey would actually see. The difference from
  # build_cdc_panel() is purely in how that same ambiguous observation gets
  # coded for msm: here it is naively carried forward as if it were exact
  # (obs_type = 2 throughout), instead of being flagged as uncertain
  # (obs_type = 3 at discharge in the IC-corrected panel).
  #
  # This same panel is reused for BOTH "Weighted + IC-naive" (subject.weights
  # supplied at the msm() fitting stage) and "Unweighted + IC-naive"
  # (subject.weights omitted).
  #
  # No patient can be censored under this scheme -- every discharge is
  # observed, so the discharge row below is unconditionally coded as an
  # exact transition to state 3 at exit_time. los_status is not needed
  # here at all (it is always 1) and has been dropped from this panel.

  df <- data.frame(
    id                     = PPS_sub$id,
    age                    = PPS_sub$age,
    sex                    = PPS_sub$sex,
    los_time               = PPS_sub$los_time,
    exit_time              = PPS_sub$exit_time,
    prevalent_infection    = PPS_sub$prevalent_infection,
    PPS_date               = PPS_sub$PPS_date,
    cens_time_to_infection = ifelse(PPS_sub$prevalent_infection == 1, PPS_sub$hai_time, NA_real_)
  )

  df$time     <- entry_time
  df$state    <- 1L
  df$obs_type <- 2L

  df1          <- df
  df1$time     <- ifelse(df1$prevalent_infection == 1,
                          df1$cens_time_to_infection,
                          df1$PPS_date)
  df1$state    <- ifelse(df1$prevalent_infection == 1, 2L, 1L)
  df1$obs_type <- 2L

  # Same defensive clamp as build_cdc_panel(): PPS_date is sampled
  # independently and could in principle fall <= entry_time.
  df1$time <- ifelse(df1$prevalent_infection == 0 & df1$time <= entry_time,
                      entry_time + 1e-4,
                      df1$time)

  # Discharge row: naive coding -- the apparent PPS-observed state is
  # carried forward as if exact, no obs_type = 3.
  df2          <- df
  df2$time     <- df2$exit_time
  df2$state    <- 3L
  df2$obs_type <- 2L

  cdc_msm_noic <- rbind(df, df1, df2) %>% arrange(id)

  cdc_msm_noic$weight <- 1 / cdc_msm_noic$los_time
  cdc_msm_noic
}

run_one_cdc <- function(data, run_id, entry_time, pps_date_min, pps_date_max, qmat) {

  data$PPS_date <- runif(nrow(data), min = pps_date_min, max = pps_date_max)  # SIMPLIFIED: data-driven max
  PPS_0 <- data[data$los_time > data$PPS_date, ]

  PPS_0$prevalent_infection <- ifelse(
    PPS_0$hai_status == 1 & PPS_0$hai_time < PPS_0$PPS_date, 1L, 0L
  )

  cdc_msm_ic   <- build_cdc_panel(PPS_0, entry_time)
  cdc_msm_noic <- build_cdc_panel_no_ic(PPS_0, entry_time)

  # SIMPLIFIED: no censor = 99 argument -- the IC-corrected panel no longer
  # contains a state-99 placeholder (see build_cdc_panel()).
  tra_ic   <- crudeinits.msm(state ~ time, subject = id, data = cdc_msm_ic,   qmatrix = qmat)
  tra_noic <- crudeinits.msm(state ~ time, subject = id, data = cdc_msm_noic, qmatrix = qmat)

  # Unweighted + IC-corrected
  fit_uw_ic  <- msm(state ~ time, subject = id, data = cdc_msm_ic,
                     qmatrix = tra_ic, obstype = obs_type)
  # Weighted + IC-corrected
  fit_w_ic   <- msm(state ~ time, subject = id, data = cdc_msm_ic,
                     qmatrix = tra_ic, obstype = obs_type, subject.weights = weight)
  # Weighted + IC-naive
  fit_w_noic <- msm(state ~ time, subject = id, data = cdc_msm_noic,
                     qmatrix = tra_noic, obstype = obs_type, subject.weights = weight)
  # Unweighted + IC-naive -- same panel as Weighted + IC-naive
  # (cdc_msm_noic / tra_noic), but subject.weights omitted.
  fit_uw_noic <- msm(state ~ time, subject = id, data = cdc_msm_noic,
                      qmatrix = tra_noic, obstype = obs_type)

  fit_uw_ic_cov <- tryCatch(
    msm(state ~ time, subject = id, data = cdc_msm_ic, qmatrix = tra_ic, obstype = obs_type,
        covariates = ~ age + sex),
    error = function(e) NULL
  )
  fit_w_ic_cov <- tryCatch(
    msm(state ~ time, subject = id, data = cdc_msm_ic, qmatrix = tra_ic, obstype = obs_type,
        subject.weights = weight, covariates = ~ age + sex),
    error = function(e) NULL
  )
  fit_w_noic_cov <- tryCatch(
    msm(state ~ time, subject = id, data = cdc_msm_noic, qmatrix = tra_noic, obstype = obs_type,
        subject.weights = weight, covariates = ~ age + sex),
    error = function(e) NULL
  )
  fit_uw_noic_cov <- tryCatch(
    msm(state ~ time, subject = id, data = cdc_msm_noic, qmatrix = tra_noic, obstype = obs_type,
        covariates = ~ age + sex),
    error = function(e) NULL
  )

  get_hr <- function(fit) {
    if (is.null(fit)) return(NULL)
    tryCatch(hazard.msm(fit), error = function(e) NULL)
  }

  list(
    run        = run_id,
    Q_uw_ic    = qmatrix.msm(fit_uw_ic,   ci = "none"),
    Q_w_ic     = qmatrix.msm(fit_w_ic,    ci = "none"),
    Q_w_noic   = qmatrix.msm(fit_w_noic,  ci = "none"),
    Q_uw_noic  = qmatrix.msm(fit_uw_noic, ci = "none"),
    hr_uw_ic   = get_hr(fit_uw_ic_cov),
    hr_w_ic    = get_hr(fit_w_ic_cov),
    hr_w_noic  = get_hr(fit_w_noic_cov),
    hr_uw_noic = get_hr(fit_uw_noic_cov),
    error      = NA_character_
  )
}

cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl,
  varlist = c("mimic_0", "qmat", "entry_time", "pps_date_min", "pps_date_max", "run_one_cdc",
              "build_cdc_panel", "build_cdc_panel_no_ic"),
  envir = environment()
)
clusterEvalQ(cl, { library(msm); library(dplyr); library(tibble) })

CDC_results <- foreach(
  i = seq_len(n_runs),
  .packages      = c("msm", "dplyr", "tibble"),
  .errorhandling = "pass"
) %dopar% {
  set.seed(200 + i)
  tryCatch(
    run_one_cdc(mimic_0, i, entry_time, pps_date_min, pps_date_max, qmat),
    error = function(e) list(run = i, Q_uw_ic = NULL, Q_w_ic = NULL, Q_w_noic = NULL, Q_uw_noic = NULL,
                             hr_uw_ic = NULL, hr_w_ic = NULL, hr_w_noic = NULL, hr_uw_noic = NULL,
                             error = conditionMessage(e))
  )
}

stopCluster(cl)
registerDoSEQ()

n_err_cdc <- sum(sapply(CDC_results, function(x) !is.na(x$error)))
message(sprintf("CDC complete. Failed runs: %d", n_err_cdc))

cdc_hz <- bind_rows(lapply(CDC_results, function(res) {
  if (!is.na(res$error)) return(NULL)
  data.frame(
    run            = res$run,
    hz_12_uw_noic  = safe_q(res$Q_uw_noic, "1","2"), hz_13_uw_noic  = safe_q(res$Q_uw_noic, "1","3"), hz_23_uw_noic  = safe_q(res$Q_uw_noic, "2","3"),
    hz_12_uw_ic    = safe_q(res$Q_uw_ic,   "1","2"), hz_13_uw_ic    = safe_q(res$Q_uw_ic,   "1","3"), hz_23_uw_ic    = safe_q(res$Q_uw_ic,   "2","3"),
    hz_12_w_noic   = safe_q(res$Q_w_noic,  "1","2"), hz_13_w_noic   = safe_q(res$Q_w_noic,  "1","3"), hz_23_w_noic   = safe_q(res$Q_w_noic,  "2","3"),
    hz_12_w_ic     = safe_q(res$Q_w_ic,    "1","2"), hz_13_w_ic     = safe_q(res$Q_w_ic,    "1","3"), hz_23_w_ic     = safe_q(res$Q_w_ic,    "2","3")
  )
}))

# Left-to-right plot order per requested colors:
#   Unweighted+IC-naive (red), Unweighted+IC-corrected (orange),
#   Weighted+IC-naive (yellow), Weighted+IC-corrected (green)
cdc_method_levels <- c("Unweighted + IC-naive", "Unweighted + IC-corrected",
                        "Weighted + IC-naive", "Weighted + IC-corrected")

cdc_hz_long <- cdc_hz %>%
  pivot_longer(cols = hz_12_uw_noic:hz_23_w_ic, names_to = c("transition", "wt"),
               names_pattern = "hz_([0-9]+)_(uw_noic|uw_ic|w_noic|w_ic)", values_to = "hazard") %>%
  mutate(
    weighting  = case_when(
      wt == "uw_noic" ~ "Unweighted + IC-naive",
      wt == "uw_ic"   ~ "Unweighted + IC-corrected",
      wt == "w_noic"  ~ "Weighted + IC-naive",
      wt == "w_ic"    ~ "Weighted + IC-corrected"
    ),
    weighting  = factor(weighting, levels = cdc_method_levels),
    x_pos      = case_when(wt == "uw_noic" ~ 1, wt == "uw_ic" ~ 2, wt == "w_noic" ~ 3, wt == "w_ic" ~ 4),
    transition = factor(transition, levels = c("12","13","23"))
  )

cdc_hr_long <- bind_rows(lapply(CDC_results, function(res) {
  if (!is.na(res$error)) return(NULL)
  bind_rows(lapply(c("age","sex"), function(cov) {
    bind_rows(lapply(list(c("1","2","12"), c("1","3","13"), c("2","3","23")), function(tr) {
      data.frame(
        run = res$run, covariate = cov, transition = tr[3],
        hr_uw_noic = extract_hr_val(res$hr_uw_noic, cov, tr[1], tr[2]),
        hr_uw_ic   = extract_hr_val(res$hr_uw_ic,   cov, tr[1], tr[2]),
        hr_w_noic  = extract_hr_val(res$hr_w_noic,  cov, tr[1], tr[2]),
        hr_w_ic    = extract_hr_val(res$hr_w_ic,    cov, tr[1], tr[2])
      )
    }))
  }))
})) %>%
  pivot_longer(c(hr_uw_noic, hr_uw_ic, hr_w_noic, hr_w_ic), names_to = "wt", values_to = "hr_val") %>%
  mutate(
    weighting  = case_when(
      wt == "hr_uw_noic" ~ "Unweighted + IC-naive",
      wt == "hr_uw_ic"   ~ "Unweighted + IC-corrected",
      wt == "hr_w_noic"  ~ "Weighted + IC-naive",
      wt == "hr_w_ic"    ~ "Weighted + IC-corrected"
    ),
    weighting  = factor(weighting, levels = cdc_method_levels),
    x_pos      = case_when(wt == "hr_uw_noic" ~ 1, wt == "hr_uw_ic" ~ 2, wt == "hr_w_noic" ~ 3, wt == "hr_w_ic" ~ 4),
    covariate  = factor(covariate,  levels = c("age","sex")),
    transition = factor(transition, levels = c("12","13","23"))
  )

cdc_ratio_long <- cdc_hz %>%
  mutate(
    ratio_uw_noic = hz_23_uw_noic / hz_13_uw_noic,
    ratio_uw_ic   = hz_23_uw_ic   / hz_13_uw_ic,
    ratio_w_noic  = hz_23_w_noic  / hz_13_w_noic,
    ratio_w_ic    = hz_23_w_ic    / hz_13_w_ic
  ) %>%
  pivot_longer(c(ratio_uw_noic, ratio_uw_ic, ratio_w_noic, ratio_w_ic), names_to = "wt", values_to = "ratio_val") %>%
  mutate(
    weighting = case_when(
      wt == "ratio_uw_noic" ~ "Unweighted + IC-naive",
      wt == "ratio_uw_ic"   ~ "Unweighted + IC-corrected",
      wt == "ratio_w_noic"  ~ "Weighted + IC-naive",
      wt == "ratio_w_ic"    ~ "Weighted + IC-corrected"
    ),
    weighting = factor(weighting, levels = cdc_method_levels),
    x_pos     = case_when(wt == "ratio_uw_noic" ~ 1, wt == "ratio_uw_ic" ~ 2, wt == "ratio_w_noic" ~ 3, wt == "ratio_w_ic" ~ 4)
  )

mean_cdc_hz_uw_noic <- colMeans(cdc_hz[, c("hz_12_uw_noic","hz_13_uw_noic","hz_23_uw_noic")], na.rm = TRUE)
mean_cdc_hz_uw_ic   <- colMeans(cdc_hz[, c("hz_12_uw_ic","hz_13_uw_ic","hz_23_uw_ic")],       na.rm = TRUE)
mean_cdc_hz_w_noic  <- colMeans(cdc_hz[, c("hz_12_w_noic","hz_13_w_noic","hz_23_w_noic")],    na.rm = TRUE)
mean_cdc_hz_w_ic    <- colMeans(cdc_hz[, c("hz_12_w_ic","hz_13_w_ic","hz_23_w_ic")],          na.rm = TRUE)

panel_levels_cdc <- c("Full Cohort (Reference)", cdc_method_levels)

spp_cdc <- bind_rows(
  compute_spp_exact(ref_Q["1","2"], ref_Q["1","3"], ref_Q["2","3"],
                    time_grid_since_entry) %>% mutate(Panel = "Full Cohort (Reference)"),
  compute_spp_exact(mean_cdc_hz_uw_noic["hz_12_uw_noic"], mean_cdc_hz_uw_noic["hz_13_uw_noic"], mean_cdc_hz_uw_noic["hz_23_uw_noic"],
                    time_grid_since_entry) %>% mutate(Panel = "Unweighted + IC-naive"),
  compute_spp_exact(mean_cdc_hz_uw_ic["hz_12_uw_ic"], mean_cdc_hz_uw_ic["hz_13_uw_ic"], mean_cdc_hz_uw_ic["hz_23_uw_ic"],
                    time_grid_since_entry) %>% mutate(Panel = "Unweighted + IC-corrected"),
  compute_spp_exact(mean_cdc_hz_w_noic["hz_12_w_noic"], mean_cdc_hz_w_noic["hz_13_w_noic"], mean_cdc_hz_w_noic["hz_23_w_noic"],
                    time_grid_since_entry) %>% mutate(Panel = "Weighted + IC-naive"),
  compute_spp_exact(mean_cdc_hz_w_ic["hz_12_w_ic"], mean_cdc_hz_w_ic["hz_13_w_ic"], mean_cdc_hz_w_ic["hz_23_w_ic"],
                    time_grid_since_entry) %>% mutate(Panel = "Weighted + IC-corrected")
) %>% mutate(Panel = factor(Panel, levels = panel_levels_cdc))

cdc_x_labels <- c("Unweighted\n(IC-naive)", "Unweighted\n(IC-corrected)",
                   "Weighted\n(IC-naive)", "Weighted\n(IC-corrected)")

fig_cdc_A <- make_hz_boxplot(cdc_hz_long, ref_hz_df, fill_colors_cdc, "CDC",
                              "100 CDC runs — MIMIC cohort (simplified: no LOS filter, no admin censoring)",
                              x_breaks = c(1, 2, 3, 4), x_labels = cdc_x_labels, x_limits = c(0.4, 4.6))
print(fig_cdc_A)

fig_cdc_B <- plot_spp(spp_cdc, "CDC — Stacked Probability Plot (mean hazards across 100 runs)")
print(fig_cdc_B)

fig_cdc_C <- make_hr_boxplot(cdc_hr_long, ref_hr_df, fill_colors_cdc, "CDC",
                              "100 CDC runs — MIMIC cohort (simplified)",
                              x_breaks = c(1, 2, 3, 4), x_labels = cdc_x_labels, x_limits = c(0.4, 4.6))
print(fig_cdc_C)

fig_cdc_D <- make_ratio_boxplot(cdc_ratio_long, ref_results$ratio_est,
                                ref_results$ratio_CI_L, ref_results$ratio_CI_U,
                                fill_colors_cdc, "CDC", "100 CDC runs — MIMIC cohort (simplified)",
                                x_breaks = c(1, 2, 3, 4), x_labels = cdc_x_labels, x_limits = c(0.4, 4.6))
print(fig_cdc_D)

cat("SECTION 4 COMPLETE\n")
cat("============================================================\n\n")

############################################################
# HISTOGRAM: SAMPLED PATIENTS PER PPS DRAW
############################################################

pps_n_sampled <- sapply(seq_len(n_runs), function(i) {
  set.seed(123 + i)
  pps_dates <- runif(nrow(mimic_0), min = pps_date_min, max = pps_date_max)  # SIMPLIFIED: data-driven max
  sum(mimic_0$los_time > pps_dates)
})

hist(pps_n_sampled,
     main = "Distribution of PPS sample size across 100 runs",
     xlab = "Number of sampled patients",
     breaks = 15)

############################################################
# OBJECTS AVAILABLE AFTER RUNNING THIS SCRIPT:
#
#   mimic_0       : analysis cohort (post consistency-exclusion cleaning,
#                   no LOS filter, no admin censoring)
#   pps_date_max  : data-driven upper bound (max observed LOS) used for
#                   survey-date sampling and the SPP plot horizon
#   ref_results   : full-cohort reference MSM fits, Q, ratio, HRs
#   PPS_results / CDC_results : lists of 100 per-protocol run results
#   pps_hz / cdc_hz            : wide hazard data frames (100 runs each)
#
#   FIGURES:
#   fig_pps_A / fig_cdc_A : hazard rate boxplots (lambda12, 13, 23)
#   fig_pps_B / fig_cdc_B : stacked probability plots (mean hazards)
#   fig_pps_C / fig_cdc_C : covariate HR boxplots (age, sex)
#   fig_pps_D / fig_cdc_D : ratio lambda23/lambda13 boxplots
#
#   CDC reports FOUR methods: Unweighted+IC-naive, Unweighted+IC-corrected,
#   Weighted+IC-naive, Weighted+IC-corrected.
############################################################

############################################################
# EXCESS LENGTH OF STAY DUE TO HAI (Wolkewitz et al. 2017, Table 1)
#
# Adapted to this 3-state illness-death model (no death/discharge split,
# state 3 = combined discharge endpoint):
#   excess LOS = (1 / (lambda12 + lambda13)) * (lambda13 / lambda23 - 1)
############################################################

excess_los <- function(h12, h13, h23) {
  (1 / (h12 + h13)) * (h13 / h23 - 1)
}

cat("--- Excess length of stay due to HAI (days) ---\n")
cat(sprintf("  Full Cohort (Reference)     : %.2f\n",
            excess_los(ref_Q["1","2"], ref_Q["1","3"], ref_Q["2","3"])))
cat(sprintf("  Unweighted + IC-naive       : %.2f\n",
            excess_los(mean_cdc_hz_uw_noic["hz_12_uw_noic"], mean_cdc_hz_uw_noic["hz_13_uw_noic"], mean_cdc_hz_uw_noic["hz_23_uw_noic"])))
cat(sprintf("  Unweighted + IC-corrected   : %.2f\n",
            excess_los(mean_cdc_hz_uw_ic["hz_12_uw_ic"], mean_cdc_hz_uw_ic["hz_13_uw_ic"], mean_cdc_hz_uw_ic["hz_23_uw_ic"])))
cat(sprintf("  Weighted + IC-naive         : %.2f\n",
            excess_los(mean_cdc_hz_w_noic["hz_12_w_noic"], mean_cdc_hz_w_noic["hz_13_w_noic"], mean_cdc_hz_w_noic["hz_23_w_noic"])))
cat(sprintf("  Weighted + IC-corrected     : %.2f\n",
            excess_los(mean_cdc_hz_w_ic["hz_12_w_ic"], mean_cdc_hz_w_ic["hz_13_w_ic"], mean_cdc_hz_w_ic["hz_23_w_ic"])))
