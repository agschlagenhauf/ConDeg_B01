################################################################################
# Tasting and drink pleasantness ConDeg Paper Bode et al. https://doi.org/10.31234/osf.io/vnzqu_v1
# Purpose: Analyse the tasting and drink pleasentness data
# Author: Erik Lukas Bode (erik-lukas.bode@charite.de) with inital work from Milena Musial and Claudia Ebrahimi
################################################################################

# -----------------------------------------------------------------------------#
# Packages and global settings
# -----------------------------------------------------------------------------#

rm(list = ls())

libs <- c(
  "lme4", "lmerTest", "lmerTest", "tidyverse", "stringr", "ggplot2", "cowplot",
  "viridis", "Hmisc", "lubridate", "ez", "ggpubr", "rstatix",
  "moments", "coin", "robustlmm", "ggdist", "dplyr", "patchwork", "car", "emmeans"
)
sapply(libs, require, character.only = TRUE)

# raincloudplots (optional, only needed if plots are used later)
#if (!require(remotes)) install.packages("remotes")
#remotes::install_github("erocoar/gghalves")
#remotes::install_github("jorvlan/raincloudplots")
library(raincloudplots)

# Standard error helper (used in plots / summaries)
se <- function(x) sd(x, na.rm = TRUE) / sqrt(length(x) - sum(is.nan(x)))

theme_set(theme_cowplot())
options(digits = 10)

# ------------------------------------------------------------------------------
# Set your own working directory / data folder here.
# NOTE: raw and processed data are not included with this script (see data
# availability statement); this code documents the analysis pipeline as run
# on the original dataset.
# ------------------------------------------------------------------------------

EXCL_DIR <- "path/to/final_RData_excluded"

# ------------------------------------------------------------------------------
# Load excluded datasets
# ------------------------------------------------------------------------------

load(file.path(EXCL_DIR, "taste.ratings.excluded.RData"))
load(file.path(EXCL_DIR, "tasting.excluded.RData"))
load(file.path(EXCL_DIR, "redcap.excluded.RData"))

# -----------------------------------------------------------------------------#
# Prepare group-specific dataframes (0 = HC, 1 = AUD)
# -----------------------------------------------------------------------------#

taste.ratings.excluded.hc  <- taste.ratings.excluded[taste.ratings.excluded$aud_group == 0, ]
taste.ratings.excluded.aud <- taste.ratings.excluded[taste.ratings.excluded$aud_group == 1, ]

tasting.excluded.hc  <- tasting.excluded[tasting.excluded$aud_group == 0, ]
tasting.excluded.aud <- tasting.excluded[tasting.excluded$aud_group == 1, ]

redcap.excluded.hc  <- redcap.excluded[redcap.excluded$aud_group == 0, ]
redcap.excluded.aud <- redcap.excluded[redcap.excluded$aud_group == 1, ]

# Ensure gender is numeric everywhere (for descriptives / tables)
# Coding: 1 = female, 2 = male (divers not present in excluded sample)
redcap.excluded$gender     <- as.numeric(redcap.excluded$gender)
redcap.excluded.hc$gender  <- as.numeric(redcap.excluded.hc$gender)
redcap.excluded.aud$gender <- as.numeric(redcap.excluded.aud$gender)

# -----------------------------------------------------------------------------#
# Demographics (descriptives + group comparisons)
# -----------------------------------------------------------------------------#
# -----------------------------#
# Sample size per group
# -----------------------------#
table(redcap.excluded$aud_group)  # 0 = HC, 1 = AUD

# -----------------------------#
# Age
# -----------------------------#
summary(redcap.excluded$age)
sd(redcap.excluded$age, na.rm = TRUE)

# Age difference between groups (Welch t-test)
t.test(age ~ aud_group, data = redcap.excluded)  # Welch default

# -----------------------------#
# Gender
# -----------------------------#
# Coding: 1 = female, 2 = male
table(redcap.excluded$gender)

# Gender counts by group
gender_tab <- table(redcap.excluded$aud_group, redcap.excluded$gender)
gender_tab

# Group difference test (Chi-square)
chisq.test(gender_tab)

# -----------------------------#
# AUD symptom count (aud_sum)
# -----------------------------#
summary(redcap.excluded$aud_sum)
sd(redcap.excluded$aud_sum, na.rm = TRUE)

# Group difference (Welch t-test)
t.test(aud_sum ~ aud_group, data = redcap.excluded)  # Welch default

# -----------------------------#
# AUDIT sum (audit_sum)
# -----------------------------#
summary(redcap.excluded$audit_sum)
sd(redcap.excluded$audit_sum, na.rm = TRUE)

# Group difference (Welch t-test)
t.test(audit_sum ~ aud_group, data = redcap.excluded)  # Welch default


# -----------------------------#
# PSS (Perceived Stress Scale)
# -----------------------------#
# pss_sum is computed from restricted REDCap PSS exports (raw item-level
# responses are not included, see data availability
# statement) and merged into redcap.excluded prior to this script.

# Descriptives per group
redcap.excluded %>%
  group_by(aud_group) %>%
  summarise(
    n    = sum(!is.na(pss_sum)),
    mean = mean(pss_sum, na.rm = TRUE),
    sd   = sd(pss_sum, na.rm = TRUE),
    min  = min(pss_sum, na.rm = TRUE),
    max  = max(pss_sum, na.rm = TRUE)
  )

# Welch t-test
t.test(pss_sum ~ aud_group, data = redcap.excluded)

# -----------------------------------------------------------------------------#
# Group-wise descriptives (HC / AUD separately)
# -----------------------------------------------------------------------------#

# -----------------------------#
# HC only
# -----------------------------#
summary(redcap.excluded.hc$age)
sd(redcap.excluded.hc$age, na.rm = TRUE)

table(redcap.excluded.hc$gender)  # gender counts in HC

table(redcap.excluded.hc$smoker, useNA = "ifany")    # smoker counts in HC

summary(redcap.excluded.hc$aud_sum)
sd(redcap.excluded.hc$aud_sum, na.rm = TRUE)

summary(redcap.excluded.hc$audit_sum)
sd(redcap.excluded.hc$audit_sum, na.rm = TRUE)

# -----------------------------#
# AUD only
# -----------------------------#
summary(redcap.excluded.aud$age)
sd(redcap.excluded.aud$age, na.rm = TRUE)

table(redcap.excluded.aud$gender)  # gender counts in AUD

table(redcap.excluded.aud$smoker, useNA = "ifany")   # smoker counts in AUD

summary(redcap.excluded.aud$aud_sum)
sd(redcap.excluded.aud$aud_sum, na.rm = TRUE)

summary(redcap.excluded.aud$audit_sum)
sd(redcap.excluded.aud$audit_sum, na.rm = TRUE)


# ------------------------------------------------------------------------------
# Smoking status by group
# ------------------------------------------------------------------------------

# Contingency table
smoking_tab <- table(redcap.excluded$aud_group,
                     redcap.excluded$smoker)

smoking_tab

# ------------------------------------------------------------------------------
# Group difference test (Chi-square or Fisher if needed)
# ------------------------------------------------------------------------------

smoking_test <- chisq.test(smoking_tab)

# If expected counts are small, use Fisher's exact test instead
if (any(smoking_test$expected < 5)) {
  smoking_test <- fisher.test(smoking_tab)
  message("Used Fisher's exact test (small expected cell counts).")
} else {
  message("Used Pearson's Chi-square test.")
}

smoking_test
# -----------------------------------------------------------------------------#
# (0) Tasting phase: ratings across all drinks (main model)
# Goal:
# - LMM: rating ~ timepoint * type * aud_group + age (z)
# - random effects: compare candidate structures, keep chosen model
# -----------------------------------------------------------------------------#

# Set your own output paths here before running
fig_dir <- "path/to/your/figures"
out_dir <- "path/to/your/output"

# -----------------------------------------------------------------------------#
# Prepare variables 
# -----------------------------------------------------------------------------#

# aud_group as factor for effect coding
taste.ratings.excluded$aud_group <- as.factor(taste.ratings.excluded$aud_group)

# effect coding (keep as in original script)
contrasts(taste.ratings.excluded$type)      <- c(-0.5, 0.5)
contrasts(taste.ratings.excluded$aud_group) <- c(-0.5, 0.5)

# timepoint is already a factor with levels:
# tasting, pre, block1, block2, block3, block4, block5
# Convert to 0..6 safely using factor order
taste.ratings.excluded$timepoint <- as.numeric(taste.ratings.excluded$timepoint) - 1L

# sanity check
stopifnot(all(na.omit(taste.ratings.excluded$timepoint) %in% 0:6))
table(taste.ratings.excluded$timepoint, useNA = "ifany")

# age z-standardized (within available sample)
taste.ratings.excluded$agez <- scale(taste.ratings.excluded$age)

# -----------------------------------------------------------------------------#
# Candidate LMMs (same fixed effects, increasing random effect complexity)
# -----------------------------------------------------------------------------#

model1 <- lmer(
  rating ~ timepoint * type * aud_group + agez + (1 | participant),
  data      = taste.ratings.excluded,
  na.action = na.omit
)
summary(model1)

model2 <- lmer(
  rating ~ timepoint * type * aud_group + agez + (timepoint | participant),
  data      = taste.ratings.excluded,
  na.action = na.omit
)
summary(model2)

# Main model (best)
model3 <- lmer(
  rating ~ timepoint * type * aud_group + agez + (timepoint + type | participant),
  data      = taste.ratings.excluded,
  na.action = na.omit,
  REML      = TRUE,
  control   = lmerControl(optimizer = "Nelder_Mead")
)
summary(model3)

# Convergence notes 
summary(model3)$optinfo$conv$lme4$messages

# Model comparison across candidate structures
anova(model1, model2, model3) #model 3 best wins

# Profile likelihood CIs 
ci_prof <- confint(model3, method = "profile")
ci_prof

# -----------------------------------------------------------------------------#
# Post-hoc: simple slopes for timepoint by type with emmeans
# -----------------------------------------------------------------------------#

slopes_by_type <- emtrends(model3, ~ type, var = "timepoint")
slopes_by_type
confint(slopes_by_type)
pairs(slopes_by_type)

# -----------------------------------------------------------------------------#
# Model fit: R2 (Nakagawa) for reporting
# -----------------------------------------------------------------------------#

library(performance)
r2 <- r2_nakagawa(model3)
r2

# -----------------------------------------------------------------------------#
# Alternative random effects structure
# NOTE: converges but singular, therefore not used
# -----------------------------------------------------------------------------#

model4 <- lmer(
  rating ~ timepoint * type * aud_group + agez + (timepoint * type | participant),
  data      = taste.ratings.excluded,
  na.action = na.omit
)
summary(model4)

model4_bobyqa <- lmer(
  rating ~ timepoint * type * aud_group + agez + (timepoint * type | participant),
  data      = taste.ratings.excluded,
  na.action = na.omit,
  control   = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)
summary(model4_bobyqa)


# -----------------------------------------------------------------------------#
# Assumption checks for model3 (diagnostics)
# -----------------------------------------------------------------------------#

library(DHARMa)
library(influence.ME)

sim_res <- simulateResiduals(model3)
plot(sim_res)
testUniformity(sim_res)
testDispersion(sim_res)
testOutliers(sim_res)

# Visual inspection
qqnorm(resid(model3)); qqline(resid(model3))
plot(model3, which = 1)

# Multicollinearity (VIF)
vif_values <- car::vif(model3)
vif_values

# Influential observations (participant level)
infl <- influence(model3, group = "participant")
plot(infl, which = "cook")
summary(infl)

# -----------------------------------------------------------------------------#
# Robustness check: robustlmm (confirmatory robustness)
# NOTE: robust model used to confirm classical LMM results, not primary inference
# -----------------------------------------------------------------------------#

library(robustlmm)

m_rob <- rlmer(
  rating ~ timepoint * type * aud_group + agez + (timepoint + type | participant),
  data      = taste.ratings.excluded,
  na.action = na.omit,
  method    = "DASvar",
  init      = model3,
  max.iter  = 400,
  rel.tol   = 1e-7
)
summary(m_rob)

# Wald CIs based on robust SEs 
co <- summary(m_rob)$coefficients
est <- co[, "Estimate"]
se  <- co[, "Std. Error"]
z <- qnorm(0.975)

wald_ci <- data.frame(
  term     = rownames(co),
  estimate = est,
  se       = se,
  CI_low   = est - z * se,
  CI_high  = est + z * se,
  row.names = NULL
)
wald_ci

# -----------------------------------------------------------------------------#
# "Hybrid" robust p-values (Geniole-Paper, see supplement robust modelling): df from lmerTest + robust t-values
# -----------------------------------------------------------------------------#

# Classical LMM fit for df extraction (keep separate, df only)
model3_df <- lmer(
  rating ~ timepoint * type * aud_group + agez + (timepoint + type | participant),
  data      = taste.ratings.excluded,
  na.action = na.omit,
  REML      = TRUE,
  control   = lmerControl(optimizer = "Nelder_Mead")
)

sm_df <- summary(model3_df)
df_tab <- data.frame(
  term = rownames(sm_df$coefficients),
  df   = sm_df$coefficients[, "df"],
  row.names = NULL
)

rob_tab <- as.data.frame(coef(summary(m_rob)))
rob_tab$term <- rownames(rob_tab)

out <- merge(rob_tab, df_tab, by = "term", all.x = TRUE)
out$p_wald_satt <- 2 * pt(-abs(out[["t value"]]), df = out$df)

hybrid_summary <- out[, c("term", "Estimate", "Std. Error", "t value", "df", "p_wald_satt")]
names(hybrid_summary) <- c("term", "estimate_rob", "se_rob", "t_rob", "df_from_lmer", "p_wald_satt")
hybrid_summary

# -----------------------------------------------------------------------------#
# Save mean pleasantness rating per person across timepoints for each reinforcer
# -----------------------------------------------------------------------------#

# 1) Compute reinforcer-specific mean per participant
taste.ratings.excluded <- taste.ratings.excluded %>%
  group_by(participant, type) %>%
  mutate(
    mean_pleasantness = mean(rating, na.rm = TRUE)
  ) %>%
  ungroup()

# 2) Keep long format but INCLUDE type
df_pleasantness <- taste.ratings.excluded %>%
  dplyr::select(participant, type, mean_pleasantness)

# 3) Save it
save(df_pleasantness, file = file.path(out_dir, "df_pleasantness.RData"))