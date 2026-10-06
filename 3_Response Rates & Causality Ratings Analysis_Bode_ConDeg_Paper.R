################################################################################
# Response rate / contingency rating analysis 
# ConDeg Paper Bode et al. https://doi.org/10.31234/osf.io/vnzqu_v1
# Author: Erik Lukas Bode (erik-lukas.bode@charite.de) with inital work from Milena Musial and Claudia Ebrahimi
#
# Goal:
# - compute response rates per ΔP and reinforcer type
# - main LMM on response_rate (DV)
# - robustness checks (DHARMa + robustlmm)
################################################################################

# -----------------------------------------------------------------------------#
# Packages and global settings
# -----------------------------------------------------------------------------#

rm(list = ls())

libs <- c(
  "dplyr", "ggplot2", "tidyverse", "ez", "schoRsch", "lme4", "knitr",
  "gridExtra", "paletteer", "gt", "ggtext", "WRS2", "moments", "rstatix",
  "cowplot", "ARTool", "MASS", "lmerTest", "robustlmm", "ggpubr",
  "emmeans", "performance", "DHARMa", "car", "influence.ME", "lmerTest"
)
sapply(libs, require, character.only = TRUE)

se <- function(x) sd(x, na.rm = TRUE) / sqrt(length(x) - sum(is.nan(x)))
theme_set(theme_cowplot())

# -----------------------------------------------------------------------------#
# Load final excluded data + saved random effects from taste model
# -----------------------------------------------------------------------------#

# Set your own paths here.
# NOTE: raw and processed data are not included with this script (see data
# availability statement); this code documents the analysis pipeline as run
# on the original dataset.

EXCL_DIR <- "path/to/final_RData_excluded"
out_dir  <- "path/to/model_outputs"
fig_dir  <- "path/to/figures"


load(file.path(EXCL_DIR,"taste.ratings.excluded.RData"))
load(file.path(EXCL_DIR,"tasting.excluded.RData"))
load(file.path(EXCL_DIR,"redcap.excluded.RData"))
load(file.path(EXCL_DIR,"cd.task.RData"))
load(file.path(EXCL_DIR,"cd.causality.RData"))
# plaesantness
load(file.path(EXCL_DIR,"df_pleasantness.RData"))


# -----------------------------------------------------------------------------#
# Harmonize ID types across datasets
# Goal: participant IDs numeric everywhere for joins
# -----------------------------------------------------------------------------#

cd.causality$participant          <- as.numeric(as.character(cd.causality$participant))
cd.task$participant               <- as.numeric(as.character(cd.task$participant))
taste.ratings.excluded$participant<- as.numeric(as.character(taste.ratings.excluded$participant))
tasting.excluded$participant      <- as.numeric(as.character(tasting.excluded$participant))
redcap.excluded$participant       <- as.numeric(as.character(redcap.excluded$participant))
df_pleasantness$participant        <- as.numeric(as.character(df_pleasantness$participant))


# -----------------------------------------------------------------------------#
# Optional: split CD datasets by group (only needed for descriptive plots)
# -----------------------------------------------------------------------------#

cd.task.hc      <- cd.task[cd.task$aud_group == 0, ]
cd.task.aud     <- cd.task[cd.task$aud_group == 1, ]
cd.causality.hc <- cd.causality[cd.causality$aud_group == 0, ]
cd.causality.aud<- cd.causality[cd.causality$aud_group == 1, ]

# -----------------------------------------------------------------------------#
# Response rates per ΔP
# Definition:
# - x1 = number of choices==1
# - x2 = number of valid choices (1 + 0)
# - response_rate = x1/x2
# -----------------------------------------------------------------------------#

#-----------------------response rates---------------------

# calculate response rates per delta(p)

cd.rr <- cd.task %>% group_by(participant, DP, type, Version, aud_sum, aud_group) %>% 
  dplyr::summarize(x1 = sum(choice == 1),
                   x2 = sum(choice == 1) + sum(choice == 0)) %>% 
  ungroup() %>% 
  mutate(response_rate = x1 / x2) %>% 
  dplyr::select(-c(x1, x2)) # %>% 

cd.rr.hc <- cd.task.hc %>% group_by(participant, DP, type, Version) %>% 
  dplyr::summarize(x1 = sum(choice == 1),
                   x2 = sum(choice == 1) + sum(choice == 0)) %>% 
  ungroup() %>% 
  mutate(response_rate = x1 / x2) %>% 
  dplyr::select(-c(x1, x2)) # %>% 

cd.rr.aud <- cd.task.aud %>% group_by(participant, DP, type, Version) %>% 
  dplyr::summarize(x1 = sum(choice == 1),
                   x2 = sum(choice == 1) + sum(choice == 0)) %>% 
  ungroup() %>% 
  mutate(response_rate = x1 / x2) %>% 
  dplyr::select(-c(x1, x2)) # %>% 


# -----------------------------------------------------------------------------#
# Prepare predictors / contrasts for modelling
# NOTE: set contrasts ONCE here (avoid re-setting later)
# -----------------------------------------------------------------------------#

cd.rr$aud_group <- as.factor(cd.rr$aud_group)

contrasts(cd.rr$type)      <- c(-0.5, 0.5)
contrasts(cd.rr$aud_group) <- c(-0.5, 0.5)

#make pleasantness one row per type x participant
df_pleasantness <- df_pleasantness %>%
  distinct(participant, type, mean_pleasantness)

df_pleasantness <- df_pleasantness %>% 
  mutate(
    mean_pleasantness_z = as.numeric(scale(mean_pleasantness))
  )


# Merge covariates age and mean pleasantness rating
cd.rr <- left_join(cd.rr, redcap.excluded[, c("participant", "age")], by = "participant")
cd.rr <- cd.rr %>%
  left_join(
    df_pleasantness %>% dplyr::select(participant, type, mean_pleasantness, mean_pleasantness_z),
    by = c("participant", "type")
  )

# age z-standardized
cd.rr$agez <- scale(cd.rr$age)

cd.rr %>%
  group_by(type, DP) %>%
  summarise(
    mean_rr = mean(response_rate, na.rm = TRUE),
    sd_rr = sd(response_rate, na.rm = TRUE),
    min_rr = min(response_rate, na.rm = TRUE),
    max_rr = max(response_rate, na.rm = TRUE)
  )

# -----------------------------------------------------------------------------#
# Main LMM: response_rate ~ DP * type * aud_group (+ covariates / preference)
# Random effects: compare candidate structures
# -----------------------------------------------------------------------------#

model1a <- lmer(response_rate ~ DP * type * aud_group + agez + (1 | participant),
                data = cd.rr, na.action = na.omit)
summary(model1a)

model2a <- lmer(response_rate ~ DP * type * aud_group + agez + (DP | participant),
                data = cd.rr, na.action = na.omit)
summary(model2a)

model3a <- lmer(response_rate ~ DP * type * aud_group + agez + (DP + type | participant),
                data = cd.rr, na.action = na.omit)
summary(model3a)

# best model
model4a <- lmer(
  response_rate ~ DP * type * aud_group + agez + (DP * type | participant),
  data = cd.rr,
  na.action = na.omit,
  REML = TRUE)

summary(model4a)
isSingular(model4a)

# Model comparison
anova(model1a, model2a, model3a, model4a)

# Profile CIs for used model (time intensive)
ci_prof_4a <- confint(model4a, method = "profile")
ci_prof_4a

# Post-hoc: simple slopes of ΔP within each reinforcer type
slopes_DP_type <- emtrends(model4a, ~ type, var = "DP")
slopes_DP_type
confint(slopes_DP_type)
pairs(slopes_DP_type)  # optional, OK for supplement

# R2 (Nakagawa) for used model
r2_4a <- performance::r2_nakagawa(model4a)
r2_4a

# -----------------------------------------------------------------------------#
# Diagnostics: DHARMa + influence
# -----------------------------------------------------------------------------#

sim_res <- DHARMa::simulateResiduals(model4a)
plot(sim_res)
DHARMa::testUniformity(sim_res)
DHARMa::testDispersion(sim_res)
DHARMa::testOutliers(sim_res)

qqnorm(resid(model4a)); qqline(resid(model4a))
plot(model4a, which = 1)

vif_values <- car::vif(model4a)
vif_values

infl <- influence.ME::influence(model4a, group = "participant")
plot(infl, which = "cook")
summary(infl)

# -----------------------------------------------------------------------------#
# Robustness check: robustlmm
# NOTE: robustness only, primary inference from classical LMM
# -----------------------------------------------------------------------------#

m_rob4a <- robustlmm::rlmer(
  response_rate ~ DP * type * aud_group + agez + (DP * type | participant),
  data = cd.rr,
  na.action = na.omit,
  method = "DASvar",
  init = model4a,
  max.iter = 400,
  rel.tol = 1e-7
)
summary(m_rob4a)

# Wald CIs based on robust SEs
co <- summary(m_rob4a)$coefficients
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

# Hybrid robust p-values (df from lmerTest; robust t from rlmer)
model4a_df <- lmer(
  response_rate ~ DP * type * aud_group + agez + (DP * type | participant),
  data = cd.rr,
  na.action = na.omit,
  REML = TRUE
)

sm_df <- summary(model4a_df)
df_tab <- data.frame(
  term = rownames(sm_df$coefficients),
  df   = sm_df$coefficients[, "df"],
  row.names = NULL
)

rob_tab <- as.data.frame(coef(summary(m_rob4a)))
rob_tab$term <- rownames(rob_tab)

out <- merge(rob_tab, df_tab, by = "term", all.x = TRUE)
out$p_wald_satt <- 2 * pt(-abs(out[["t value"]]), df = out$df)

hybrid_summary <- out[, c("term","Estimate","Std. Error","t value","df","p_wald_satt")]
names(hybrid_summary) <- c("term","estimate_rob","se_rob","t_rob","df_from_lmer","p_wald_satt")
hybrid_summary

# Check model with pleasantess for control
#sanity checks before
cd.rr %>% summarise(sum(is.na(mean_pleasantness_z))) #no NAs in mean pleasantness

cd.rr_sorted <- cd.rr %>% #visually checked that mean for alk is same as before 
  arrange(
    participant,
    desc(type == "Alk"),  # Alk rows first
    DP                   # then DP from low to high
  )

model4a_check <- lmer(
  response_rate ~ DP * type * aud_group + mean_pleasantness_z + agez + (DP * type | participant),
  data = cd.rr,
  na.action = na.omit,
  REML = TRUE
)
summary(model4a_check)
isSingular(model4a_check)

# Profile CIs for used model (time intensive)
ci_prof_4a_check <- confint(model4a_check, method = "profile")
ci_prof_4a_check

# R2 (Nakagawa) for used model
r2_4a_check <- performance::r2_nakagawa(model4a_check)
r2_4a_check

# Model comparison (for Model comparison set REML to FALSE first)

anova(
  update(model4a, REML = FALSE), ##model with pleasantness is not increasing model fit, also all effects remained unchanged, therefore not reported besides model fit
  update(model4a_check, REML = FALSE)
)  
# -----------------------------------------------------------------------------#
# revision additional checks
# -----------------------------------------------------------------------------#
#check how block order is saved within the cd.task
#blockV is the condition, runV is the block order (confusing)
cd.task %>%
  dplyr::select(participant, type, DP, runV, Version) %>%
  distinct() %>%
  arrange(participant, type, runV) %>%
  head(30)
##left join ruNV by participant, type, DP to make sure the order stays exactly the same
# Pull runV into cd.rr 
runV_lookup <- cd.task %>%
  dplyr::select(participant, type, DP, runV) %>%
  distinct()

cd.rr <- cd.rr %>%
  left_join(runV_lookup, by = c("participant", "type", "DP"))

# Verify it joined cleanly - should be no NAs
cd.rr %>% summarise(na_runV = sum(is.na(runV))) # no NAs

# Check participant 3 in both dataframes side by side
##perfectly aligned in both DFs
cd.task %>%
  filter(participant == 3) %>%
  dplyr::select(participant, type, DP, runV) %>%
  distinct() %>%
  arrange(type, runV)

cd.rr %>%
  filter(participant == 3) %>%
  dplyr::select(participant, type, DP, runV) %>%
  arrange(type, runV)

# Robustness check: block order (runV) as fixed effect
cd.rr$runVc <- scale(cd.rr$runV, center = TRUE, scale = FALSE) #center order to have -2 -1 0 1 2

#with the three way interaction
model4a_order <- lmer(
  response_rate ~ DP * type * aud_group +
    DP * runVc * type +
    agez +
    (DP * type | participant),
  data = cd.rr,
  na.action = na.omit,
  REML = TRUE
)
summary(model4a_order)
isSingular(model4a_order)

# Profile CIs
ci_prof_4a_order <- confint(model4a_order, method = "profile")
ci_prof_4a_order

# R2 (Nakagawa)
r2_4a_order <- performance::r2_nakagawa(model4a_order)
r2_4a_order

# Model comparison (REML = FALSE first)
anova(
  update(model4a, REML = FALSE),
  update(model4a_order, REML = FALSE)
)

# Random slope version for completeness (expected to fail)
model4a_order_ran <- lmer(
  response_rate ~ DP * type * aud_group +
    agez +
    runVc +
    runVc:type +
    runVc:DP +
    (DP * type + runVc | participant),
  data = cd.rr,
  na.action = na.omit,
  REML = TRUE
)
isSingular(model4a_order_ran)


###check for high responder / high baseline responder
# Quartile split of mean response rate
q25 <- quantile(mean_rr$mean_rr, 0.25)
q75 <- quantile(mean_rr$mean_rr, 0.75)

mean_rr <- mean_rr %>%
  mutate(rr_group = case_when(
    mean_rr <= q25 ~ "low_responder",
    mean_rr >= q75 ~ "high_responder",
    TRUE ~ "middle"
  ))

# Join aud_group for descriptives
mean_rr_groups <- mean_rr %>%
  left_join(cd.rr %>% distinct(participant, aud_group), by = "participant")

# Count per group and AUD status
mean_rr_groups %>%
  filter(rr_group != "middle") %>%
  count(rr_group, aud_group)

# Descriptives per group
mean_rr_groups %>%
  filter(rr_group != "middle") %>%
  group_by(rr_group) %>%
  summarise(
    n = n(),
    mean_rr = mean(mean_rr),
    sd_rr = sd(mean_rr),
    min_rr = min(mean_rr),
    max_rr = max(mean_rr)
  )

# -----------------------------------------------------------------------------#
# Exploratory model in high responders only (no group factor, simplified random effects)
# -----------------------------------------------------------------------------#

# Get high responder participant IDs
high_resp_ids <- mean_rr %>%
  filter(rr_group == "high_responder") %>%
  pull(participant)

# Filter cd.rr to high responders only
cd.rr.high <- cd.rr %>%
  filter(participant %in% high_resp_ids)

cat("N high responders:", length(unique(cd.rr.high$participant)), "\n")

# Simplified model — no group factor, random intercept only given small n
model_high_resp <- lmer(
  response_rate ~ DP * type + agez + (1 | participant),
  data = cd.rr.high,
  na.action = na.omit,
  REML = TRUE
)
summary(model_high_resp)
isSingular(model_high_resp)

# Profile CIs
ci_high_resp <- confint(model_high_resp, method = "profile")
ci_high_resp

# R2
r2_high_resp <- performance::r2_nakagawa(model_high_resp)
r2_high_resp

# -----------------------------------------------------------------------------#
# Save final response rate model outputs 
# -----------------------------------------------------------------------------#

save(model4a, file = file.path(out_dir, "cd_rr_model4a.RData"))
save(m_rob4a, file = file.path(out_dir, "cd_rr_model4a_robust.RData"))
write.csv(hybrid_summary, file.path(out_dir, "cd_rr_model4a_hybrid_pvals.csv"), row.names = FALSE)
write.csv(wald_ci, file.path(out_dir, "cd_rr_model4a_robust_wald_ci.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------#
# Causality ratings (CD causality)
# Goal:
# - harmonize rating scaling across dummycausal conditions
# - derive a single causality measure: rating2 - rating3 (dummycausal == 4)
# -----------------------------------------------------------------------------#

# -----------------------------#
# Sanity checks (optional)
# -----------------------------#
table(cd.causality$dummycausal, cd.causality$blockV)
table(cd.causality$dummycausal, cd.causality$participant)


# -----------------------------#
# Create derived causality variable (dummycausal == 4)
# rating = rating(dummy 2) - rating(dummy 3)
# NOTE: keeping your current logic
# -----------------------------#
df.causality4 <- dplyr::filter(cd.causality, dummycausal == 1)
df.causality4$dummycausal <- 4

df.causality4$rating <- cd.causality$rating[cd.causality$dummycausal == 2] -
  cd.causality$rating[cd.causality$dummycausal == 3]

# -----------------------------#
# Check raw scale range per item
# -----------------------------#
tapply(cd.causality$rating, cd.causality$dummycausal, range, na.rm = TRUE)
tapply(cd.causality$rating, cd.causality$dummycausal, summary)


# -----------------------------#
# Confirm the linear relationship holds exactly (rating.rec == transform(rating))
# -----------------------------#
check1 <- cd.causality$rating[cd.causality$dummycausal == 1] * 200 - 100
all.equal(check1, cd.causality$rating.rec[cd.causality$dummycausal == 1])

check2 <- cd.causality$rating[cd.causality$dummycausal == 2] * 100
all.equal(check2, cd.causality$rating.rec[cd.causality$dummycausal == 2])

check3 <- cd.causality$rating[cd.causality$dummycausal == 3] * 100
all.equal(check3, cd.causality$rating.rec[cd.causality$dummycausal == 3])

# -----------------------------#
# Check derived ΔP variable (dummycausal == 4) range
# -----------------------------#
range(df.causality4$rating, na.rm = TRUE)       # should be within [-1, 1]


# covariate
df.causality4$agez <- scale(df.causality4$age)
df.causality4 <- df.causality4 %>%
  left_join(
    df_pleasantness %>% dplyr::select(participant, type, mean_pleasantness, mean_pleasantness_z),
    by = c("participant", "type")
  )

#quick sanity check NA
df.causality4 %>% summarise( 
  n_missing_z = sum(is.na(mean_pleasantness_z))
)

# -----------------------------#
# Prepare contrasts / join predictors
# -----------------------------#
df.causality4$aud_group <- as.factor(df.causality4$aud_group)

contrasts(df.causality4$type)      <- c(-0.5, 0.5)
contrasts(df.causality4$aud_group) <- c(-0.5, 0.5)

# -----------------------------------------------------------------------------#
# MODEL ROAD (aud_group predictor)
# -----------------------------------------------------------------------------#

# (1) minimal random intercept model (singular fit) -> best model
model1g <- lmer(
  rating ~ DP * type * aud_group + agez + (1 | participant),
  data      = df.causality4,
  na.action = na.omit
)
summary(model1g)

# check fixed effects without random effects (check fixed effects because of singular fit)
model1g_base_check <- lm(rating ~ DP * type * aud_group + agez,
                         data = df.causality4, na.action = na.omit)
summary(model1g_base_check)

# (2) try adding random slope for DP (singular / does not converge)
model2g <- lmer(
  rating ~ DP * type * aud_group + agez + (DP | participant),
  data      = df.causality4,
  na.action = na.omit,
  control   = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)
summary(model2g)
isSingular(model2g)

# (3) try DP + type random slopes (singular / does not converge)
model3g_bobyqa <- lmer(
  rating ~ DP * type * aud_group + agez + (DP + type | participant),
  data      = df.causality4,
  na.action = na.omit,
  control   = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)
summary(model3g_bobyqa)
isSingular(model3g_bobyqa)

# (3c) remove correlation among random effects (singular / does not converge)
model3g_uncorr_bobyqa <- lmer(
  rating ~ DP * type * aud_group + agez + (0 + DP + type || participant),
  data      = df.causality4,
  na.action = na.omit,
  control   = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)
summary(model3g_uncorr_bobyqa)
isSingular(model3g_uncorr_bobyqa)

###adding quadratic term
df.causality4 <- df.causality4 %>%
  mutate(
    DP_c  = scale(DP, center = TRUE, scale = FALSE),
    DP_c2 = DP_c^2
  )

model1g_quad <- lmer(
  rating ~ (DP + DP_c2) * type * aud_group + agez +
    (1 | participant),
  data      = df.causality4,
  REML      = TRUE,
  na.action = na.omit
)
summary(model1g_quad)

anova(
  update(model1g, REML = FALSE),
  update(model1g_quad, REML = FALSE)
) ###model 1g without quatradtic term is the better fit.

# Model comparison overview (road visible)
anova(model1g, model2g, model3g_bobyqa, model3g_uncorr_bobyqa)

# Profile likelihood CIs (FINAL model)
ci_prof_1g <- confint(model1g, method = "profile")
ci_prof_1g

# Simple slopes of ΔP within each group (FINAL model)
slopes_DP_group <- emtrends(model1g, ~ aud_group, var = "DP")
slopes_DP_group
confint(slopes_DP_group)
pairs(slopes_DP_group)  # optional (supplement)

# R2 (FINAL model)
r2_1g <- performance::r2_nakagawa(model1g)
r2_1g


# -----------------------------------------------------------------------------#
# Diagnostics (FINAL model only)
# -----------------------------------------------------------------------------#

sim_res_1g <- DHARMa::simulateResiduals(model1g)
plot(sim_res_1g)
DHARMa::testUniformity(sim_res_1g)
DHARMa::testDispersion(sim_res_1g)
DHARMa::testOutliers(sim_res_1g)

qqnorm(resid(model1g)); qqline(resid(model1g))
plot(model1g, which = 1)

vif_values_1g <- car::vif(model1g)
vif_values_1g

infl_1g <- influence.ME::influence(model1g, group = "participant")
plot(infl_1g, which = "cook")
summary(infl_1g)

# -----------------------------------------------------------------------------#
# Robustness check (FINAL model only, aud_group)
# -----------------------------------------------------------------------------#

m_rob_1g <- robustlmm::rlmer(
  rating ~ DP * type * aud_group + agez + (1 | participant),
  data      = df.causality4,     
  na.action = na.omit,
  method    = "DAStau",
  init      = model1g,
  max.iter  = 400,
  rel.tol   = 1e-7
)
summary(m_rob_1g)

# robust Wald CIs (FINAL robust model)
co <- summary(m_rob_1g)$coefficients
est <- co[, "Estimate"]
se  <- co[, "Std. Error"]
z <- qnorm(0.975)

wald_ci_1g <- data.frame(
  term     = rownames(co),
  estimate = est,
  se       = se,
  CI_low   = est - z * se,
  CI_high  = est + z * se,
  row.names = NULL
)
wald_ci_1g

# Hybrid robust p-values (df from lmerTest model; robust t from rlmer)
model1g_df <- lmer(
  rating ~ DP * type * aud_group + agez + (1 | participant),
  data      = df.causality4,
  REML      = TRUE,
  na.action = na.omit
)

sm_df <- summary(model1g_df)
df_tab <- data.frame(
  term = rownames(sm_df$coefficients),
  df   = sm_df$coefficients[, "df"],
  row.names = NULL
)

rob_tab <- as.data.frame(coef(summary(m_rob_1g)))
rob_tab$term <- rownames(rob_tab)

out <- merge(rob_tab, df_tab, by = "term", all.x = TRUE)
out$p_wald_satt <- 2 * pt(-abs(out[["t value"]]), df = out$df)

hybrid_summary_1g <- out[, c("term","Estimate","Std. Error","t value","df","p_wald_satt")]
names(hybrid_summary_1g) <- c("term","estimate_rob","se_rob","t_rob","df_from_lmer","p_wald_satt")
hybrid_summary_1g

# check for pleasantness add mean taste rating
model1g_check <- lmer(
  rating ~ DP * type * aud_group + mean_pleasantness_z + agez + (1 | participant),
  data      = df.causality4,
  REML      = TRUE,
  na.action = na.omit
)
summary(model1g_check)
isSingular(model1g_check)

# Profile CIs for used model (time intensive)
ci_prof_1g_check <- confint(model1g_check, method = "profile")
ci_prof_1g_check

# R2 (Nakagawa) for used model
r2_1g_check <- performance::r2_nakagawa(model1g_check)
r2_1g_check

anova(
  update(model1g, REML = FALSE),
  update(model1g_check, REML = FALSE)
) ###model 1g without quatradtic term is the better fit.


# -----------------------------------------------------------------------------#
# Save final causality model outputs 
# -----------------------------------------------------------------------------#

save(model1g,    file = file.path(out_dir, "cd_causality_model1g.RData"))
save(m_rob_1g,   file = file.path(out_dir, "cd_causality_model1g_robust.RData"))

write.csv(hybrid_summary_1g,
          file.path(out_dir, "cd_causality_model1g_hybrid_pvals.csv"),
          row.names = FALSE)

write.csv(wald_ci_1g,
          file.path(out_dir, "cd_causality_model1g_robust_wald_ci.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------#
# (S) Supplement: Response rate predicted by causality judgment (CJ)
# Goal:
# - Merge mean causality judgment (CJ) onto response-rate data (cd.rr)
# - LMM: response_rate ~ CJz * aud_group * type + covariates + preference proxies
# -----------------------------------------------------------------------------#
# -----------------------------#
# Prepare CJ predictor (participant × DP × type)
# NOTE: df.causality4 must exist (constructed in causality section)
# -----------------------------#
cd.rr.cj <- cd.rr %>%
  left_join(
    df.causality4 %>%
      group_by(participant, DP, type, aud_group) %>%
      summarise(CJ = mean(rating, na.rm = TRUE), .groups = "drop"),
    by = c("participant", "DP", "type", "aud_group")
  )

# effect coding as in the other models
cd.rr.cj$type      <- factor(cd.rr.cj$type)
cd.rr.cj$aud_group <- factor(cd.rr.cj$aud_group)
contrasts(cd.rr.cj$type)      <- c(-0.5, 0.5)
contrasts(cd.rr.cj$aud_group) <- c(-0.5, 0.5)

# z-score CJ globally
cd.rr.cj$CJz <- scale(cd.rr.cj$CJ, center = TRUE, scale = TRUE)

# -----------------------------#
# FINAL model (supplement): RR by CJp
# -----------------------------#
model5rcj <- lmer(
  response_rate ~ CJz * aud_group * type +
    agez +
    (1 + CJz * type | participant),
  data      = cd.rr.cj,
  na.action = na.omit,
  REML      = TRUE
)
summary(model5rcj)
isSingular(model5rcj)

# Profile likelihood CIs (time intensive)
ci_prof_5rcj <- confint(model5rcj, method = "profile")
ci_prof_5rcj

# R2 for reporting (supplement)
r2_5rcj <- performance::r2_nakagawa(model5rcj, tolerance = 1e-10)
r2_5rcj

# -----------------------------#
# Diagnostics 
# -----------------------------#
sim_res_5rcj <- DHARMa::simulateResiduals(model5rcj)
plot(sim_res_5rcj)
DHARMa::testUniformity(sim_res_5rcj)
DHARMa::testDispersion(sim_res_5rcj)
DHARMa::testOutliers(sim_res_5rcj)

qqnorm(resid(model5rcj)); qqline(resid(model5rcj))
plot(model5rcj, which = 1)

vif_values_5rcj <- car::vif(model5rcj)
vif_values_5rcj

infl_5rcj <- influence.ME::influence(model5rcj, group = "participant")
plot(infl_5rcj, which = "cook")
summary(infl_5rcj)

# -----------------------------#
# Robustness check 
# -----------------------------#
m_5rcj <- robustlmm::rlmer(
  response_rate ~ CJz * aud_group * type +
    agez +
    (1 + CJz * type | participant),
  data      = cd.rr.cj,
  na.action = na.omit,
  method    = "DASvar",
  init      = model5rcj,
  max.iter  = 400,
  rel.tol   = 1e-7
)
summary(m_5rcj)

# robust Wald CIs
co <- summary(m_5rcj)$coefficients
est <- co[, "Estimate"]
se  <- co[, "Std. Error"]
z <- qnorm(0.975)

wald_ci_5rcj <- data.frame(
  term     = rownames(co),
  estimate = est,
  se       = se,
  CI_low   = est - z * se,
  CI_high  = est + z * se,
  row.names = NULL
)
wald_ci_5rcj

# hybrid robust p-values (df from lmerTest model, robust t from rlmer)
model5rcj_df <- lmer(
  response_rate ~ CJz * aud_group * type +
    agez +
    (1 + CJz * type | participant),
  data      = cd.rr.cj,
  na.action = na.omit,
  REML      = TRUE
)

sm_df <- summary(model5rcj_df)
df_tab <- data.frame(
  term = rownames(sm_df$coefficients),
  df   = sm_df$coefficients[, "df"],
  row.names = NULL
)

rob_tab <- as.data.frame(coef(summary(m_5rcj)))
rob_tab$term <- rownames(rob_tab)

out <- merge(rob_tab, df_tab, by = "term", all.x = TRUE)
out$p_wald_satt <- 2 * pt(-abs(out[["t value"]]), df = out$df)

hybrid_summary_5rcj <- out[, c("term","Estimate","Std. Error","t value","df","p_wald_satt")]
names(hybrid_summary_5rcj) <- c("term","estimate_rob","se_rob","t_rob","df_from_lmer","p_wald_satt")
hybrid_summary_5rcj

#check for pleasantness effect
model5rcj_check <- lmer(
  response_rate ~ CJz * aud_group * type +
    agez + mean_pleasantness_z +
    (1 + CJz * type | participant),
  data      = cd.rr.cj,
  na.action = na.omit,
  REML      = TRUE
)
summary(model5rcj_check)
isSingular(model5rcj_check)

##adapt it
# Profile CIs for used model (time intensive)
ci_prof_5rcj_check <- confint(model5rcj_check, method = "profile")
ci_prof_5rcj_check

# R2 (Nakagawa) for used model
r2_5rcj_check <- performance::r2_nakagawa(model5rcj_check, tolerance = 1e-10)
r2_5rcj_check

anova(
  update(model5rcj, REML = FALSE),
  update(model5rcj_check, REML = FALSE)
) 

# Save 
save(model5rcj, file = file.path(out_dir, "cd_rr_cj_model5rcj.RData"))
save(m_5rcj,    file = file.path(out_dir, "cd_rr_cj_model5rcj_robust.RData"))
write.csv(hybrid_summary_5rcj, file.path(out_dir, "cd_rr_cj_model5rcj_hybrid_pvals.csv"), row.names = FALSE)
write.csv(wald_ci_5rcj,        file.path(out_dir, "cd_rr_cj_model5rcj_robust_wald_ci.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------#
# Ratio score analysis (Contingency degradation)
#
# STRUCTURE:
# A) Ratio score ΔP .6 → .3  
# B) Ratio score ΔP .6 → .0  
# C) Supplement Figure: barplots side-by-side (ΔP .6→.3 and ΔP .6→.0)
#
# NOTE:
# Ratio score definition:
#   ratio = R(.6) / [R(.6) + R(degraded)]
# where "lower" is .3 or 0.0 depending on contrast
# -----------------------------------------------------------------------------#

# helper: standardize labels once (type + group)
harmonize_ratio_labels <- function(df) {
  df %>%
    mutate(
      type = dplyr::recode(as.character(type),
                           "Alk"  = "alcohol",
                           "Saft" = "juice",
                           .default = as.character(type)),
      aud_group = dplyr::case_when(
        aud_group %in% c("HC", "AUD") ~ aud_group,
        aud_group %in% c(0, "0")      ~ "HC",
        aud_group %in% c(1, "1")      ~ "AUD",
        TRUE                          ~ as.character(aud_group)
      )
    )
}

# helper: Cohen's dz for paired designs
cohens_dz <- function(x, y) {
  d <- x - y
  mean(d, na.rm = TRUE) / sd(d, na.rm = TRUE)
}

# -----------------------------------------------------------------------------#
# A) Ratio score ΔP .6 → .3
# -----------------------------------------------------------------------------#

df_ratio_03 <- cd.rr %>%
  filter(DP %in% c(0.6, 0.3)) %>%
  harmonize_ratio_labels()

# participant-level stable vars (aud_group, agez)
demo_vars_03 <- df_ratio_03 %>%
  group_by(participant) %>%
  summarise(
    aud_group = first(aud_group),
    agez      = first(agez),
    .groups   = "drop"
  )

# compute ratio per participant × type
df_ratio_wide_03 <- df_ratio_03 %>%
  dplyr::select(participant, type, DP, response_rate) %>%
  tidyr::pivot_wider(names_from = DP, values_from = response_rate, names_prefix = "DP_") %>%
  mutate(
    ratio_score = DP_0.6 / (DP_0.6 + DP_0.3)
  ) %>%
  tidyr::drop_na(ratio_score) %>%            # require both DP values
  left_join(demo_vars_03, by = "participant") %>%
  mutate(
    type      = factor(type, levels = c("juice", "alcohol")),
    aud_group = factor(aud_group, levels = c("HC", "AUD"))
  )

summary(df_ratio_wide_03$ratio_score)  # sanity check

# paired t-test alcohol vs juice (within-subject)
df_ratio_paired_03 <- df_ratio_wide_03 %>%
  dplyr::select(participant, type, ratio_score) %>%
  tidyr::pivot_wider(names_from = type, values_from = ratio_score) %>%
  tidyr::drop_na(juice, alcohol)

tt_03 <- t.test(df_ratio_paired_03$alcohol, df_ratio_paired_03$juice, paired = TRUE)
tt_03

dz_03 <- cohens_dz(df_ratio_paired_03$alcohol, df_ratio_paired_03$juice)
dz_03

# Means/SDs for reporting 
overall_03 <- df_ratio_wide_03 %>%
  group_by(type) %>%
  summarise(
    M  = mean(ratio_score, na.rm = TRUE),
    SD = sd(ratio_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(M = sprintf("%.3f", M), SD = sprintf("%.3f", SD))
overall_03

# -----------------------------------------------------------------------------#
# B) Ratio score ΔP .6 → .0 
# -----------------------------------------------------------------------------#

df_ratio_00 <- cd.rr %>%
  filter(DP %in% c(0.6, 0.0)) %>%
  harmonize_ratio_labels() %>%
  mutate(
    # stable DP label for pivoting (0.0 -> "0", 0.6 -> "0.6")
    DP_lbl = case_when(abs(DP) < 1e-12 ~ "0", TRUE ~ as.character(DP))
  )

demo_vars_00 <- df_ratio_00 %>%
  group_by(participant) %>%
  summarise(
    aud_group = first(aud_group),
    agez      = first(agez),
    .groups   = "drop"
  )

df_ratio_wide_00 <- df_ratio_00 %>%
  dplyr::select(participant, type, DP_lbl, response_rate) %>%
  tidyr::pivot_wider(names_from = DP_lbl, values_from = response_rate, names_prefix = "DP_") %>%
  mutate(
    ratio_score = DP_0.6 / (DP_0.6 + DP_0)
  ) %>%
  tidyr::drop_na(ratio_score) %>%
  left_join(demo_vars_00, by = "participant") %>%
  mutate(
    type      = factor(type, levels = c("juice", "alcohol")),
    aud_group = factor(aud_group, levels = c("HC", "AUD"))
  )

summary(df_ratio_wide_00$ratio_score)  # sanity check

df_ratio_paired_00 <- df_ratio_wide_00 %>%
  dplyr::select(participant, type, ratio_score) %>%
  tidyr::pivot_wider(names_from = type, values_from = ratio_score) %>%
  tidyr::drop_na(juice, alcohol)

tt_00 <- t.test(df_ratio_paired_00$alcohol, df_ratio_paired_00$juice, paired = TRUE)
tt_00

dz_00 <- cohens_dz(df_ratio_paired_00$alcohol, df_ratio_paired_00$juice)
dz_00

# Means/SDs for reporting (SUPPLEMENT)
overall_00 <- df_ratio_wide_00 %>%
  group_by(type) %>%
  summarise(
    M  = mean(ratio_score, na.rm = TRUE),
    SD = sd(ratio_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(M = sprintf("%.3f", M), SD = sprintf("%.3f", SD))
overall_00


# -----------------------------------------------------------------------------#
# C) Supplement Figure: barplots side-by-side (ΔP .6→.3 and ΔP .6→.0)
# -----------------------------------------------------------------------------#

# ensure consistent order
df_ratio_wide_03$type <- factor(df_ratio_wide_03$type, levels = c("juice","alcohol"))
df_ratio_wide_00$type <- factor(df_ratio_wide_00$type, levels = c("juice","alcohol"))

# helper to make the barplot (mean ± SD + significance label)
make_ratio_bar <- function(df_sum, p_value, title_str) {
  p_lab <- if (p_value < .001) "***" else if (p_value < .01) "**" else if (p_value < .05) "*" else "n.s."
  y_max <- max(df_sum$M + df_sum$SD, na.rm = TRUE) + 0.03
  
  ggplot(df_sum, aes(x = type, y = M, fill = type)) +
    geom_col(width = 0.6, color = "black") +
    geom_errorbar(aes(ymin = M - SD, ymax = M + SD), width = 0.15) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey55") +
    annotate("segment", x = 1, xend = 2, y = y_max, yend = y_max) +
    annotate("segment", x = 1, xend = 1, y = y_max, yend = y_max - 0.02) +
    annotate("segment", x = 2, xend = 2, y = y_max, yend = y_max - 0.02) +
    annotate("text", x = 1.5, y = y_max + 0.015, label = p_lab, size = 3) +
    scale_x_discrete(labels = c(juice = "Juice", alcohol = "Alcohol")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_fill_manual(values = c("juice" = "#0072B2", "alcohol" = "#D55E00")) +
    labs(x = NULL, y = "Ratio score", title = title_str) +
    theme_classic(base_size = 7) +
    theme(
      legend.position = "none",
      axis.text  = element_text(size = 6),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# summary tables for bars
sum_03 <- df_ratio_wide_03 %>%
  group_by(type) %>%
  summarise(M = mean(ratio_score, na.rm = TRUE),
            SD = sd(ratio_score, na.rm = TRUE),
            .groups = "drop")

sum_00 <- df_ratio_wide_00 %>%
  group_by(type) %>%
  summarise(M = mean(ratio_score, na.rm = TRUE),
            SD = sd(ratio_score, na.rm = TRUE),
            .groups = "drop")

p_bar_03 <- make_ratio_bar(sum_03, tt_03$p.value, expression(Delta * "P 0.6 " %->% " 0.3"))
p_bar_00 <- make_ratio_bar(sum_00, tt_00$p.value, expression(Delta * "P 0.6 " %->% " 0.0"))

Supp_Ratio <- ggpubr::ggarrange(
  p_bar_03, p_bar_00,
  ncol = 2,
  labels = c("a)", "b)"),
  font.label = list(size = 7, face = "bold"),
  align = "hv"
)

Supp_Ratio

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

message("Saving figures to: ", normalizePath(fig_dir))

# save (supp figure) - explicit path (no setwd)
ggsave(file.path(fig_dir, "Supp_RatioScores.png"),
       Supp_Ratio, width = 160, height = 80, units = "mm", dpi = 800)

ggsave(file.path(fig_dir, "Supp_RatioScores.tiff"),
       Supp_Ratio, width = 160, height = 80, units = "mm", dpi = 800, compression = "lzw")

ggsave(file.path(fig_dir, "Supp_RatioScores.pdf"),
       Supp_Ratio, width = 160, height = 80, units = "mm", device = cairo_pdf)