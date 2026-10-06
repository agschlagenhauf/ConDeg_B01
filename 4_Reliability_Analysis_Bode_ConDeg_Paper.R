################################################################################
# reliability ConDeg Paper Bode et al. https://doi.org/10.31234/osf.io/vnzqu_v1
# Author: Erik Lukas Bode (erik-lukas.bode@charite.de) 
#
# Reported method: stability of contingency sensitivity was assessed using
# odd-even split-half reliability (Pronk et al., 2022) on individual DeltaP
# slopes, using a permuted approach (5,000 iterations, stratified by
# participant x reinforcer type x DeltaP level) with Spearman-Brown correction.
################################################################################

rm(list = ls())

libs <- c("dplyr", "tidyr", "stringr", "ggplot2", "tibble")
sapply(libs, require, character.only = TRUE)

# ------------------------------------------------------------------------------
# Set your own paths here.
# NOTE: raw and processed data are not included with this script (see data
# availability statement); this code documents the analysis pipeline as run
# on the original dataset.
# ------------------------------------------------------------------------------

EXCL_DIR <- "path/to/final_RData_excluded"

# ------------------------------------------------------------------------------
# Load frozen / excluded datasets
# ------------------------------------------------------------------------------
load(file.path(EXCL_DIR, "redcap.excluded.RData"))        # redcap.excluded
load(file.path(EXCL_DIR, "cd.task.RData"))                # cd.task
load(file.path(EXCL_DIR, "df_pleasantness.RData"))        # df_pleasantness

# ------------------------------------------------------------------------------
to_num_id <- function(x) as.numeric(as.character(x))
cd.task$participant          <- to_num_id(cd.task$participant)
redcap.excluded$participant  <- to_num_id(redcap.excluded$participant)
df_pleasantness$participant  <- to_num_id(df_pleasantness$participant)

# ------------------------------------------------------------------------------
# Define reinforcer condition from blockV
# blockV 1 = juice, 2 = alcohol (double check later)
# ------------------------------------------------------------------------------
cd.task <- cd.task %>%
  mutate(
    reinforcer_type = case_when(
      blockV == 1 ~ "juice",
      blockV == 2 ~ "alcohol",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(cd.task$reinforcer_type))) {
  warning("Some rows have blockV not in {1,2}. Check blockV coding.")
}

# ------------------------------------------------------------------------------
# Pleasantness: one row per participant × reinforcer_type
# ------------------------------------------------------------------------------
df_pleasantness_clean <- df_pleasantness %>%
  distinct(participant, type, mean_pleasantness) %>%
  mutate(
    mean_pleasantness_z = as.numeric(scale(mean_pleasantness)),
    reinforcer_type = as.character(type)  # assumes type uses same labels, e.g., "alcohol"/"juice"
  ) %>%
  dplyr::select(participant, reinforcer_type, mean_pleasantness, mean_pleasantness_z)

# ------------------------------------------------------------------------------
# BIN-level response rates
# Unit for later splitting: participant × reinforcer_type × runV × trialV
# DP is carried along, and should be constant within runV (usually)
# response_rate = count(choice==1) / count(valid choice in {0,1})
# ------------------------------------------------------------------------------
cd.rr_bin <- cd.task %>%
  group_by(participant, aud_group, reinforcer_type, runV, trialV, DP) %>%
  summarise(
    x1 = sum(choice == 1, na.rm = TRUE),
    x2 = sum(choice %in% c(0, 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    response_rate = ifelse(x2 > 0, x1 / x2, NA_real_)
  ) %>%
  dplyr::select(-x1, -x2)

# ------------------------------------------------------------------------------
# Merge covariates: age + pleasantness
# ------------------------------------------------------------------------------
cd.rr_bin <- cd.rr_bin %>%
  left_join(
    redcap.excluded %>% dplyr::select(participant, age),
    by = "participant"
  ) %>%
  left_join(
    df_pleasantness_clean,
    by = c("participant", "reinforcer_type")
  ) %>%
  mutate(
    agez = as.numeric(scale(age)),
    aud_group = as.factor(aud_group),
    reinforcer_type = factor(reinforcer_type, levels = c("juice", "alcohol"))
  )

# ------------------------------------------------------------------------------
# Quick diagnostics prints
# ------------------------------------------------------------------------------
message("Prepared cd.rr_bin")
message("Rows: ", nrow(cd.rr_bin))
message("Participants: ", n_distinct(cd.rr_bin$participant))
message("Conditions: ", paste(levels(cd.rr_bin$reinforcer_type), collapse = ", "))
message("DP values: ", paste(sort(unique(cd.rr_bin$DP)), collapse = ", "))

################################################################################
# ConDeg: Permutated (no-replacement) split-half reliability of ΔP slope
# Stratified within participant × reinforcer_type × DP
# Outputs:
# - per-condition mean r and Spearman–Brown corrected r_SB
# - overall (both conditions pooled) mean r and r_SB
################################################################################

# --------------------------------------------------------------------------
# 0) Sanity check: even bin counts within strata
# --------------------------------------------------------------------------
uneven <- cd.rr_bin %>%
  group_by(participant, reinforcer_type, DP) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n %% 2 != 0)

# Optional deeper check (usually unnecessary if above is OK)
cd.rr_bin %>%
  group_by(participant, reinforcer_type, DP, runV) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n %% 2 != 0) %>%
  print(n = Inf)

if (nrow(uneven) > 0) {
  warning("Uneven bin counts in some strata — permutated split will be imbalanced!")
  print(uneven)
} else {
  message("All strata have even bin counts. Permutated split OK.")
}

# --------------------------------------------------------------------------
# 1) Permutated splits
# --------------------------------------------------------------------------
set.seed(123)
n_iter <- 5000

r_perm_juice   <- numeric(n_iter)
r_perm_alcohol <- numeric(n_iter)
r_perm_overall <- numeric(n_iter)

# will hold slopes from the last iteration for scatter plot
slopes_last <- NULL

get_slope <- function(df) {
  df <- df %>% filter(!is.na(response_rate), !is.na(DP))
  if (nrow(df) < 5 || n_distinct(df$DP) < 2) return(NA_real_)
  coef(lm(response_rate ~ DP, data = df))[["DP"]]
}

for (i in 1:n_iter) {
  
  # 1) Stratified permutated split within participant × reinforcer_type × DP
  df_split <- cd.rr_bin %>%
    group_by(participant, reinforcer_type, DP) %>%
    mutate(
      n_g = n(),
      split_half = sample(rep(c("A", "B"), each = n_g / 2))
    ) %>%
    ungroup() %>%
    select(-n_g)
  
  # 2) Slopes per participant × condition × half
  slopes_iter <- df_split %>%
    group_by(participant, reinforcer_type, split_half) %>%
    summarise(
      slope = get_slope(pick(everything())),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = split_half, values_from = slope)
  
  # store LAST iteration slopes for plotting
  if (i == n_iter) slopes_last <- slopes_iter
  
  # 3) Correlations per condition + overall
  sj <- slopes_iter %>% filter(reinforcer_type == "juice")
  sa <- slopes_iter %>% filter(reinforcer_type == "alcohol")
  
  r_perm_juice[i]   <- cor(sj$A, sj$B, use = "complete.obs")
  r_perm_alcohol[i] <- cor(sa$A, sa$B, use = "complete.obs")
  r_perm_overall[i] <- cor(slopes_iter$A, slopes_iter$B, use = "complete.obs")
}

# --------------------------------------------------------------------------
# 2) Summaries + Spearman–Brown correction
# --------------------------------------------------------------------------
spearman_brown <- function(r) 2 * r / (1 + r)

summarise_r <- function(x) {
  m <- mean(x, na.rm = TRUE)
  tibble(
    mean_r = m,
    r_SB   = spearman_brown(m),
    p2.5   = quantile(x, 0.025, na.rm = TRUE),
    p97.5  = quantile(x, 0.975, na.rm = TRUE)
  )
}

results_perm <- bind_rows(
  summarise_r(r_perm_juice)   %>% mutate(reinforcer_type = "juice"),
  summarise_r(r_perm_alcohol) %>% mutate(reinforcer_type = "alcohol"),
  summarise_r(r_perm_overall) %>% mutate(reinforcer_type = "overall")
) %>%
  select(reinforcer_type, mean_r, r_SB, p2.5, p97.5)

print(results_perm)

# --------------------------------------------------------------------------
# 3) Plots
# --------------------------------------------------------------------------

# Plot 1: distributions of split-half correlations across permutations
df_r <- tibble(
  iter = 1:n_iter,
  juice   = r_perm_juice,
  alcohol = r_perm_alcohol,
  overall = r_perm_overall
) %>%
  pivot_longer(-iter, names_to = "series", values_to = "r")

ggplot(df_r, aes(x = r)) +
  geom_histogram(bins = 40) +
  facet_wrap(~series, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Split-half correlation r",
    y = "Count",
    title = "Permutated split-half reliability, distribution across splits"
  )

# Plot 2: density version
ggplot(df_r, aes(x = r)) +
  geom_density() +
  facet_wrap(~series, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Split-half correlation r",
    y = "Density",
    title = "Permutated split-half reliability, density across splits"
  )

# Plot 3: A vs B participant slopes from the LAST permutation split
ggplot(slopes_last, aes(x = A, y = B)) +
  geom_point(alpha = .8, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~reinforcer_type) +
  theme_classic() +
  labs(
    x = "ΔP slope (half A, last split)",
    y = "ΔP slope (half B, last split)",
    title = "Permutated split-half: participant slopes (last split)"
  )

################################################################################
# Notes:
# - mean_r is the average split-half correlation across random splits
# - r_SB is Spearman–Brown corrected reliability for the full-length task
# - p2.5/p97.5 are percentiles of the split-half r distribution (not SB-corrected)
################################################################################