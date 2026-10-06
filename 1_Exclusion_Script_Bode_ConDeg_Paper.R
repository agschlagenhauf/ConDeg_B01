################################################################################
# Exclusion Script ConDeg Paper Bode et al. https://doi.org/10.31234/osf.io/vnzqu_v1
# Author: Erik Lukas Bode (erik-lukas.bode@charite.de)
# 
# Purpose:
# - Load pre-merged WP2 cohort files
# - Harmonize participant ID types and tasting timepoint coding
# - Define participant-level exclusions (predefined + missingness)
# - Apply exclusions consistently across all modality datasets
# - Save filtered datasets using legacy filenames expected by downstream scripts
#
# Notes:
# - Input files are curated ("*_104.RData") (misleading name its 106 unique IDs in total but mismatching)
# - Participant IDs can differ slightly across modality files due to missing data
################################################################################

library(dplyr)
rm(list=ls())

# ------------------------------------------------------------------------------
# Set your own working directory / data folder here.
# NOTE: raw and processed data are not included with this script (see data
# availability statement); this code documents the analysis pipeline as run
# on the original dataset.
# ------------------------------------------------------------------------------

DATA_DIR <- "path/to/R_dfs/final_merged_pilot+wp2"

# ------------------------------------------------------------------------------
# Loads (same relative paths everywhere)
# ------------------------------------------------------------------------------

load(file.path(DATA_DIR, "tasting_104.RData"))
load(file.path(DATA_DIR, "taste_ratings_complete_104.RData"))
load(file.path(DATA_DIR, "cd_104.RData"))
load(file.path(DATA_DIR, "redcap_104.RData"))


# Rename to short, explicit objects (KEEP NAMES)
cd.causality  <- cd.causality.comb
cd.task       <- cd.task.comb
taste.ratings <- taste.ratings.comb
tasting       <- tasting.comb
redcap        <- redcap.comb

# -----------------------------------------------------------------------------#
# Harmonize ID types and tasting timepoint coding
# -----------------------------------------------------------------------------#

# Ensure participant IDs are numeric across all datasets
cd.causality$participant   <- as.numeric(as.character(cd.causality$participant))
cd.task$participant        <- as.numeric(as.character(cd.task$participant))
taste.ratings$participant  <- as.numeric(as.character(taste.ratings$participant))
tasting$participant        <- as.numeric(as.character(tasting$participant))
redcap$participant         <- as.numeric(as.character(redcap$participant))

# Recode "number" to start at 0 for baseline tasting (KEEP LOGIC)
# Original coding: 1..7  ->  recoded: 0..6
taste.ratings <- taste.ratings %>%
  mutate(number = dplyr::recode(
    number,
    "1" = 0, "2" = 1, "3" = 2, "4" = 3, "5" = 4, "6" = 5, "7" = 6
  ))

cat("\nUnique participants in taste.ratings:", length(unique(taste.ratings$participant)), "\n")

# -----------------------------------------------------------------------------#
# Document modality-specific ID mismatches (informative only)
# -----------------------------------------------------------------------------#
redcap_ids <- unique(redcap$participant)
taste_ids  <- unique(taste.ratings$participant)

missing_in_taste <- setdiff(redcap_ids, taste_ids)  # present in redcap, absent in taste.ratings
extra_in_taste   <- setdiff(taste_ids, redcap_ids)  # present in taste.ratings, absent in redcap

cat("\nParticipants present in REDCap but absent in taste.ratings:\n")
print(sort(missing_in_taste))

cat("\nParticipants present in taste.ratings but absent in REDCap:\n")
print(sort(extra_in_taste))


# ------------------------------------------------------------------------------
# Harmonize gender coding across projects
#
# Main studies:
#   1 = female, 2 = male, 3 = divers
# Pilot study:
#   1 = male, 2 = female, 3 = divers
#
# We recode ONLY for pilot == 1 so that the harmonized variable matches:
#   1 = female, 2 = male, 3 = divers
# ------------------------------------------------------------------------------

# assume:
#   gender is numeric/integer with values 1,2,3 (may include NA)
redcap$gender <- as.numeric(redcap$gender)
#   1 for pilot, 0 for not pilot

# keep original for traceability
redcap <- redcap %>%
  mutate(
    gender_raw = gender,
    gender = case_when(
      pilot == 1 & gender_raw == 1 ~ 2,
      pilot == 1 & gender_raw == 2 ~ 1,
      TRUE ~ gender_raw
    ),
    gender_f = factor(
      gender,
      levels = c(1, 2, 3),
      labels = c("weiblich", "männlich", "divers")
    )
  )

# Check pilot vs non-pilot distributions BEFORE and AFTER
redcap %>%
  dplyr::count(pilot, gender_raw) %>%
  tidyr::complete(pilot, gender_raw = c(1,2,3), fill = list(n = 0)) %>%
  dplyr::arrange(pilot, gender_raw)

redcap %>%
  dplyr::count(pilot, gender) %>%
  tidyr::complete(pilot, gender = c(1,2,3), fill = list(n = 0)) %>%
  dplyr::arrange(pilot, gender)


# -----------------------------------------------------------------------------#
# Define participant-level exclusions (KEEP LOGIC + NAMES)
# -----------------------------------------------------------------------------#

# Predefined taste-based criteria at baseline tasting (number == 0)
excluded_by_criterion_1 <- unique(
  taste.ratings$participant[taste.ratings$number == 0 & taste.ratings$rating < 0.5]
)

excluded_by_criterion_2 <- unique(
  taste.ratings$participant[taste.ratings$number == 0 & abs(taste.ratings$diff_AS) > 0.4]
)

ids.excluded_predefined <- union(excluded_by_criterion_1, excluded_by_criterion_2)

# Additional exclusions due to missing / incomplete data (KEEP LIST)
ids.excluded_missingness <- c(12390, 2, 5, 7, 18)

# Final master exclusion list
ids.excluded <- union(ids.excluded_predefined, ids.excluded_missingness)

# Diagnostics (audit trail)
cat("\nIDs excluded by Criterion 1 (baseline rating < 0.5):\n")
print(sort(excluded_by_criterion_1))

cat("\nIDs excluded by Criterion 2 (baseline abs(diff_AS) > 0.4):\n")
print(sort(excluded_by_criterion_2))

cat("\nAdditional exclusions due to missing/incomplete data:\n")
print(sort(ids.excluded_missingness))

cat("\nFinal combined exclusion list (participant-level):\n")
print(sort(ids.excluded))

cat("\nTotal excluded IDs (participant-level master list):", length(ids.excluded), "\n")

# -----------------------------------------------------------------------------#
# Save FULL dataset (before taste-based exclusions)
# Only missingness exclusions applied (ids.excluded_missingness)
# Used for the full-responder causality and response-rate models reported
# in the manuscript
# -----------------------------------------------------------------------------#

out_dir <- "path/to/R_dfs/final_RData_excluded"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Apply only missingness exclusions to get "full" datasets
cd.task.full      <- dplyr::filter(cd.task,      !participant %in% ids.excluded_missingness)
cd.causality.full <- dplyr::filter(cd.causality, !participant %in% ids.excluded_missingness)
redcap.full       <- dplyr::filter(redcap,        !participant %in% ids.excluded_missingness)
taste.ratings.full <- dplyr::filter(taste.ratings, !participant %in% ids.excluded_missingness)

cat("\nUnique participants in FULL dataset (missingness only):\n")
cat("  cd.task.full     :", length(unique(cd.task.full$participant)), "\n")
cat("  cd.causality.full:", length(unique(cd.causality.full$participant)), "\n")
cat("  redcap.full      :", length(unique(redcap.full$participant)), "\n")
cat("  taste.ratings.full:", length(unique(taste.ratings.full$participant)), "\n")

save(cd.task.full,       file = file.path(out_dir, "cd.task.full.RData"))
save(cd.causality.full,  file = file.path(out_dir, "cd.causality.full.RData"))
save(redcap.full,        file = file.path(out_dir, "redcap.full.RData"))
save(taste.ratings.full, file = file.path(out_dir, "taste.ratings.full.RData"))

# Also save the excluded IDs for reference
ids.excluded_taste <- ids.excluded_predefined  # taste-based only

save(ids.excluded_taste,
     ids.excluded_missingness,
     ids.excluded,
     file = file.path(out_dir, "exclusion_ids.RData"))

cat("\nSaved full-sample datasets:\n",
    " - cd.task.full.RData\n",
    " - cd.causality.full.RData\n",
    " - redcap.full.RData\n",
    " - taste.ratings.full.RData\n",
    " - exclusion_ids.RData\n", sep = "")

# -----------------------------------------------------------------------------#
# Save TASTE-EXCLUDED-ONLY dataset
# Only the participants removed by taste criteria (not missingness)
# Used for the primary/full/taste-excluded-only sensitivity comparison
# reported in the supplement (response rate and causality models only —
# not used for the pleasantness/taste-rating analysis)
# -----------------------------------------------------------------------------#

cd.task.taste_excl_only      <- dplyr::filter(cd.task,      participant %in% ids.excluded_predefined)
cd.causality.taste_excl_only <- dplyr::filter(cd.causality, participant %in% ids.excluded_predefined)
redcap.taste_excl_only       <- dplyr::filter(redcap,       participant %in% ids.excluded_predefined)

cat("\nUnique participants in TASTE-EXCLUDED-ONLY dataset:\n")
cat("  cd.task.taste_excl_only      :", length(unique(cd.task.taste_excl_only$participant)), "\n")
cat("  cd.causality.taste_excl_only :", length(unique(cd.causality.taste_excl_only$participant)), "\n")
cat("  redcap.taste_excl_only       :", length(unique(redcap.taste_excl_only$participant)), "\n")

# Group split in excluded-only (informative for power assessment)
redcap.taste_excl_only %>% count(aud_group)

save(cd.task.taste_excl_only,      file = file.path(out_dir, "cd.task.taste_excl_only.RData"))
save(cd.causality.taste_excl_only, file = file.path(out_dir, "cd.causality.taste_excl_only.RData"))
save(redcap.taste_excl_only,       file = file.path(out_dir, "redcap.taste_excl_only.RData"))

cat("\nSaved taste-excluded-only datasets:\n",
    " - cd.task.taste_excl_only.RData\n",
    " - cd.causality.taste_excl_only.RData\n",
    " - redcap.taste_excl_only.RData\n", sep = "")

# -----------------------------------------------------------------------------#
# Apply exclusions to each dataset
# -----------------------------------------------------------------------------#

# Cache BEFORE counts (so the "before" print is actually before)
n_taste_before <- length(unique(taste.ratings$participant))
n_tasting_before <- length(unique(tasting$participant))
n_redcap_before <- length(unique(redcap$participant))
n_cd_task_before <- length(unique(cd.task$participant))
n_cd_caus_before <- length(unique(cd.causality$participant))

# Apply participant-level exclusion list
taste.ratings.excluded <- dplyr::filter(taste.ratings, !participant %in% ids.excluded)
tasting.excluded       <- dplyr::filter(tasting,       !participant %in% ids.excluded)
redcap.excluded        <- dplyr::filter(redcap,        !participant %in% ids.excluded)

# Legacy behavior: overwrite cd.task and cd.causality in place
cd.task      <- dplyr::filter(cd.task,      !participant %in% ids.excluded)
cd.causality <- dplyr::filter(cd.causality, !participant %in% ids.excluded)

# Dataset-specific N checks
cat("\nUnique participants BEFORE exclusion:\n")
cat("  taste.ratings:", n_taste_before, "\n")
cat("  tasting      :", n_tasting_before, "\n")
cat("  redcap       :", n_redcap_before, "\n")
cat("  cd.task      :", n_cd_task_before, "\n")
cat("  cd.causality :", n_cd_caus_before, "\n")

cat("\nUnique participants AFTER exclusion:\n")
cat("  taste.ratings:", length(unique(taste.ratings.excluded$participant)), "\n")
cat("  tasting      :", length(unique(tasting.excluded$participant)), "\n")
cat("  redcap       :", length(unique(redcap.excluded$participant)), "\n")
cat("  cd.task      :", length(unique(cd.task$participant)), "\n")
cat("  cd.causality :", length(unique(cd.causality$participant)), "\n")

# -----------------------------------------------------------------------------#
# Save filtered datasets
# -----------------------------------------------------------------------------#

save(taste.ratings.excluded, file = file.path(out_dir, "taste.ratings.excluded.RData"))
save(tasting.excluded,       file = file.path(out_dir, "tasting.excluded.RData"))
save(redcap.excluded,        file = file.path(out_dir, "redcap.excluded.RData"))

# Legacy filenames expected by downstream scripts
save(cd.task,      file = file.path(out_dir, "cd.task.RData"))
save(cd.causality, file = file.path(out_dir, "cd.causality.RData"))

cat("\nSaved excluded datasets:\n",
    " - taste.ratings.excluded.RData\n",
    " - tasting.excluded.RData\n",
    " - redcap.excluded.RData\n",
    " - cd.task.RData\n",
    " - cd.causality.RData\n", sep = "")