################################################################################
# Supplement: Response Rate Model & Causality Rating Model Full Sample
# Author: Erik Lukas Bode (erik-lukas.bode@charite.de) 
# ConDeg Paper Bode et al. https://doi.org/10.31234/osf.io/vnzqu_v1
################################################################################

rm(list = ls())

libs <- c(
  "dplyr", "ggplot2", "tidyverse", "lme4", "lmerTest",
  "emmeans", "performance", "cowplot", "ggpubr", "ragg"
)
sapply(libs, require, character.only = TRUE)

se <- function(x) sd(x, na.rm = TRUE) / sqrt(length(x) - sum(is.nan(x)))
theme_set(theme_cowplot())

# -----------------------------------------------------------------------------#
# Set your own paths here.
# NOTE: raw and processed data are not included with this script (see data
# availability statement); this code documents the analysis pipeline as run
# on the original dataset.
# -----------------------------------------------------------------------------#
EXCL_DIR <- "path/to/final_RData_excluded"
fig_dir  <- "path/to/figures"

# -----------------------------------------------------------------------------#
# Load datasets
# -----------------------------------------------------------------------------#
load(file.path(EXCL_DIR, "cd.task.RData"))                  # primary (n=81)
load(file.path(EXCL_DIR, "cd.task.full.RData"))             # full (n=101)
load(file.path(EXCL_DIR, "cd.task.taste_excl_only.RData"))  # excluded only (n=20)
load(file.path(EXCL_DIR, "redcap.excluded.RData"))
load(file.path(EXCL_DIR, "redcap.full.RData"))
load(file.path(EXCL_DIR, "redcap.taste_excl_only.RData"))

# -----------------------------------------------------------------------------#
# Helper: compute response rates + prepare predictors
# -----------------------------------------------------------------------------#
prepare_rr <- function(cd_task, redcap_df) {
  
  # Harmonize IDs
  cd_task$participant  <- as.numeric(as.character(cd_task$participant))
  redcap_df$participant <- as.numeric(as.character(redcap_df$participant))
  
  # Response rates per DP x type x participant
  cd_rr <- cd_task %>%
    group_by(participant, DP, type, Version, aud_sum, aud_group) %>%
    dplyr::summarize(
      x1 = sum(choice == 1),
      x2 = sum(choice == 1) + sum(choice == 0),
      .groups = "drop"
    ) %>%
    mutate(response_rate = x1 / x2) %>%
    dplyr::select(-c(x1, x2))
  
  # Merge age
  cd_rr <- left_join(cd_rr, redcap_df[, c("participant", "age")], by = "participant")
  
  # Contrasts + agez
  cd_rr$aud_group <- as.factor(cd_rr$aud_group)
  contrasts(cd_rr$type)      <- c(-0.5, 0.5)
  contrasts(cd_rr$aud_group) <- c(-0.5, 0.5)
  cd_rr$agez <- as.numeric(scale(cd_rr$age))
  
  cd_rr
}

cd.rr.primary   <- prepare_rr(cd.task,               redcap.excluded)
cd.rr.full      <- prepare_rr(cd.task.full,           redcap.full)
cd.rr.excl_only <- prepare_rr(cd.task.taste_excl_only, redcap.taste_excl_only)

cat("N primary   :", length(unique(cd.rr.primary$participant)), "\n")
cat("N full      :", length(unique(cd.rr.full$participant)), "\n")
cat("N excl only :", length(unique(cd.rr.excl_only$participant)), "\n")


# Group split in excluded only
cd.rr.excl_only %>%
  distinct(participant, aud_group) %>%
  count(aud_group)

# -----------------------------------------------------------------------------#
# Main LMM — primary sample (reference)
# -----------------------------------------------------------------------------#
model4a_primary <- lmer(
  response_rate ~ DP * type * aud_group + agez + (DP * type | participant),
  data      = cd.rr.primary,
  na.action = na.omit,
  REML      = TRUE
)
summary(model4a_primary)
isSingular(model4a_primary)
r2_primary <- performance::r2_nakagawa(model4a_primary)
r2_primary

# Post-hoc slopes
slopes_primary <- emtrends(model4a_primary, ~ type, var = "DP")
slopes_primary
confint(slopes_primary)

# -----------------------------------------------------------------------------#
# Main LMM — full sample (sensitivity)
# -----------------------------------------------------------------------------#
model4a_full <- lmer(
  response_rate ~ DP * type * aud_group + agez + (DP * type | participant),
  data      = cd.rr.full,
  na.action = na.omit,
  REML      = TRUE
)
summary(model4a_full)
isSingular(model4a_full)
r2_full <- performance::r2_nakagawa(model4a_full)
r2_full

# Profile CIs
ci_prof_full <- confint(model4a_full, method = "profile")
ci_prof_full

slopes_full <- emtrends(model4a_full, ~ type, var = "DP")
slopes_full
confint(slopes_full)

# -----------------------------------------------------------------------------#
# Main LMM — excluded only (n=20, simplified random effects)
# NOTE: aud_group interaction unreliable (14 HC / 6 AUD)
# -----------------------------------------------------------------------------#

# Simplified model — random intercept only given small n
model4a_excl_only <- lmer(
  response_rate ~ DP * type * aud_group + agez + (1 | participant),
  data      = cd.rr.excl_only,
  na.action = na.omit,
  REML      = TRUE
)
summary(model4a_excl_only)
isSingular(model4a_excl_only)
r2_excl_only <- performance::r2_nakagawa(model4a_excl_only)
r2_excl_only

# -----------------------------------------------------------------------------#
# Global style (matches main figures)
# -----------------------------------------------------------------------------#
text_size    <- 8
lw_line      <- 0.25
lw_err       <- 0.18
pt_small     <- 0.6
panel_margin <- margin(4, 6, 4, 6)

base_panel_theme <- theme_classic(base_size = text_size) +
  theme(
    text         = element_text(family = "Arial", size = text_size),
    axis.text    = element_text(family = "Arial", size = text_size),
    axis.title   = element_text(family = "Arial", size = text_size),
    legend.text  = element_text(family = "Arial", size = text_size),
    legend.title = element_text(family = "Arial", size = text_size),
    plot.title   = element_text(family = "Arial", size = text_size, hjust = 0.5),
    plot.margin  = panel_margin
  )

col_scale <- scale_color_manual(
  values = c("#00AFBB", "#E7B800"),
  name   = "Reinforcer",
  labels = c("Juice", "Alcohol"),
  drop   = FALSE
)

# -----------------------------------------------------------------------------#
# Plot helper — HC + AUD panels per dataset
# -----------------------------------------------------------------------------#
make_rr_panels <- function(cd_rr, hc_title, aud_title,
                           y_limits = c(0.10, 0.50)) {
  
  cd_agg <- cd_rr %>%
    group_by(type, DP, aud_group) %>%
    summarise(
      mean.rr = mean(response_rate, na.rm = TRUE),
      SE      = se(response_rate),
      .groups = "drop"
    ) %>%
    mutate(type = factor(type, levels = c("Saft", "Alk")))
  
  cd_HC  <- filter(cd_agg, aud_group == 0)
  cd_AUD <- filter(cd_agg, aud_group == 1)
  
  p_HC <- ggplot(cd_HC, aes(DP, mean.rr, color = type)) +
    geom_line(linewidth = lw_line) +
    geom_point(size = pt_small) +
    geom_errorbar(aes(ymin = mean.rr - SE, ymax = mean.rr + SE),
                  width = 0.05, linewidth = lw_err) +
    ggtitle(hc_title) +
    xlab(expression(Delta * P)) + ylab("Response rate") +
    scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
    scale_y_continuous(limits = y_limits) +
    col_scale + base_panel_theme +
    theme(legend.position = "none")
  
  p_AUD <- ggplot(cd_AUD, aes(DP, mean.rr, color = type)) +
    geom_line(linewidth = lw_line) +
    geom_point(size = pt_small) +
    geom_errorbar(aes(ymin = mean.rr - SE, ymax = mean.rr + SE),
                  width = 0.05, linewidth = lw_err) +
    ggtitle(aud_title) +
    xlab(expression(Delta * P)) + ylab("Response rate") +
    scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
    scale_y_continuous(limits = y_limits) +
    col_scale + base_panel_theme +
    theme(legend.position = "none")
  
  list(HC = p_HC, AUD = p_AUD)
}

# Single panel for excluded only (no group split)
make_rr_single <- function(cd_rr, title, y_limits = c(0.10, 0.60)) {
  
  cd_agg <- cd_rr %>%
    group_by(type, DP) %>%
    summarise(
      mean.rr = mean(response_rate, na.rm = TRUE),
      SE      = se(response_rate),
      .groups = "drop"
    ) %>%
    mutate(type = factor(type, levels = c("Saft", "Alk")))
  
  ggplot(cd_agg, aes(DP, mean.rr, color = type)) +
    geom_line(linewidth = lw_line) +
    geom_point(size = pt_small) +
    geom_errorbar(aes(ymin = mean.rr - SE, ymax = mean.rr + SE),
                  width = 0.05, linewidth = lw_err) +
    ggtitle(title) +
    xlab(expression(Delta * P)) + ylab("Response rate") +
    scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
    scale_y_continuous(limits = y_limits) +
    col_scale + base_panel_theme +
    theme(legend.position = "bottom")
}

# -----------------------------------------------------------------------------#
# Generate plots
# -----------------------------------------------------------------------------#
n_hc_primary <- sum(cd.rr.primary$aud_group == 0 &
                      !duplicated(cd.rr.primary$participant))
n_aud_primary <- sum(cd.rr.primary$aud_group == 1 &
                       !duplicated(cd.rr.primary$participant))
n_hc_full <- sum(cd.rr.full$aud_group == 0 &
                   !duplicated(cd.rr.full$participant))
n_aud_full <- sum(cd.rr.full$aud_group == 1 &
                    !duplicated(cd.rr.full$participant))

panels_primary <- make_rr_panels(
  cd.rr.primary,
  paste0("HC — Primary (n = ", n_hc_primary, ")"),
  paste0("AUD — Primary (n = ", n_aud_primary, ")")
)

panels_full <- make_rr_panels(
  cd.rr.full,
  paste0("HC — Full (n = ", n_hc_full, ")"),
  paste0("AUD — Full (n = ", n_aud_full, ")")
)

panel_excl_only <- make_rr_single(
  cd.rr.excl_only,
  paste0("Taste-excluded only (n = 20; 14 HC / 6 AUD)")
)

# Shared legend
legend_plot <- ggplot(
  data.frame(type = c("Saft", "Alk"), y = c(1, 1), x = c(1, 1)),
  aes(x, y, color = type)
) +
  geom_line() +
  col_scale +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = text_size))

shared_legend <- ggpubr::get_legend(legend_plot)

# Assemble rows
row1 <- ggarrange(panels_primary$HC, panels_primary$AUD,
                  ncol = 2, legend = "none")
row2 <- ggarrange(panels_full$HC, panels_full$AUD,
                  ncol = 2, legend = "none")
row3 <- ggarrange(panel_excl_only, ncol = 1, legend = "bottom")

fig_sensitivity_rr <- ggarrange(
  row1, row2, row3,
  nrow   = 3,
  labels = c("Primary sample", "Full sample", "Taste-excluded only"),
  font.label = list(size = text_size, face = "bold", family = "Arial")
)

fig_sensitivity_rr

# -----------------------------------------------------------------------------#
# Save
# -----------------------------------------------------------------------------#
ggsave(
  filename = file.path(fig_dir, "FigureS_sensitivity_RR.png"),
  plot     = fig_sensitivity_rr,
  width    = 180, height = 240, units = "mm",
  dpi      = 800, device = ragg::agg_png
)

ggsave(
  filename = file.path(fig_dir, "FigureS_sensitivity_RR.tiff"),
  plot     = fig_sensitivity_rr,
  width    = 180, height = 240, units = "mm",
  dpi      = 800, device = ragg::agg_tiff,
  compression = "lzw"
)

ggsave(
  filename = file.path(fig_dir, "FigureS_sensitivity_RR.pdf"),
  plot     = fig_sensitivity_rr,
  width    = 180, height = 240, units = "mm",
  device   = cairo_pdf
)

cat("\nDone! Saved sensitivity response rate figures.\n")

# -----------------------------------------------------------------------------#
# Causality sensitivity analysis
# -----------------------------------------------------------------------------#

load(file.path(EXCL_DIR, "cd.causality.RData"))
load(file.path(EXCL_DIR, "cd.causality.full.RData"))
load(file.path(EXCL_DIR, "cd.causality.taste_excl_only.RData"))
load(file.path(EXCL_DIR, "df_pleasantness.RData"))

cd.causality$participant                 <- as.numeric(as.character(cd.causality$participant))
cd.causality.full$participant            <- as.numeric(as.character(cd.causality.full$participant))
cd.causality.taste_excl_only$participant <- as.numeric(as.character(cd.causality.taste_excl_only$participant))
df_pleasantness$participant              <- as.numeric(as.character(df_pleasantness$participant))

df_pleasantness <- df_pleasantness %>%
  distinct(participant, type, mean_pleasantness) %>%
  mutate(mean_pleasantness_z = as.numeric(scale(mean_pleasantness)))

# -----------------------------------------------------------------------------#
# Primary sample
# -----------------------------------------------------------------------------#
cd.causality <- cd.causality %>%
  arrange(participant, type, DP, dummycausal)

df.caus4.primary <- dplyr::filter(cd.causality, dummycausal == 1) %>%
  arrange(participant, type, DP)
df.caus4.primary$dummycausal <- 4
df.caus4.primary$rating <-
  dplyr::filter(cd.causality, dummycausal == 2) %>%
  arrange(participant, type, DP) %>% pull(rating) -
  dplyr::filter(cd.causality, dummycausal == 3) %>%
  arrange(participant, type, DP) %>% pull(rating)

df.caus4.primary$agez <- as.numeric(scale(df.caus4.primary$age))
df.caus4.primary <- df.caus4.primary %>%
  left_join(
    df_pleasantness %>% dplyr::select(participant, type, mean_pleasantness, mean_pleasantness_z),
    by = c("participant", "type")
  )
df.caus4.primary$type      <- factor(df.caus4.primary$type, levels = c("Saft", "Alk"))
df.caus4.primary$aud_group <- as.factor(df.caus4.primary$aud_group)
contrasts(df.caus4.primary$type)      <- c(-0.5, 0.5)
contrasts(df.caus4.primary$aud_group) <- c(-0.5, 0.5)
cat("N causality primary:", length(unique(df.caus4.primary$participant)), "\n")

# -----------------------------------------------------------------------------#
# Full sample
# -----------------------------------------------------------------------------#
cd.causality.full <- cd.causality.full %>%
  arrange(participant, type, DP, dummycausal)

df.caus4.full <- dplyr::filter(cd.causality.full, dummycausal == 1) %>%
  arrange(participant, type, DP)
df.caus4.full$dummycausal <- 4
df.caus4.full$rating <-
  dplyr::filter(cd.causality.full, dummycausal == 2) %>%
  arrange(participant, type, DP) %>% pull(rating) -
  dplyr::filter(cd.causality.full, dummycausal == 3) %>%
  arrange(participant, type, DP) %>% pull(rating)

df.caus4.full$agez      <- as.numeric(scale(df.caus4.full$age))
df.caus4.full$type      <- factor(df.caus4.full$type, levels = c("Saft", "Alk"))
df.caus4.full$aud_group <- as.factor(df.caus4.full$aud_group)
contrasts(df.caus4.full$type)      <- c(-0.5, 0.5)
contrasts(df.caus4.full$aud_group) <- c(-0.5, 0.5)
cat("N causality full:", length(unique(df.caus4.full$participant)), "\n")

# -----------------------------------------------------------------------------#
# Excluded only
# -----------------------------------------------------------------------------#
cd.causality.taste_excl_only <- cd.causality.taste_excl_only %>%
  arrange(participant, type, DP, dummycausal)

df.caus4.excl_only <- dplyr::filter(cd.causality.taste_excl_only, dummycausal == 1) %>%
  arrange(participant, type, DP)
df.caus4.excl_only$dummycausal <- 4
df.caus4.excl_only$rating <-
  dplyr::filter(cd.causality.taste_excl_only, dummycausal == 2) %>%
  arrange(participant, type, DP) %>% pull(rating) -
  dplyr::filter(cd.causality.taste_excl_only, dummycausal == 3) %>%
  arrange(participant, type, DP) %>% pull(rating)

df.caus4.excl_only$agez      <- as.numeric(scale(df.caus4.excl_only$age))
df.caus4.excl_only$type      <- factor(df.caus4.excl_only$type, levels = c("Saft", "Alk"))
df.caus4.excl_only$aud_group <- as.factor(df.caus4.excl_only$aud_group)
contrasts(df.caus4.excl_only$type)      <- c(-0.5, 0.5)
contrasts(df.caus4.excl_only$aud_group) <- c(-0.5, 0.5)
cat("N causality excl only:", length(unique(df.caus4.excl_only$participant)), "\n")

# -----------------------------------------------------------------------------#
# Models
# -----------------------------------------------------------------------------#
model1g_primary <- lmer(
  rating ~ DP * type * aud_group + agez + (1 | participant),
  data = df.caus4.primary, na.action = na.omit
)
summary(model1g_primary)
isSingular(model1g_primary)
r2_caus_primary <- performance::r2_nakagawa(model1g_primary)
r2_caus_primary
slopes_caus_primary <- emtrends(model1g_primary, ~ aud_group, var = "DP")
slopes_caus_primary
confint(slopes_caus_primary)

model1g_full <- lmer(
  rating ~ DP * type * aud_group + agez + (1 | participant),
  data = df.caus4.full, na.action = na.omit
)
summary(model1g_full)
isSingular(model1g_full)

# Profile CIs
ci_prof_csfull <- confint(model1g_full, method = "profile")
ci_prof_csfull

r2_caus_full <- performance::r2_nakagawa(model1g_full)
r2_caus_full
slopes_caus_full <- emtrends(model1g_full, ~ aud_group, var = "DP")
slopes_caus_full
confint(slopes_caus_full)

# Excluded only — descriptives + simplified model (n=20, underpowered)
df.caus4.excl_only %>%
  group_by(type, DP) %>%
  summarise(
    mean_rating = mean(rating, na.rm = TRUE),
    sd_rating   = sd(rating, na.rm = TRUE),
    n           = n(),
    .groups = "drop"
  )

model1g_excl_only <- lmer(
  rating ~ DP * type * aud_group + agez + (1 | participant),
  data = df.caus4.excl_only, na.action = na.omit
)
summary(model1g_excl_only)
isSingular(model1g_excl_only)
r2_caus_excl_only <- performance::r2_nakagawa(model1g_excl_only)
r2_caus_excl_only

# -----------------------------------------------------------------------------#
# Causality plots
# -----------------------------------------------------------------------------#
make_caus_panels <- function(df4, hc_title, aud_title,
                             y_limits = c(-0.9, 0.8)) {
  caus_agg <- df4 %>%
    group_by(type, DP, aud_group) %>%
    summarise(
      m.rating = mean(rating, na.rm = TRUE),
      SE       = se(rating),
      .groups  = "drop"
    )
  
  caus_HC  <- filter(caus_agg, aud_group == 0)
  caus_AUD <- filter(caus_agg, aud_group == 1)
  
  p_HC <- ggplot(caus_HC, aes(DP, m.rating, color = type, group = type)) +
    geom_line(linewidth = lw_line) +
    geom_point(size = pt_small) +
    geom_errorbar(aes(ymin = m.rating - SE, ymax = m.rating + SE),
                  width = 0.05, linewidth = lw_err) +
    ggtitle(hc_title) +
    xlab(expression(Delta * P)) + ylab("Causality judgements") +
    scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
    scale_y_continuous(limits = y_limits) +
    col_scale + base_panel_theme +
    theme(legend.position = "none")
  
  p_AUD <- ggplot(caus_AUD, aes(DP, m.rating, color = type, group = type)) +
    geom_line(linewidth = lw_line) +
    geom_point(size = pt_small) +
    geom_errorbar(aes(ymin = m.rating - SE, ymax = m.rating + SE),
                  width = 0.05, linewidth = lw_err) +
    ggtitle(aud_title) +
    xlab(expression(Delta * P)) + ylab("Causality judgements") +
    scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
    scale_y_continuous(limits = y_limits) +
    col_scale + base_panel_theme +
    theme(legend.position = "none")
  
  list(HC = p_HC, AUD = p_AUD)
}

make_caus_single <- function(df4, title, y_limits = c(-0.9, 0.8)) {
  caus_agg <- df4 %>%
    group_by(type, DP) %>%
    summarise(
      m.rating = mean(rating, na.rm = TRUE),
      SE       = se(rating),
      .groups  = "drop"
    )
  
  ggplot(caus_agg, aes(DP, m.rating, color = type, group = type)) +
    geom_line(linewidth = lw_line) +
    geom_point(size = pt_small) +
    geom_errorbar(aes(ymin = m.rating - SE, ymax = m.rating + SE),
                  width = 0.05, linewidth = lw_err) +
    ggtitle(title) +
    xlab(expression(Delta * P)) + ylab("Causality judgements") +
    scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
    scale_y_continuous(limits = y_limits) +
    col_scale + base_panel_theme +
    theme(legend.position = "bottom")
}

n_hc_primary_c  <- sum(df.caus4.primary$aud_group == 0 &
                         !duplicated(df.caus4.primary$participant))
n_aud_primary_c <- sum(df.caus4.primary$aud_group == 1 &
                         !duplicated(df.caus4.primary$participant))
n_hc_full_c     <- sum(df.caus4.full$aud_group == 0 &
                         !duplicated(df.caus4.full$participant))
n_aud_full_c    <- sum(df.caus4.full$aud_group == 1 &
                         !duplicated(df.caus4.full$participant))

panels_caus_primary <- make_caus_panels(
  df.caus4.primary,
  paste0("HC — Primary (n = ", n_hc_primary_c, ")"),
  paste0("AUD — Primary (n = ", n_aud_primary_c, ")")
)

panels_caus_full <- make_caus_panels(
  df.caus4.full,
  paste0("HC — Full (n = ", n_hc_full_c, ")"),
  paste0("AUD — Full (n = ", n_aud_full_c, ")")
)

panel_caus_excl_only <- make_caus_single(
  df.caus4.excl_only,
  "Taste-excluded only (n = 20; 14 HC / 6 AUD)"
)

row1_c <- ggarrange(panels_caus_primary$HC, panels_caus_primary$AUD,
                    ncol = 2, legend = "none")
row2_c <- ggarrange(panels_caus_full$HC, panels_caus_full$AUD,
                    ncol = 2, legend = "none")
row3_c <- ggarrange(panel_caus_excl_only, ncol = 1, legend = "bottom")

fig_sensitivity_caus <- ggarrange(
  row1_c, row2_c, row3_c,
  nrow   = 3,
  labels = c("Primary sample", "Full sample", "Taste-excluded only"),
  font.label = list(size = text_size, face = "bold", family = "Arial")
)

fig_sensitivity_caus

# Save
ggsave(
  filename = file.path(fig_dir, "FigureS_sensitivity_causality.png"),
  plot     = fig_sensitivity_caus,
  width    = 180, height = 240, units = "mm",
  dpi      = 800, device = ragg::agg_png
)

ggsave(
  filename = file.path(fig_dir, "FigureS_sensitivity_causality.tiff"),
  plot     = fig_sensitivity_caus,
  width    = 180, height = 240, units = "mm",
  dpi      = 800, device = ragg::agg_tiff,
  compression = "lzw"
)

ggsave(
  filename = file.path(fig_dir, "FigureS_sensitivity_causality.pdf"),
  plot     = fig_sensitivity_caus,
  width    = 180, height = 240, units = "mm",
  device   = cairo_pdf
)

cat("\nDone! Saved sensitivity causality figures.\n")

# Excluded only — response rate by DP and reinforcer type, split by group
panels_excl_only <- make_rr_panels(
  cd.rr.excl_only,
  "HC — Taste-excluded (n = 14)",
  "AUD — Taste-excluded (n = 6)",
  y_limits = c(0.10, 0.60)
)

panel_excl_only_grouped <- ggarrange(
  panels_excl_only$HC, panels_excl_only$AUD,
  ncol = 2, legend = "bottom",
  common.legend = TRUE
)

panel_excl_only_grouped

# Save
ggsave(
  filename = file.path(fig_dir, "FigureS_excl_only_RR_grouped.png"),
  plot     = panel_excl_only_grouped,
  width    = 180, height = 90, units = "mm",
  dpi      = 800, device = ragg::agg_png
)