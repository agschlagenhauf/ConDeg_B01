#### Figures for Psychological Medicine
# Arial, font size 8
# Author: Erik Lukas Bode (erik-lukas.bode@charite.de) 
#########################################################################
################## Con Deg Paper Figures 2025 LB ########################
#########################################################################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggh4x)

# Devices for reliable font rendering
library(ragg)     # for PNG/TIFF
# cairo_pdf is available via grDevices, no extra package needed

## ---------- GLOBAL STYLE (shared by all figures) ----------

text_size <- 8          # axis, titles, legend -> Arial 8 as required
base_size <- text_size  # theme_classic base
strip_size <- text_size # facet strip text size

lw_line   <- 0.25       # main lines
lw_err    <- 0.18       # error bars
pt_main   <- 1.8        # bigger points (pleasantness)
pt_small  <- 0.6        # smaller points (ΔP / Fig2)
pt_supp   <- 1.2        # points in supplement figure (Response x Causality)

panel_margin <- margin(4, 6, 4, 6)

base_panel_theme <- theme_classic(base_size = base_size) +
  theme(
    text         = element_text(family = "Arial", size = text_size),
    axis.text    = element_text(family = "Arial", size = text_size),
    axis.title   = element_text(family = "Arial", size = text_size),
    legend.text  = element_text(family = "Arial", size = text_size),
    legend.title = element_text(family = "Arial", size = text_size),
    plot.title   = element_text(family = "Arial", size = text_size, hjust = 0.5),
    plot.margin  = panel_margin,
    # facet styling (used in supplement)
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", size = strip_size, face = "bold"),
    strip.placement = "outside"
  )

## ----------------------------------------------------------------------
## Figure 1: Pleasantness
## ----------------------------------------------------------------------

# Standard error helper (used in plots / summaries)
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

options(digits = 10)

taste.ratings.agg.excluded <- taste.ratings.excluded %>%
  group_by(aud_group, type, number) %>%
  summarise(
    mean.r = mean(rating, na.rm = TRUE),
    SE.r   = se(rating),
    .groups = "drop"
  )

taste_plot_HC  <- filter(taste.ratings.agg.excluded, aud_group == 0)
taste_plot_AUD <- filter(taste.ratings.agg.excluded, aud_group == 1)

x_limits <- c(0, 6.5)
y_limits <- c(0.4, 1.0)

pleasplot_HC <- ggplot(
  taste_plot_HC,
  aes(x = number, y = mean.r, color = type)
) +
  geom_line(linewidth = lw_line) +
  geom_point(size = pt_main) +
  geom_errorbar(
    aes(ymin = mean.r - SE.r, ymax = mean.r + SE.r),
    width = 0.1, linewidth = lw_err
  ) +
  ggtitle("HC") +
  xlab("Timepoint") +
  ylab("Pleasantness rating") +
  scale_x_continuous(
    breaks = 0:6, limits = x_limits,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = y_limits, breaks = seq(0.4, 1.0, 0.1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_color_manual(
    values = c("#00AFBB", "#E7B800"),
    name   = "Reinforcer",
    labels = c("Juice", "Alcohol")
  ) +
  base_panel_theme

pleasplot_AUD <- ggplot(
  taste_plot_AUD,
  aes(x = number, y = mean.r, color = type)
) +
  geom_line(linewidth = lw_line) +
  geom_point(size = pt_main) +
  geom_errorbar(
    aes(ymin = mean.r - SE.r, ymax = mean.r + SE.r),
    width = 0.1, linewidth = lw_err
  ) +
  ggtitle("AUD") +
  xlab("Timepoint") +
  ylab("Pleasantness rating") +
  scale_x_continuous(
    breaks = 0:6, limits = x_limits,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = y_limits, breaks = seq(0.4, 1.0, 0.1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_color_manual(
    values = c("#00AFBB", "#E7B800"),
    name   = "Reinforcer",
    labels = c("Juice", "Alcohol")
  ) +
  base_panel_theme

combined_plot <- ggarrange(
  pleasplot_HC, pleasplot_AUD,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)

combined_plot

## ----------------------------------------------------------------------
## Save outputs (use ragg/cairo for stable Arial rendering)
## ----------------------------------------------------------------------

# PNG (recommended for quick checks)
ggsave(
  filename = "pleasantness_group_comparison.png",
  plot     = combined_plot,
  width    = 180, height = 120, units = "mm",
  dpi      = 800,
  device   = ragg::agg_png
)

# TIFF (journal-friendly)
ggsave(
  filename = "pleasantness_group_comparison.tiff",
  plot     = combined_plot,
  width    = 180, height = 120, units = "mm",
  dpi      = 800,
  device   = ragg::agg_tiff,
  compression = "lzw"
)

# PDF (vector, best for print)
ggsave(
  filename = "pleasantness_group_comparison.pdf",
  plot     = combined_plot,
  width    = 180, height = 120, units = "mm",
  device   = cairo_pdf
)

## ----------------------------------------------------------------------
## Figure 2: A) Response (ΔP), B) Causality (ΔP)  [A+B only]
## Robust legend: NO get_legend(), let ggarrange collect it (common.legend)
## ----------------------------------------------------------------------

# Make sure type is a factor with stable order (adjust labels if needed)
cd.rr <- cd.rr %>%
  mutate(type = factor(type, levels = c("Saft", "Alk")))

df.causality4 <- df.causality4 %>%
  mutate(type = factor(type, levels = c("Saft", "Alk")))

# Aggregate response rate
cd.rr.agg <- cd.rr %>%
  group_by(type, DP, aud_group) %>%
  summarise(
    mean.rr = mean(response_rate, na.rm = TRUE),
    SE      = se(response_rate),
    .groups = "drop"
  )

cd_HC  <- cd.rr.agg %>% filter(aud_group == 0)
cd_AUD <- cd.rr.agg %>% filter(aud_group == 1)

# Aggregate causality judgments
causality.rating.agg <- df.causality4 %>%
  group_by(type, DP, dummycausal, aud_group) %>%
  summarise(
    m.rating = mean(rating, na.rm = TRUE),
    SE       = se(rating),
    .groups  = "drop"
  )

caus_HC  <- causality.rating.agg %>% filter(aud_group == 0)
caus_AUD <- causality.rating.agg %>% filter(aud_group == 1)

y_resp_limits <- c(0.10, 0.50)

# Column titles (stable pt->mm conversion for annotate size)
title_size_mm <- text_size * 0.352778

title_HC <- ggplot() +
  annotate(
    "text", x = 0.5, y = 0.5, label = "HC",
    fontface = "bold", size = title_size_mm, family = "Arial"
  ) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

title_AUD <- ggplot() +
  annotate(
    "text", x = 0.5, y = 0.5, label = "AUD",
    fontface = "bold", size = title_size_mm, family = "Arial"
  ) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

col_titles <- ggarrange(title_HC, title_AUD, ncol = 2)

# Common color scale (keeps legend consistent, prevents dropping)
col_scale <- scale_color_manual(
  values = c("#00AFBB", "#E7B800"),
  name   = "Reinforcer",
  labels = c("Juice", "Alcohol"),
  drop   = FALSE
)

# ---------- A) Response rate ----------
p_a_HC <- ggplot(cd_HC, aes(DP, mean.rr, color = type)) +
  geom_line(linewidth = lw_line) +
  geom_point(size = pt_small) +
  geom_errorbar(
    aes(ymin = mean.rr - SE, ymax = mean.rr + SE),
    width = 0.05, linewidth = lw_err
  ) +
  xlab(expression(Delta * P)) +
  ylab("Response rate") +
  scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  scale_y_continuous(limits = y_resp_limits) +
  col_scale +
  base_panel_theme +
  theme(legend.position = "none")

p_a_AUD <- (p_a_HC %+% cd_AUD) + theme(legend.position = "none")
row_A <- ggarrange(p_a_HC, p_a_AUD, ncol = 2)

# ---------- B) Causality judgements ----------
p_b_HC <- ggplot(caus_HC, aes(DP, m.rating, color = type)) +
  geom_line(linewidth = lw_line) +
  geom_point(size = pt_small) +
  geom_errorbar(
    aes(ymin = m.rating - SE, ymax = m.rating + SE),
    width = 0.05, linewidth = lw_err
  ) +
  xlab(expression(Delta * P)) +
  ylab("Causality judgements") +
  scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  col_scale +
  base_panel_theme +
  theme(legend.position = "none")

# AUD panel: keep legend ON here (bottom), and let ggarrange collect it
p_b_AUD <- (p_b_HC %+% caus_AUD) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(
    override.aes = list(linewidth = lw_line, size = pt_small)
  ))

row_B <- ggarrange(
  p_b_HC, p_b_AUD,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)

# Assemble Figure A and Figure B (legend already included in row_B)
Figure_A <- ggarrange(col_titles, row_A, ncol = 1, heights = c(0.12, 1))
Figure_B <- ggarrange(col_titles, row_B, ncol = 1, heights = c(0.12, 1.20))

Figure_AB <- ggarrange(
  Figure_A, Figure_B,
  ncol = 1,
  labels = c("A", "B"),
  font.label = list(size = text_size, face = "bold", family = "Arial")
)

Figure_AB

# -------------------------
# Save
# -------------------------
fig_dir <- file.path(ROOT, "Figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

message("Saving figures to: ", normalizePath(fig_dir))

ggsave(
  filename = file.path(fig_dir, "Figure_AB_combined.tiff"),
  plot     = Figure_AB,
  width    = 180,
  height   = 120,
  units    = "mm",
  dpi      = 800,
  device   = ragg::agg_tiff,
  compression = "lzw"
)

ggsave(
  filename = file.path(fig_dir, "Figure_AB_combined.png"),
  plot     = Figure_AB,
  width    = 180,
  height   = 120,
  units    = "mm",
  dpi      = 800,
  device   = ragg::agg_png
)

ggsave(
  filename = file.path(fig_dir, "Figure_AB_combined.pdf"),
  plot     = Figure_AB,
  width    = 180,
  height   = 120,
  units    = "mm",
  device   = cairo_pdf
)

## ----------------------------------------------------------------------
## Supplement Figure: Response x Causality
## Panel A: facet by group (color/shape = reinforcer)
## Panel B: facet by reinforcer (color/shape = group)
## ----------------------------------------------------------------------
merged_plot_data <- left_join(
  causality.rating.agg,
  cd.rr.agg,
  by = c("DP", "type", "aud_group")
)

merged_plot_data <- merged_plot_data %>%
  mutate(
    aud_group = factor(aud_group, levels = c("HC","AUD")),
    type      = factor(type, levels = c("Juice","Alcohol"))
  )

y_lim <- range(
  merged_plot_data$mean.rr - merged_plot_data$SE.y,
  merged_plot_data$mean.rr + merged_plot_data$SE.y,
  na.rm = TRUE
)

x_lim <- range(
  merged_plot_data$m.rating - merged_plot_data$SE.x,
  merged_plot_data$m.rating + merged_plot_data$SE.x,
  na.rm = TRUE
)

plot_A_facet_group_bothSE <- ggplot(
  merged_plot_data,
  aes(x = m.rating, y = mean.rr, color = type, shape = type)
) +
  geom_errorbarh(
    aes(xmin = m.rating - SE.x, xmax = m.rating + SE.x),
    height = 0, linewidth = lw_err, alpha = 0.8
  ) +
  geom_errorbar(
    aes(ymin = mean.rr - SE.y, ymax = mean.rr + SE.y),
    width = 0, linewidth = lw_err, alpha = 0.8
  ) +
  geom_line(aes(group = type), linewidth = lw_line, alpha = 0.7) +
  geom_point(size = pt_supp) +
  ggh4x::facet_wrap2(~ aud_group, nrow = 1, axes = "all") +
  scale_color_manual(
    values = c("Juice" = "#00AFBB", "Alcohol" = "#E7B800"),
    name = "Reinforcer"
  ) +
  scale_shape_manual(
    values = c("Juice" = 16, "Alcohol" = 17),
    name = "Reinforcer"
  ) +
  labs(
    x = "Causality judgements",
    y = "Response rate (responses/sec)"
  ) +
  scale_x_continuous(limits = x_lim) +
  scale_y_continuous(limits = y_lim) +
  base_panel_theme +
  theme(legend.box = "horizontal")

plot_B_facet_reinforcer_bothSE <- ggplot(
  merged_plot_data,
  aes(x = m.rating, y = mean.rr, color = aud_group, shape = aud_group)
) +
  geom_errorbarh(
    aes(xmin = m.rating - SE.x, xmax = m.rating + SE.x),
    height = 0, linewidth = lw_err, alpha = 0.8
  ) +
  geom_errorbar(
    aes(ymin = mean.rr - SE.y, ymax = mean.rr + SE.y),
    width = 0, linewidth = lw_err, alpha = 0.8
  ) +
  geom_line(aes(group = aud_group), linewidth = lw_line, alpha = 0.7) +
  geom_point(size = pt_supp) +
  ggh4x::facet_wrap2(~ type, nrow = 1, axes = "all") +
  scale_color_manual(
    values = c("HC" = "#0072B2", "AUD" = "#D55E00"),
    name = "Group"
  ) +
  scale_shape_manual(
    values = c("HC" = 16, "AUD" = 17),
    name = "Group"
  ) +
  labs(
    x = "Causality judgements",
    y = "Response rate (responses/sec)"
  ) +
  scale_x_continuous(limits = x_lim) +
  scale_y_continuous(limits = y_lim) +
  base_panel_theme +
  theme(legend.box = "horizontal")

# -------------------------------------------------
# Legends under EACH panel (A gets Reinforcer, B gets Group)
# -------------------------------------------------

plot_A_noleg <- plot_A_facet_group_bothSE + theme(legend.position = "none")
plot_B_noleg <- plot_B_facet_reinforcer_bothSE + theme(legend.position = "none")

legend_reinf <- get_legend(
  plot_A_facet_group_bothSE +
    theme(legend.position = "bottom", legend.box = "horizontal",
          text = element_text(family = "Arial"))
)

legend_group <- get_legend(
  plot_B_facet_reinforcer_bothSE +
    theme(legend.position = "bottom", legend.box = "horizontal",
          text = element_text(family = "Arial"))
)

panel_A_block <- ggarrange(
  plot_A_noleg, legend_reinf,
  ncol = 1,
  heights = c(1, 0.18)
)

panel_B_block <- ggarrange(
  plot_B_noleg, legend_group,
  ncol = 1,
  heights = c(1, 0.18)
)

panel_A_block <- panel_A_block + theme(plot.margin = margin(5, 5, 5, 15))
panel_B_block <- panel_B_block + theme(plot.margin = margin(5, 5, 5, 15))

Figure_C_supp <- ggarrange(
  panel_A_block,
  panel_B_block,
  ncol = 1,
  labels = c("A", "B"),
  font.label = list(size = text_size, face = "bold", family = "Arial"),
  heights = c(1, 1),
  align = "v",
  label.x = -0.001,
  label.y = 0.99
)
Figure_C_supp

ggsave(
  "Figure_C_supp_AB.tiff",
  Figure_C_supp,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "Figure_C_supp_AB.pdf",
  Figure_C_supp,
  width = 180,
  height = 120,
  units = "mm",
  device = cairo_pdf
)
