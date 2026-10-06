# ConDeg_B01: Contingency Degradation in AUD — Analysis Code
Analysis Code for the Psychological Medicine Paper "Alcohol reinforcers attenuate goal-directed control in individuals with and without alcohol use disorder"

Data availability

The data are not included because of data-protection restrictions on human-subjects data. They are available from the corresponding author on reasonable request, subject to institutional and ethics approval. None of the scripts run without the data.

Requirements

R (version 4.3.3) with: tidyverse, lme4, lmerTest, emmeans, performance, DHARMa, influence.ME, robustlmm, car, rstatix, ez, coin, moments, Hmisc, irr, lubridate, ggplot2, cowplot, ggpubr, patchwork, ggdist, raincloudplots, viridis, ragg.

The repo assumes no project structure. Before running a script, set the input/output paths at its top to your local copy of the data.

Scripts (run in order)
1_Exclusion_Script_Wp2_LB.R: Harmonizes IDs, timepoints and gender coding, then applies exclusions (baseline taste rating < 0.5 or |diff_AS| > 0.4; missing data). Output: three samples in R_dfs/final_RData_excluded/: main (all exclusions), full (missingness only), and taste-excluded only.
2_Taste_Rating_Bode_ConDeg_Paper.R: Demographics and group comparisons, plus the LMM on taste pleasantness. Requires pss_sum to be merged into redcap.excluded first. Output: df_pleasantness.RData, used as a covariate in script 3.
3_Response_Rates&Causality_Ratings_Analysis_Bode_ConDeg_Paper.R: Main LMMs for response rates and causality ratings, with diagnostics, robustness and sensitivity checks. Also covers the ratio scores and the supplementary CJ analysis. Output: model objects and robust p-value/CI tables.
4_Reliability_Analysis_ConDeg.R: Permuted odd-even split-half reliability of individual ΔP slopes (5,000 iterations, Spearman–Brown corrected).
5_Response_Rates&Causality_Ratings_AnalysisFULLSAMPLENOEXCLUSIONS.R: Supplementary comparison of the main, full and taste-excluded samples. Output: supplementary figures.
Psychological_Medicine_Figure_Script.R: Generates the manuscript figures. Run scripts 2 and 3 in the same R session beforehand, because this script uses objects they create in the workspace.

Conventions
rating means different things by script. In script 2 it is taste pleasantness. In scripts 3 and 5 it is the causality rating, rating(dummycausal==2) − rating(dummycausal==3), modeled on the raw 0–1 scale.
ΔP vs. CJ. ΔP is the programmed contingency. CJ/CJz is the participant's explicit contingency judgment, P(O|A) − P(O|¬A), which ranges from −1 to 1.

License
MIT — see LICENSE.
