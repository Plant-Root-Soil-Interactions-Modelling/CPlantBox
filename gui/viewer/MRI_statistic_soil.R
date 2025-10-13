#Soil
library("datarium")
library(tidyverse)
library(ggpubr)
library(rstatix)

Plant <- c("5", "6", "7", "8", "5", "6", "7", "8", "5", "6", "7", "8")
Annotation <- c("M", "M", "M", "M", "Mp", "Mp", "Mp", "Mp", "A",  "A",  "A",  "A")
TRL <- c(269.24,	488.34,	68.56,	76.77,	293.15,	480.69,	76.02,	75.61,	276.57,	473.02,	64.37,	68.02)
RR <-c(83.36,	96.13,	80.66,	91.39,	90.76,	94.62,	89.44,	90.01,	85.63,	93.11,	75.73,	80.98)
RLD <-c(0.59,	1.07,	0.15,	0.17,	0.64,	1.06,	0.17,	0.17,	0.61,	1.04,	0.14,	0.15)
HMD <-c(0.73,	0.54,	1.45,	1.37,	0.70,	0.55,	1.38,	1.38,	0.72,	0.55,	1.50,	1.46)
Tips <-c(75,	192,	35,	15,	109,	199,	39,	15,	124,	267,	34,	17)
LTipsone <-c(60, 45, 34, 13, 59, 44, 37, 12, 37, 39, 11, 13)
LTipstwo <-c(14, 137, 0, 1, 49, 146, 1, 2, 64, 127, 18, 3)
LTipsthree <-c(0, 9, 0, 0, 0, 8, 0, 0, 20, 73, 1, 0)
mean_radius <-c( 0.0231,	0.0223,	0.0288,	0.0305,	0.0228,	0.023,	0.0332,	0.0337,	0.0243,	0.022,	0.0284,	0.0312)
Krs_c <-c( 0.00318,	0.00603,	0.0019,	0.00223,	0.00342,	0.0057,	0.00228,	0.00242,	0.0033,	0.0055,	0.00177,	0.00205)
Krs_v <-c(0.0131,	0.0163,	0.0133,	0.0109,	0.0142,	0.0161,	0.0139,	0.0115,	0.0125,	0.0144,	0.00865,	0.0104)
Network_conductance_c <-c( 1.1811E-05,	1.2348E-05,	2.7713E-05,	2.90478E-05,	1.16664E-05,	1.1858E-05,	2.99909E-05,	3.20063E-05,	1.19319E-05,	1.16274E-05,	2.74973E-05,	3.01382E-05)
Network_Conductance_v <-c(4.86555E-05,	3.33784E-05,	0.000193991,	0.000141983,	4.84394E-05,	3.34935E-05,	0.000182839,	0.000152096,	4.51965E-05,	3.04427E-05,	0.000134379,	0.000152896)
zSUFc <-c(-9.28,	-9.726,	-3.35,	-3.466,	-9.036,	-9.531,	-3.26,	-3.455,	-9.356,	-9.751,	-3.075,	-3.618)
zSUFv <-c(-6.64,	-5.378,	-2.964,	-2.171,	-6.303,	-5.021,	-2.916,	-2.156,	-6.044,	-5.813,	-1.748,	-2.59)

Soil <- data.frame(Plant, Annotation, TRL, RR, RLD, HMD, Tips, LTipsone, LTipstwo, LTipsthree, mean_radius, Krs_c, Krs_v, Network_conductance_c, Network_conductance_v, zSUFc, zSUFv)
Soil_ <-as_tibble(Soil)
Soil_ <- Soil_ %>%
  convert_as_factor(Plant, Annotation)
print (Soil_)

# only for M and M+
Plant_time <- c("5", "6", "7", "8", "5", "6", "7", "8")
Annotation_time <- c("M", "M", "M", "M", "Mp", "Mp", "Mp", "Mp")
Annotation_rate <- c(7.67,	6.20,	4.30,	5.02, 9.37,	8.61,	5.06,	6.52)
# RATES SOIL
#Annotation_rate <-c(7.67, 6.20, 4.30, 5.02, 9.37, 8.61, 5.06, 6.52)
Soil_time <- data.frame(Plant_time, Annotation_time, Annotation_rate)
Soil_time <-as_tibble(Soil_time)
Soil_time <- Soil_time %>%
  convert_as_factor(Plant_time, Annotation_time)
print (Soil_time)

# TRL
bxp1 <- ggboxplot(Soil_, x = "Annotation", y = "TRL", add = "point", title="TRL Soil")
bxp1

Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(TRL)

Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(TRL, type = "mean_sd")

Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(TRL)



ggqqplot(Soil_, "TRL", facet.by = "Annotation", title="TRL Soil") #graue Balken sind 95% Konfidenzintervall

res.aov <- anova_test(data = Soil_, dv = TRL, wid = Plant, within = Annotation, detailed = TRUE)
get_anova_table(res.aov)

pwc_TRL <- Soil_ %>%
  pairwise_t_test(
    TRL ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided",
    detailed = TRUE
  )
pwc_TRL


pwc_TRL <- pwc_TRL %>% add_xy_position(x = "Annotation")
bxp1 + 
  stat_pvalue_manual(pwc_TRL) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_TRL)
  )


# Recovery rate
bxp2 <- ggboxplot(Soil_, x = "Annotation", y = "RR", add = "point", title="Recovery rate Soil")
bxp2
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(RR)

Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(RR, type = "mean_sd")


Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(RR)


ggqqplot(Soil_, "RR", facet.by = "Annotation", title="Recovery rate Soil")
res.aov <- anova_test(data = Soil_, dv = RR, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_RR <- Soil_ %>%
  pairwise_t_test(
    RR ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_RR

pwc_RR <- pwc_RR %>% add_xy_position(x = "Annotation")
bxp2 + 
  stat_pvalue_manual(pwc_RR) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_RR)
  )

# RLD
bxp3 <- ggboxplot(Soil_, x = "Annotation", y = "RLD", add = "point", title="RLD Soil")
bxp3
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(RLD)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(RLD, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(RLD)

ggqqplot(Soil_, "RLD", facet.by = "Annotation", title="RLD Soil")
res.aov <- anova_test(data = Soil_, dv = RLD, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_RLD <- Soil_ %>%
  pairwise_t_test(
    RLD ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_RLD

pwc_RLD <- pwc_RLD %>% add_xy_position(x = "Annotation")
bxp3 + 
  stat_pvalue_manual(pwc_RLD) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_RLD)
  )


# HMD
bxp4 <- ggboxplot(Soil_, x = "Annotation", y = "HMD", add = "point", title="HMD Soil")
bxp4
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(HMD)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(HMD, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(HMD)

ggqqplot(Soil_, "HMD", facet.by = "Annotation", title="HMD Soil")
res.aov <- anova_test(data = Soil_, dv = HMD, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_HMD <- Soil_ %>%
  pairwise_t_test(
    HMD ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_HMD

pwc_HMD <- pwc_HMD %>% add_xy_position(x = "Annotation")
bxp4 + 
  stat_pvalue_manual(pwc_HMD) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_HMD)
  )

# Tips
bxp5 <- ggboxplot(Soil_, x = "Annotation", y = "Tips", add = "point", title="Tips Soil")
bxp5
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(Tips)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Tips, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(Tips)

ggqqplot(Soil_, "Tips", facet.by = "Annotation", title="Tips Soil")
res.aov <- anova_test(data = Soil_, dv = Tips, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Tips <- Soil_ %>%
  pairwise_t_test(
    Tips ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_Tips

pwc_Tips <- pwc_Tips %>% add_xy_position(x = "Annotation")
bxp5 + 
  stat_pvalue_manual(pwc_Tips) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_Tips)
  )


# mean_radius
bxp6 <- ggboxplot(Soil_, x = "Annotation", y = "mean_radius", add = "point", title="mean_radius Soil")
bxp6
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(mean_radius)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(mean_radius, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(mean_radius)

ggqqplot(Soil_, "mean_radius", facet.by = "Annotation", title="mean_radius Soil")
res.aov <- anova_test(data = Soil_, dv = mean_radius, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_mean_radius <- Soil_ %>%
  pairwise_t_test(
    mean_radius ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_mean_radius

pwc_mean_radius <- pwc_mean_radius %>% add_xy_position(x = "Annotation")
bxp6 + 
  stat_pvalue_manual(pwc_mean_radius) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_mean_radius)
  )
# Krs_c
bxp7 <- ggboxplot(Soil_, x = "Annotation", y = "Krs_c", add = "point", title="Krs_c Soil")
bxp7
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(Krs_c)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Krs_c, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(Krs_c)

ggqqplot(Soil_, "Krs_c", facet.by = "Annotation", title="Krs_c Soil")
res.aov <- anova_test(data = Soil_, dv = Krs_c, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Krs_c <- Soil_ %>%
  pairwise_t_test(
    Krs_c ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_Krs_c

pwc_Krs_c <- pwc_Krs_c %>% add_xy_position(x = "Annotation")
bxp7 + 
  stat_pvalue_manual(pwc_Krs_c) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_Krs_c)
  )


# Krs_v
bxp8 <- ggboxplot(Soil_, x = "Annotation", y = "Krs_v", add = "point", title="Krs_v Soil")
bxp8
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(Krs_v)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Krs_v, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(Krs_v)

ggqqplot(Soil_, "Krs_v", facet.by = "Annotation", title="Krs_v Soil")
res.aov <- anova_test(data = Soil_, dv = Krs_v, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Krs_v <- Soil_ %>%
  pairwise_t_test(
    Krs_v ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_Krs_v

pwc_Krs_v <- pwc_Krs_v %>% add_xy_position(x = "Annotation")
bxp8 + 
  stat_pvalue_manual(pwc_Krs_v) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_Krs_v)
  )



# Network_conductance_c
bxp9 <- ggboxplot(Soil_, x = "Annotation", y = "Network_conductance_c", add = "point", title="Network_conductance_c Soil")
bxp9
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(Network_conductance_c)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Network_conductance_c, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(Network_conductance_c)

ggqqplot(Soil_, "Network_conductance_c", facet.by = "Annotation", title="Network_conductance_c Soil")
res.aov <- anova_test(data = Soil_, dv = Network_conductance_c, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Network_conductance_c <- Soil_ %>%
  pairwise_t_test(
    Network_conductance_c ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_Network_conductance_c

pwc_Network_conductance_c <- pwc_Network_conductance_c %>% add_xy_position(x = "Annotation")
bxp9 + 
  stat_pvalue_manual(pwc_Network_conductance_c) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_Network_conductance_c)
  )


# Network_conductance_v
bxp10 <- ggboxplot(Soil_, x = "Annotation", y = "Network_conductance_v", add = "point", title="Network_conductance_v Soil")
bxp10
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(Network_conductance_v)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Network_conductance_v, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(Network_conductance_v)

ggqqplot(Soil_, "Network_conductance_v", facet.by = "Annotation", title="Network_conductance_v Soil")
res.aov <- anova_test(data = Soil_, dv = Network_conductance_v, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Network_conductance_v <- Soil_ %>%
  pairwise_t_test(
    Network_conductance_v ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_Network_conductance_v

pwc_Network_conductance_v <- pwc_Network_conductance_v %>% add_xy_position(x = "Annotation")
bxp10 + 
  stat_pvalue_manual(pwc_Network_conductance_v) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_Network_conductance_v)
  )

# zSUFc
bxp11 <- ggboxplot(Soil_, x = "Annotation", y = "zSUFc", add = "point", title="zSUFc Soil")
bxp11
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(zSUFc)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(zSUFc, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(zSUFc)

ggqqplot(Soil_, "zSUFc", facet.by = "Annotation", title="zSUFc Soil")
res.aov <- anova_test(data = Soil_, dv = zSUFc, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_zSUFc <- Soil_ %>%
  pairwise_t_test(
    zSUFc ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_zSUFc

pwc_zSUFc <- pwc_zSUFc %>% add_xy_position(x = "Annotation")
bxp11 + 
  stat_pvalue_manual(pwc_zSUFc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_zSUFc)
  )


# zSUFv
bxp12 <- ggboxplot(Soil_, x = "Annotation", y = "zSUFv", add = "point", title="zSUFv Soil")
bxp12
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(zSUFv)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(zSUFv, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(zSUFv)

ggqqplot(Soil_, "zSUFv", facet.by = "Annotation", title="zSUFv Soil")
res.aov <- anova_test(data = Soil_, dv = zSUFv, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_zSUFv <- Soil_ %>%
  pairwise_t_test(
    zSUFv ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_zSUFv

pwc_zSUFv <- pwc_zSUFv %>% add_xy_position(x = "Annotation")
bxp12 + 
  stat_pvalue_manual(pwc_zSUFv) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_zSUFv)
  )

# Annotation rate
bxp13 <- ggboxplot(Soil_time, x = "Annotation_time", y = "Annotation_rate", add = "point", title="Annotation rate Soil")
bxp13

Soil_time %>%
  group_by(Annotation_time) %>%
  identify_outliers(Annotation_rate)

Soil_time %>%
  group_by(Annotation_time) %>%
  get_summary_stats(Annotation_rate, type = "mean_sd")

Soil_time %>%
  group_by(Annotation_time) %>%
  shapiro_test(Annotation_rate)

ggqqplot(Soil_time, "Annotation_rate", facet.by = "Annotation_time")

res.aov <- anova_test(data = Soil_time, dv = Annotation_rate, wid = Plant_time, within = Annotation_time)
get_anova_table(res.aov)

pwc_Annotation_rate <- Soil_time %>%
  pairwise_t_test(
    Annotation_rate ~ Annotation_time, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_Annotation_rate

pwc_Annotation_rate <- pwc_Annotation_rate %>% add_xy_position(x = "Annotation")
bxp13 + 
  stat_pvalue_manual(pwc_Annotation_rate) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_Annotation_rate)
  )

Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Tips, type = "common")


res.fried <- Soil_ %>% friedman_test(Tips ~ Annotation |Plant)
res.fried

pwc_WC_Tips <- Soil_ %>%
  wilcox_test(Tips ~ Annotation, paired = TRUE, p.adjust.method = "bonferroni")
pwc_WC_Tips

pwc_WC_Tips <- pwc_WC_Tips %>% add_xy_position(x = "Annotation")
ggboxplot(Soil_, x = "Annotation", y = "Tips", add = "point") +
  stat_pvalue_manual(pwc_WC_Tips, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.fried,  detailed = TRUE),
    caption = get_pwc_label(pwc_WC_Tips)
  )

# LTipsone
bxp14 <- ggboxplot(Soil_, x = "Annotation", y = "LTipsone", add = "point", title="LTipsone Soil_")
bxp14
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(LTipsone)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipsone, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipsone)

ggqqplot(Soil_, "LTipsone", facet.by = "Annotation", title="LTipsone Soil_")
res.aov <- anova_test(data = Soil_, dv = LTipsone, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_LTipsone <- Soil_ %>%
  pairwise_t_test(
    LTipsone ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_LTipsone

pwc_LTipsone <- pwc_LTipsone %>% add_xy_position(x = "Annotation")
bxp14 + 
  stat_pvalue_manual(pwc_LTipsone) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_LTipsone)
  )


# LTipstwo
bxp15 <- ggboxplot(Soil_, x = "Annotation", y = "LTipstwo", add = "point", title="LTipstwo Soil_")
bxp15
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(LTipstwo)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipstwo, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipstwo)

ggqqplot(Soil_, "LTipstwo", facet.by = "Annotation", title="LTipstwo Soil_")
res.aov <- anova_test(data = Soil_, dv = LTipstwo, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_LTipstwo <- Soil_ %>%
  pairwise_t_test(
    LTipstwo ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_LTipstwo

pwc_LTipstwo <- pwc_LTipstwo %>% add_xy_position(x = "Annotation")
bxp15 + 
  stat_pvalue_manual(pwc_LTipstwo) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_LTipstwo)
  )

# LTipsthree
bxp16 <- ggboxplot(Soil_, x = "Annotation", y = "LTipsthree", add = "point", title="LTipsthree Soil_")
bxp16
Soil_ %>%
  group_by(Annotation) %>%
  identify_outliers(LTipsthree)
Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipsthree, type = "mean_sd")
Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipsthree)

ggqqplot(Soil_, "LTipsthree", facet.by = "Annotation", title="LTipsthree Soil_")
res.aov <- anova_test(data = Soil_, dv = LTipsthree, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_LTipsthree <- Soil_ %>%
  pairwise_t_test(
    LTipsthree ~ Annotation, paired = TRUE,
    p.adjust.method = "holm",
    alternative ="two.sided"
  )
pwc_LTipsthree

pwc_LTipsthree <- pwc_LTipsthree %>% add_xy_position(x = "Annotation")
bxp16 + 
  stat_pvalue_manual(pwc_LTipsthree) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_LTipsthree)
  )



Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipsone, type = "common")


Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipsone)

res.fried <- Soil_ %>% friedman_test(LTipsone ~ Annotation |Plant)
res.fried

pwc_WC_LTipsone <- Soil_ %>%
  wilcox_test(LTipsone ~ Annotation, paired = TRUE, p.adjust.method = "bonferroni")
pwc_WC_LTipsone

pwc_WC_LTipsone <- pwc_WC_LTipsone %>% add_xy_position(x = "Annotation")
ggboxplot(Soil_, x = "Annotation", y = "LTipsone", add = "point") +
  stat_pvalue_manual(pwc_WC_LTipsone, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.fried,  detailed = TRUE),
    caption = get_pwc_label(pwc_WC_LTipsone)
  )

Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipstwo, type = "common")

Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipstwo) 

res.fried <- Soil_ %>% friedman_test(LTipstwo ~ Annotation |Plant)
res.fried

pwc_WC_LTipstwo <- Soil_ %>%
  wilcox_test(LTipstwo ~ Annotation, paired = TRUE, p.adjust.method = "bonferroni")
pwc_WC_LTipstwo

pwc_WC_LTipstwo <- pwc_WC_LTipstwo %>% add_xy_position(x = "Annotation")
ggboxplot(Soil_, x = "Annotation", y = "LTipstwo", add = "point") +
  stat_pvalue_manual(pwc_WC_LTipstwo, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.fried,  detailed = TRUE),
    caption = get_pwc_label(pwc_WC_LTipstwo)
  )


Soil_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipsthree, type = "common")

Soil_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipsthree)

res.fried <- Soil_ %>% friedman_test(LTipsthree ~ Annotation |Plant)
res.fried

pwc_WC_LTipsthree <- Soil_ %>%
  wilcox_test(LTipsthree ~ Annotation, paired = TRUE, p.adjust.method = "bonferroni")
pwc_WC_LTipsthree

pwc_WC_LTipsthree <- pwc_WC_LTipsthree %>% add_xy_position(x = "Annotation")
ggboxplot(Soil_, x = "Annotation", y = "LTipsthree", add = "point") +
  stat_pvalue_manual(pwc_WC_LTipsthree, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.fried,  detailed = TRUE),
    caption = get_pwc_label(pwc_WC_LTipsthree)
  )
