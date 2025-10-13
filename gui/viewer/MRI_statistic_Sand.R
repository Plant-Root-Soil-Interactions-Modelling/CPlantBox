#SAND
library(tidyverse)
library(ggpubr)
library(rstatix)

Plant <- c("1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "4")
Annotation <- c("M", "M", "M", "M", "Mp", "Mp", "Mp", "Mp", "A",  "A",  "A",  "A")
TRL <- c(64.57, 148.28, 98.23, 55.93, 84.09, 185.73, 97.46, 74.13, 56.29, 138.63, 84.678, 58.82)
RR <-c(47.13,	63.64,	93.55,	57.66,	61.38,	79.71,	92.82,	76.42,	41.09,	59.50,	80.65,	60.64)
RLD <-c(0.141,	0.325,	0.215,	0.122,	0.18,	0.41,	0.21,	0.16,	0.12,	0.30,	0.19,	0.13)
HMD <-c(1.50,	0.99,	1.21,	1.61,	1.31,	0.88,	1.22,	1.40,	1.60,	1.02,	1.31,	1.57)
Tips <-c(17,	37,	18,	15,	38,	55,	19,	32,	33,	111,	36,	31)
LTipsone <-c(8, 32, 14, 13, 17, 50, 14, 29, 13, 35, 12, 21)
LTipstwo <-c(8, 4, 3, 1, 15, 4, 4, 2, 16, 40, 9, 6)
LTipsthree <-c(0, 0, 0, 0, 5, 0, 0, 0, 3, 26, 12, 3)
mean_radius <-c(0.0248,	0.0232,	0.0251,	0.0305,	0.0329,	0.0298,	0.0336,	0.0391,	0.0359,	0.028,	0.0318,	0.0337)
Krs_c <-c(1.52E-03,	2.51E-03,	2.49E-03,	1.62E-03,	2.38E-03,	3.20E-03,	3.11E-03,	2.50E-03,	1.87E-03,	3.19E-03,	2.60E-03,	1.83E-03)
Krs_v <-c(7.28E-03,	1.20E-02,	1.11E-02,	9.50E-03,	1.07E-02,	1.39E-02,	1.21E-02,	1.36E-02,	9.04E-03,	1.71E-02,	5.61E-03,	1.01E-02)
Network_conductance_c <-c(2.35E-05,	1.69E-05,	2.53E-05,	2.90E-05,	2.83E-05,	1.72E-05,	3.19E-05,	3.37E-05,	3.32E-05,	2.30E-05,	3.07E-05,	3.11E-05)
Network_conductance_v <-c(1.13E-04,	8.09E-05,	1.13E-04,	1.70E-04,	1.27E-04,	7.48E-05,	1.24E-04,	1.83E-04,	1.61E-04,	1.23E-04,	6.63E-05,	1.72E-04)
zSUFc<-c( -4.05,	-5.74,	-3.51,	-4.88,	-4.16,	-5.40,	-3.39,	-4.08,	-4.17,	-5.80,	-3.76,	-4.74)
zSUFv<-c(	-3.48,	-4.68,	-1.95,	-3.91,	-3.48,	-3.78,	-1.83,	-3.13,	-3.31,	-4.25,	-2.71,	-3.90)



Sand <- data.frame(Plant, Annotation, TRL, RR, RLD, HMD, Tips, LTipsone, LTipstwo, LTipsthree, mean_radius, Krs_c, Krs_v, Network_conductance_c, Network_conductance_v, zSUFc, zSUFv)
Sand_ <-as_tibble(Sand)
Sand_ <- Sand_ %>%
  convert_as_factor(Plant, Annotation)
print (Sand_)

# only for M and M+
Plant_time <- c("1", "2", "3", "4", "1", "2", "3", "4")
Annotation_time <- c("M", "M", "M", "M", "Mp", "Mp", "Mp", "Mp")
Annotation_rate <- c(2.207,	2.5217,	4.519,	3.9111,	4.856,	5.962,	8.032,	7.127)
Sand_time <- data.frame(Plant_time, Annotation_time, Annotation_rate)
Sand_time <-as_tibble(Sand_time)
Sand_time <- Sand_time %>%
  convert_as_factor(Plant_time, Annotation_time)
print (Sand_time)

# TRL
bxp1 <- ggboxplot(Sand_, x = "Annotation", y = "TRL", add = "point", title="TRL Sand")
bxp1

Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(TRL)

Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(TRL, type = "mean_sd")

Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(TRL)



ggqqplot(Sand_, "TRL", facet.by = "Annotation", title="TRL Sand") #graue Balken sind 95% Konfidenzintervall

res.aov <- anova_test(data = Sand_, dv = TRL, wid = Plant, within = Annotation, detailed = TRUE)
get_anova_table(res.aov)

pwc_TRL <- Sand_ %>%
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
bxp2 <- ggboxplot(Sand_, x = "Annotation", y = "RR", add = "point", title="Recovery rate sand")
bxp2
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(RR)

Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(RR, type = "mean_sd")


Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(RR)


ggqqplot(Sand_, "RR", facet.by = "Annotation", title="Recovery rate sand")
res.aov <- anova_test(data = Sand_, dv = RR, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_RR <- Sand_ %>%
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
bxp3 <- ggboxplot(Sand_, x = "Annotation", y = "RLD", add = "point", title="RLD sand")
bxp3
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(RLD)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(RLD, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(RLD)

ggqqplot(Sand_, "RLD", facet.by = "Annotation", title="RLD sand")
res.aov <- anova_test(data = Sand_, dv = RLD, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_RLD <- Sand_ %>%
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
bxp4 <- ggboxplot(Sand_, x = "Annotation", y = "HMD", add = "point", title="HMD sand")
bxp4
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(HMD)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(HMD, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(HMD)

ggqqplot(Sand_, "HMD", facet.by = "Annotation", title="HMD sand")
res.aov <- anova_test(data = Sand_, dv = HMD, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_HMD <- Sand_ %>%
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
bxp5 <- ggboxplot(Sand_, x = "Annotation", y = "Tips", add = "point", title="Tips sand")
bxp5
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(Tips)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Tips, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(Tips)

ggqqplot(Sand_, "Tips", facet.by = "Annotation", title="Tips sand")
res.aov <- anova_test(data = Sand_, dv = Tips, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Tips <- Sand_ %>%
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
bxp6 <- ggboxplot(Sand_, x = "Annotation", y = "mean_radius", add = "point", title="mean_radius sand")
bxp6
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(mean_radius)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(mean_radius, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(mean_radius)

ggqqplot(Sand_, "mean_radius", facet.by = "Annotation", title="mean_radius sand")
res.aov <- anova_test(data = Sand_, dv = mean_radius, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_mean_radius <- Sand_ %>%
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
bxp7 <- ggboxplot(Sand_, x = "Annotation", y = "Krs_c", add = "point", title="Krs_c sand")
bxp7
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(Krs_c)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Krs_c, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(Krs_c)

ggqqplot(Sand_, "Krs_c", facet.by = "Annotation", title="Krs_c sand")
res.aov <- anova_test(data = Sand_, dv = Krs_c, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Krs_c <- Sand_ %>%
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
bxp8 <- ggboxplot(Sand_, x = "Annotation", y = "Krs_v", add = "point", title="Krs_v sand")
bxp8
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(Krs_v)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Krs_v, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(Krs_v)

ggqqplot(Sand_, "Krs_v", facet.by = "Annotation", title="Krs_v sand")
res.aov <- anova_test(data = Sand_, dv = Krs_v, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Krs_v <- Sand_ %>%
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
bxp9 <- ggboxplot(Sand_, x = "Annotation", y = "Network_conductance_c", add = "point", title="Network_conductance_c sand")
bxp9
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(Network_conductance_c)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Network_conductance_c, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(Network_conductance_c)

ggqqplot(Sand_, "Network_conductance_c", facet.by = "Annotation", title="Network_conductance_c sand")
res.aov <- anova_test(data = Sand_, dv = Network_conductance_c, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Network_conductance_c <- Sand_ %>%
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
bxp10 <- ggboxplot(Sand_, x = "Annotation", y = "Network_conductance_v", add = "point", title="Network_conductance_v sand")
bxp10
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(Network_conductance_v)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Network_conductance_v, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(Network_conductance_v)

ggqqplot(Sand_, "Network_conductance_v", facet.by = "Annotation", title="Network_conductance_v sand")
res.aov <- anova_test(data = Sand_, dv = Network_conductance_v, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_Network_conductance_v <- Sand_ %>%
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
bxp11 <- ggboxplot(Sand_, x = "Annotation", y = "zSUFc", add = "point", title="zSUFc sand")
bxp11
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(zSUFc)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(zSUFc, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(zSUFc)

ggqqplot(Sand_, "zSUFc", facet.by = "Annotation", title="zSUFc sand")
res.aov <- anova_test(data = Sand_, dv = zSUFc, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_zSUFc <- Sand_ %>%
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
bxp12 <- ggboxplot(Sand_, x = "Annotation", y = "zSUFv", add = "point", title="zSUFv sand")
bxp12
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(zSUFv)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(zSUFv, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(zSUFv)

ggqqplot(Sand_, "zSUFv", facet.by = "Annotation", title="zSUFv sand")
res.aov <- anova_test(data = Sand_, dv = zSUFv, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_zSUFv <- Sand_ %>%
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
bxp13 <- ggboxplot(Sand_time, x = "Annotation_time", y = "Annotation_rate", add = "point", title="Annotation rate sand")
bxp13

Sand_time %>%
  group_by(Annotation_time) %>%
  identify_outliers(Annotation_rate)

Sand_time %>%
  group_by(Annotation_time) %>%
  get_summary_stats(Annotation_rate, type = "mean_sd")

Sand_time %>%
  group_by(Annotation_time) %>%
  shapiro_test(Annotation_rate)

ggqqplot(Sand_time, "Annotation_rate", facet.by = "Annotation_time")

res.aov <- anova_test(data = Sand_time, dv = Annotation_rate, wid = Plant_time, within = Annotation_time)
get_anova_table(res.aov)

pwc_Annotation_rate <- Sand_time %>%
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



Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(Tips, type = "common")

Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(Tips)

res.fried <- Sand_ %>% friedman_test(Tips ~ Annotation |Plant)
res.fried

pwc_WC_Tips <- Sand_ %>%
  wilcox_test(Tips ~ Annotation, paired = TRUE, p.adjust.method = "bonferroni")
pwc_WC_Tips

pwc_WC_Tips <- pwc_WC_Tips %>% add_xy_position(x = "Annotation")
ggboxplot(Sand_, x = "Annotation", y = "Tips", add = "point") +
  stat_pvalue_manual(pwc_WC_Tips, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.fried,  detailed = TRUE),
    caption = get_pwc_label(pwc_WC_Tips)
  )


# LTipsone
bxp14 <- ggboxplot(Sand_, x = "Annotation", y = "LTipsone", add = "point", title="LTipsone sand")
bxp14
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(LTipsone)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipsone, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipsone)

ggqqplot(Sand_, "LTipsone", facet.by = "Annotation", title="LTipsone sand")
res.aov <- anova_test(data = Sand_, dv = LTipsone, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_LTipsone <- Sand_ %>%
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
bxp15 <- ggboxplot(Sand_, x = "Annotation", y = "LTipstwo", add = "point", title="LTipstwo sand")
bxp15
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(LTipstwo)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipstwo, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipstwo)

ggqqplot(Sand_, "LTipstwo", facet.by = "Annotation", title="LTipstwo Sand")
res.aov <- anova_test(data = Sand_, dv = LTipstwo, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_LTipstwo <- Sand_ %>%
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
bxp16 <- ggboxplot(Sand_, x = "Annotation", y = "LTipsthree", add = "point", title="LTipsthree Sand")
bxp16
Sand_ %>%
  group_by(Annotation) %>%
  identify_outliers(LTipsthree)
Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipsthree, type = "mean_sd")
Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipsthree)

ggqqplot(Sand_, "LTipsthree", facet.by = "Annotation", title="LTipsthree sand")
res.aov <- anova_test(data = Sand_, dv = LTipsthree, wid = Plant, within = Annotation)
get_anova_table(res.aov)

pwc_LTipsthree <- Sand_ %>%
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

Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipsone, type = "common")


Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipsone)

res.fried <- Sand_ %>% friedman_test(LTipsone ~ Annotation |Plant)
res.fried

pwc_WC_LTipsone <- Sand_ %>%
  wilcox_test(LTipsone ~ Annotation, paired = TRUE, p.adjust.method = "bonferroni")
pwc_WC_LTipsone

pwc_WC_LTipsone <- pwc_WC_LTipsone %>% add_xy_position(x = "Annotation")
ggboxplot(Sand_, x = "Annotation", y = "LTipsone", add = "point") +
  stat_pvalue_manual(pwc_WC_LTipsone, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.fried,  detailed = TRUE),
    caption = get_pwc_label(pwc_WC_LTipsone)
  )

Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipstwo, type = "common")

Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipstwo) 

res.fried <- Sand_ %>% friedman_test(LTipstwo ~ Annotation |Plant)
res.fried

pwc_WC_LTipstwo <- Sand_ %>%
  wilcox_test(LTipstwo ~ Annotation, paired = TRUE, p.adjust.method = "bonferroni")
pwc_WC_LTipstwo

pwc_WC_LTipstwo <- pwc_WC_LTipstwo %>% add_xy_position(x = "Annotation")
ggboxplot(Sand_, x = "Annotation", y = "LTipstwo", add = "point") +
  stat_pvalue_manual(pwc_WC_LTipstwo, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.fried,  detailed = TRUE),
    caption = get_pwc_label(pwc_WC_LTipstwo)
  )


Sand_ %>%
  group_by(Annotation) %>%
  get_summary_stats(LTipsthree, type = "common")

Sand_ %>%
  group_by(Annotation) %>%
  shapiro_test(LTipsthree)

res.fried <- Sand_ %>% friedman_test(LTipsthree ~ Annotation |Plant)
res.fried

pwc_WC_LTipsthree <- Sand_ %>%
  wilcox_test(LTipsthree ~ Annotation, paired = TRUE, p.adjust.method = "bonferroni")
pwc_WC_LTipsthree

pwc_WC_LTipsthree <- pwc_WC_LTipsthree %>% add_xy_position(x = "Annotation")
ggboxplot(Sand_, x = "Annotation", y = "LTipsthree", add = "point") +
  stat_pvalue_manual(pwc_WC_LTipsthree, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.fried,  detailed = TRUE),
    caption = get_pwc_label(pwc_WC_LTipsthree)
  )

