library(tidyverse)
library(readxl)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(here)

#1. Data upload
#Dataset with measured root and shoot traits
P_levels <- as.data.frame(read_excel(here("P_CA.xlsx")) %>%
                            mutate(treatment = factor(treatment))) 
P_levels <- P_levels %>%
  mutate(P_level = as.factor(str_split(plant, "_", simplify = TRUE)[, 1]))
#Dataset with modeled Krs 
krs <- read_excel(here("krs.xlsx"))
#some data manipulation to obtain the format needed for PCA
krs <- krs %>%
  pivot_longer(cols = -day, names_sep = "_", names_to = c("trait", "P_level")) %>%
  filter(day %in% c(28)) %>% #filter the days you need 
  pivot_wider(names_from = c(day, trait), values_from = value)
#Join the datasets, include the plant_names as index and eliminate unused columns
P_levels <- P_levels %>%
  left_join(krs, by = "P_level")
rownames(P_levels) <- P_levels$plant
P_levels <- P_levels %>%
  select(-c(treatment, plant)) %>%
  mutate(P_level = factor(P_level))

#2. Correlation matrix
ggcorrplot(cor(P_levels %>%
                 select(-c(P_level))), lab = TRUE) #the treatment column must be excluded
# and also any other variables that are not of interest 

#3. PCA calculation
res.pca <- PCA(P_levels %>%
                 select(-c(P_level)), graph = F) 
#select same variables (or less) as for the correlation matrix

#4. Plot the eigenvalues
get_eig(res.pca)
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

#5. Contribution of the factors to the PC axes
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE)

#6. PCA cluster (FIGURE 4)
fviz_pca_biplot(res.pca, 
                col.ind = P_levels$P_level, 
                addEllipses = TRUE,ellipse.level = .75, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Phosphor levels",
                title = "Phosphor experiment",
                palette = c("blue", "orange", "green3", "red")) +
  theme(legend.position = c(.9,.8))


