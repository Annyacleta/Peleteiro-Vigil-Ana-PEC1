# Preparación -------------------------------------------------------------
if (!require(FactoMineR)) install.packages("FactoMineR", dep=TRUE)
if (!require(factoextra)) install.packages("factoextra", dep=TRUE)
if (!require(tibble)) install.packages("tibble", dep=TRUE)
if (!require(tidyr)) install.packages("tidyr", dep=TRUE)
if (!require(patchwork)) install.packages("patchwork", dep=TRUE)

if (!require(ade4)) install.packages("ade4", dep=TRUE)
if (!require(BiocManager))
  install.packages("BiocManager", dep=TRUE)
if (!require(made4))
  BiocManager::install("made4")


# Cargo librerías ---------------------------------------------------------

library(SummarizedExperiment)
library(tidyverse)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(patchwork)

# Inspección del dataset --------------------------------------------------

# Cargo objeto SE con dataset
load("SummarizedExperiment_metabolomics.Rda")

# Inspecciono el objeto
se

# Compruebo las dimensiones
dim(se)

# Exploro 3 primeras filas (metabolitos) para todas las muestras
assay(se)[1:3, ]

# ROWDATA -metadatos de metabolitos (6 primeros)
head(rowData(se))

# COLDATA -metadatos de muestras (6 primeros)
head(colData(se))

# Matriz con datos de intensidades de picos
intensidades <- assay(se)

# Número de medidas totales
prod(dim(se))

# Número de valores faltantes
sum(is.na(intensidades))


# Análisis exploratorio de los datos --------------------------------------

# 1. Análisis univariante

## 1.1. BARPLOT para ver perfil metabolómico
# Medias de intensidad de metabolito por grupo experimental
kc_cols <- which(condition == "KC")
sc_cols <- which(condition == "SC")

kc_means <- rowMeans(intensity_matrix[, kc_cols, drop = FALSE], na.rm = TRUE)
sc_means <- rowMeans(intensity_matrix[, sc_cols, drop = FALSE], na.rm = TRUE)

#SD
kc_sd <- apply(intensity_matrix[, sc_cols, drop = FALSE], 1, sd, na.rm = TRUE)
sc_sd <- apply(intensity_matrix[, sc_cols, drop = FALSE], 1, sd, na.rm = TRUE)


# Creo dataframe largo para ggplot
df_plot <- tibble(
  metabolite = rep(rownames(intensity_matrix), times = 2),
  group = rep(c("KC", "SC"), each = nrow(intensity_matrix)),
  mean = c(kc_means, sc_means),
  sd = c(kc_sd, sc_sd)
)

#Elimino valores NA
df_plot <- na.omit(df_plot)

# Obtengo una lista única de metabolitos
unique_metabs <- unique(df_plot$metabolite)
n <- length(unique_metabs)
half <- ceiling(n / 2)

# Divido en dos subconjuntos
metabs_1 <- unique_metabs[1:half]
metabs_2 <- unique_metabs[(half + 1):n]

df_plot_1 <- df_plot %>% filter(metabolite %in% metabs_1)
df_plot_2 <- df_plot %>% filter(metabolite %in% metabs_2)


p1 <- ggplot(df_plot_1, aes(x = metabolite, y = mean, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), na.rm = TRUE) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.9),
                width = 0.4,
                na.rm = TRUE) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    x = "",
    y = "Intensidad",
    fill = "Grupo"
  ) +
  ylim(0, 5)

# Segundo gráfico
p2 <- ggplot(df_plot_2, aes(x = metabolite, y = mean, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), na.rm = TRUE) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.9),
                width = 0.4,
                na.rm = TRUE) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    x = "Metabolito",
    y = "Intensidad",
    fill = "Grupo"
  ) +
  ylim(0, 5)

p1 + p2 + plot_layout(ncol = 1)


## 1.2. BOXPLOT -> uridina, para ver dispersión -> df originaaaalll, no medias!!!
intesity_df <- data.frame(intensity_matrix)

#Elimino valores NA
intesity_df <- na.omit(intesity_df)

ggplot(intesity_df, aes(x = metabolite_name, y = intensity, fill = condition)) +
  geom_boxplot() +
  labs(x = "Metabolite",
       y = "Intensidad",
       fill = "Condition") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Análisis multivariante

