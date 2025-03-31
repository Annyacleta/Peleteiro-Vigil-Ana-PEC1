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

## 1.1. BARPLOT para ver perfil metabolómico (no mostrado en informe)
# Medias de intensidad de metabolito por grupo experimental
kc_cols <- which(condition == "KC")
sc_cols <- which(condition == "SC")

kc_means <- rowMeans(intensidades[, kc_cols, drop = FALSE], na.rm = TRUE)
sc_means <- rowMeans(intensidades[, sc_cols, drop = FALSE], na.rm = TRUE)

#SD
kc_sd <- apply(intensidades[, sc_cols, drop = FALSE], 1, sd, na.rm = TRUE)
sc_sd <- apply(intensidades[, sc_cols, drop = FALSE], 1, sd, na.rm = TRUE)


# Creo dataframe largo para ggplot
df_plot <- tibble(
  metabolite = rep(rownames(intensidades), times = 2),
  group = rep(c("KC", "SC"), each = nrow(intensidades)),
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
  ylim(0, 3.5)

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
    fill = "Condición"
  ) +
  ylim(0, 3.5)

plot_final <- p1 + p2 + plot_layout(ncol = 1)

# Guardo el gráfico
ggsave("barplot_metabolites.png", plot = plot_final, width = 8, height = 12)

## 1.2. BOXPLOT -> con df original
intesity_df <- data.frame(intensidades)

# Obtengo metadatos de columna (condición experimental)
condition <- colData(se)$condition

# Transformo df a formato largo
intesity_df <- as.data.frame(intensidades) %>%
  rownames_to_column(var = "metabolite") %>%
  pivot_longer(-metabolite, names_to = "sample", values_to = "intensity")

# Añado condición experimental a cada muestra
intesity_df$condition <- condition[match(intesity_df$sample, colnames(intensidades))]

#Elimino valores NA
intesity_df <- na.omit(intesity_df)

#Divido en 2 el df
df_boxplot_1 <- intesity_df %>% filter(metabolite %in% metabs_1)
df_boxplot_2 <- intesity_df %>% filter(metabolite %in% metabs_2)

# Boxplot 1
box1 <- ggplot(df_boxplot_1, aes(x = metabolite, y = intensity, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
    labs(x = "",
       y = "Abundancia",
       fill = "Condición") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 2)

# Boxplot 2
box2 <- ggplot(df_boxplot_2, aes(x = metabolite, y = intensity, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "",
       y = "Abundandica",
       fill = "Condición") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 4)

boxplot_final <- box1 + box2 + plot_layout(ncol = 1, heights = c(1,2))

# Guardo el gráfico
ggsave("boxplot_metabolites.png", plot = boxplot_final, width = 8, height = 12)

#View(intesity_df) #Quitar L-lactic acid, se sale demasiado!!!

# 2. Análisis multivariante

