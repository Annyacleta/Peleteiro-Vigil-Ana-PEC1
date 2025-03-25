# Preparación -------------------------------------------------------------
if (!require(FactoMineR)) install.packages("FactoMineR", dep=TRUE)
if (!require(factoextra)) install.packages("factoextra", dep=TRUE)
if (!require(ade4))
  install.packages("ade4", dep=TRUE)
if (!require(BiocManager))
  install.packages("BiocManager", dep=TRUE)
if (!require(made4))
  BiocManager::install("made4")


# Cargo librerías ---------------------------------------------------------

library(SummarizedExperiment)
library(tidyverse)


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

# Media de intensidad de metabolito por grupo


# 2. Análisis multivariante

