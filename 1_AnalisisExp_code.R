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


# Análisis exploratorio ---------------------------------------------------


# Cargar objeto SE con dataset

load("SummarizedExperiment_metabolomics.Rda")

# Inspección del objeto
se
