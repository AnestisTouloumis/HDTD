if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocCheck")

library("BiocCheck")
BiocCheck(".")
