# __author__ = "Zain Arifin"
# __affiliation__ = "University College Dublin and University California Santa Barbara"
# __email__ = "muhammad.arifin@ucdconnect.ie"
# __github__ = "github.com/ZainulArifin1"
# __date__ = "August 17, 2023"

#remove duplicated column and fix column name

# Load necessary libraries
library(tidyverse)
library(edgeR)

#input <- read_table("output/run64_k27/rawCount/raw_count_mapq30_BlClean_64.txt")
input_file <- snakemake@input[[1]]

# Read input data
input <- read_table(input_file)

input <- input %>% remove_rownames %>% column_to_rownames(var="Geneid")
## edgeR:: calcNormFactors
NormFactor <- calcNormFactors(object = input, method = "TMM")

## if you prefer to use the DESeq2 strategy use method="RLE" instead

## raw library size:
LibSize <- colSums(input)

## calculate size factors:
SizeFactors <- NormFactor * LibSize / 1000000

## Reciprocal, please read section below:   
SizeFactors.Reciprocal <- 1/SizeFactors %>% as.data.frame()


write.table(SizeFactors.Reciprocal, snakemake@output[[1]], 
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)