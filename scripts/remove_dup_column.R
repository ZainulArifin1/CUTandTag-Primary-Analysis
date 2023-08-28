# __author__ = "Zain Arifin"
# __affiliation__ = "University College Dublin and University California Santa Barbara"
# __email__ = "muhammad.arifin@ucdconnect.ie"
# __github__ = "github.com/ZainulArifin1"
# __date__ = "August 17, 2023"

#remove duplicated column and fix column name

# Load necessary libraries
library(tidyverse)

# Parse command-line arguments
input_file <- snakemake@input[[1]]

# Read input data
input <- read_table(input_file)

# Remove duplicate columns
input <- input[, !duplicated(as.list(input))]

# Rename columns
colnames(input) <- gsub(".*/", "", colnames(input))

# Write processed data to output file
write.table(input, snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE)