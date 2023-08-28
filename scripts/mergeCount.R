# __author__ = "Zain Arifin"
# __affiliation__ = "University College Dublin and University California Santa Barbara"
# __email__ = "muhammad.arifin@ucdconnect.ie"
# __github__ = "github.com/ZainulArifin1"
# __date__ = "August 17, 2023"

#remove duplicated column and fix column name

# Load necessary libraries
library(tidyverse)

# Parse command-line arguments
input_file1 <- snakemake@input[[1]]
input_file2 <- snakemake@input[[2]]

# Read input data
input1 <- read_table(input_file1)
input2 <- read_table(input_file2)


full_count <- left_join(input1, input2, by = "Geneid")

write.table(full_count , snakemake@output[[1]], 
    sep = "\t", quote = FALSE, row.names = F, col.names= T)