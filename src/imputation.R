

# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)             
#args = c("RData/Abundance_filtered_species_p10a001.rds", "RData/Data_master.rds", "LASSO", "CLR", "1", "5", "SET2", "Rdata/ML", "F41", "I10") 
cat(args, sep = "\n")


# Expected input:
#        1                  2             3     
# input_data_file      model_name       seed

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 


# Declare count matrix input
data_file <- args[1]

# Declare method for zCompositions ("p-counts"/...)
method <- args[2]

# Declare prefix for output files
output_path <- args[3]


# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("dplyr")
library("zCompositions")

# Load source scripts
source("./convenience.R")
source("./data_analysis.R")


# Read count matrix used as input (genes, species etc)
count_data <- read.table(data_file, header = TRUE)

print("Imported data!")

# Data transformations ----
#------------------------------------------------#
#                                                #
#              DATA IMPUTATION                   # 
#                                                #
#------------------------------------------------#


imputed_data <- cbind(imputed_data[1], cmultRepl(imputed_data[,2:ncol(imputed_data)], output = method))

imputed_data <- imputed_data %>% 
  mutate_if(is.numeric, round, digits=3)

print("Imputed data!")

# Output ----
#------------------------------------------------#
#                                                #
#              OUTPUT PHASE                      # 
#                                                #
#------------------------------------------------#

write.table(imputed_data, output_path, row.names = FALSE)

print("Finished script!")