

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

# Declare preprocessing step used
preprocessing <- args[2]

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

# Load source scripts
source("./convenience.R")
source("./data_analysis.R")


# Read count matrix used as input (genes, species etc)
count_data <- read.table(data_file, header = TRUE)

print("Imported data!")

# Data transformations ----
#------------------------------------------------#
#                                                #
#              DATA TRANSFORMATION               # 
#                                                #
#------------------------------------------------#

if (preprocessing == "TSS"){
  
  data_transformed <- my_TSS(count_data)

} else if (preprocessing == "CLR"){
  
  data_transformed <- my_CLR(count_data)
  
} else if (preprocessing == "ALR_optimal"){
  
  data_transformed <- my_optimal_ALR(count_data)
  
} else if (preprocessing == "ALR_worst"){
  
  data_transformed <- my_worst_ALR(count_data)
  
} else if (preprocessing == "ALR_random"){
  
  data_transformed <- my_random_ALR(count_data)
  
}

print("Successfully transformed data!")

# Output ----
#------------------------------------------------#
#                                                #
#              OUTPUT PHASE                      # 
#                                                #
#------------------------------------------------#

write.table(data_transformed, output_path, row.names = FALSE)

print("Finished script!")