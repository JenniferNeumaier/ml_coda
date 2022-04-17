

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

# Declare abundance for filtering
percent.filter <- args[2]

# Declare mean relative abundance
relabund.filter <- args[3]

# Declare prefix for output files
output_path <- args[4]


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
#              DATA FILTERING                    # 
#                                                #
#------------------------------------------------#


filtered_data <- taxa.filter(count_data, percent.filter = percent.filter, relabund.filter = relabund.filter)

print("Successfully filtered data!")

# Output ----
#------------------------------------------------#
#                                                #
#              OUTPUT PHASE                      # 
#                                                #
#------------------------------------------------#

write.table(filtered_data, output_path, row.names = FALSE)

print("Finished script!")