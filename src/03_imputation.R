

# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)            
cat(args, sep = "\n")


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
print("Importing data...")

count_data <- read.table(data_file, header = TRUE)

print("Imported data!")

# Data transformations ----
#------------------------------------------------#
#                                                #
#              DATA IMPUTATION                   # 
#                                                #
#------------------------------------------------#


imputed_data <- cbind(count_data[1], cmultRepl(count_data[,2:ncol(count_data)], output = method))

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