# TODO:
# - rewrite for a single data set


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


if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 


# Declare file path of data
file_path <- "./data/PCOS/PCOS_abundances_10.txt"

# Declare file path for metadata
metadata_file <- "./data/PCOS/PCOS_metadata.txt"

# Declare Predictor
predictor <- "PCOS_Riikka" 

# Transformation
pre_processing <- "ALR_optimal"

# Number of repeats per data set
n_repeats <- 10

# Declare model type
model_name <- "glmnet"  #glmnet or xgboost

# Declare type of loss function
loss_function <- "AUC" # AUC or RMSE

# Declare number of kfolds
kfold <- 5

# Declare number if CV times
cv_times <- 10

# Declare amount of training fraction
training_frac <- 0.8

# Declare seed 
seed <- 2022

# Declare country as holdout set
# country <- "FRA" # GER, FRA, CHI, USA, AUS


# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("stringr")
library("dplyr")
library("tidymodels")
library("mikropml")

# Load src documents
source("./src/convenience.R")
source("./src/data_analysis.R")


# Import data
data_set <- read.table(file_path, header = TRUE)
metadata <- read.table(metadata_file, header = TRUE)

print("Imported data!")

# Data transformation ----
#------------------------------------------------#
#                                                #
#              DATA TRANSFORMATION              # 
#                                                #
#------------------------------------------------#

set.seed(seed)
data_set <- transformation(data_set, pre_processing)


# Data preparation ----
#------------------------------------------------#
#                                                #
#              DATA PREPARATION                  # 
#                                                #
#------------------------------------------------#

# set variables
merge_id <- names(metadata[1])
# leftover_factors <- outersect(c(merge_id, predictor, "Country"), names(metadata))
leftover_factors <- outersect(c(merge_id, predictor), names(metadata))

# merge data used for mikrop models
data_set <- data_set %>% 
  left_join(metadata, by = merge_id) %>% 
  relocate(c(names(metadata))) %>% 
  dplyr::select(-c(merge_id, leftover_factors))

# # Split the data
# train_data <- data_set %>% 
#   filter(!grepl(country, Country)) %>% 
#   dplyr::select(-c(Country))
# 
# test_data <- data_set %>% 
#   filter(grepl(country, Country)) %>% 
#   dplyr::select(-c(Country))  
# 
# # save ids for data not belonging to selected country
# idx <- data_set %>%
#   rownames_to_column() %>%
#   filter(!grepl(country, Country)) %>%
#   `[[`("rowname") %>%
#   as.numeric()
# 
# # remove Country column
# data_set <- data_set %>% 
#   dplyr::select(-Country)

# Split the data
set.seed(2022)
data_split <- initial_split(data_set, prop = 0.8, strata = predictor) 
train_data <- training(data_split) %>% as.data.frame()
test_data <- testing(data_split) %>% as.data.frame()

# save integers that control training set
idx <- data_split$in_id

# Mikropml ----
#------------------------------------------------#
#                                                #
#                MIKROPML                         # 
#                                                #
#------------------------------------------------#

# set variables
mikrop_training <- list()
mikrop_test <- list()

counter <- seq(1,n_repeats,1)


# repeat experiment n times
for (j in counter) { 
  
  ml_result <- run_ml(data_set,
                      model_name,
                      outcome_colname = predictor, 
                      kfold = kfold,
                      cv_times = cv_times,
                      training_frac = idx,
                      seed = NA)
  
  # save results in temp lists
  if(loss_function == "AUC"){
    
    mikrop_training[j] <- ml_result$performance[1]
    mikrop_test[j] <- ml_result$performance[3]
    
  } else if (loss_function == "RMSE"){
    
    mikrop_training[j] <- ml_result$performance[1]
    mikrop_test[j] <- ml_result$performance[2]
    
  }

}



# output ----
#------------------------------------------------#
#                                                #
#               OUTPUT                           # 
#                                                #
#------------------------------------------------#

## mikropml


# curating training set
training_df <- data.frame("loss" = unlist(mikrop_training),
                          "type" = "training",
                          "transformation" = pre_processing)
     
# curating test set
test_df <- data.frame("loss" = unlist(mikrop_test),
                          "type" = "test",
                          "transformation" = pre_processing)


# combine and outout
final <- rbind(training_df, test_df)

write.table(final, paste("./out/PCOS_leaky", predictor, model_name, loss_function, pre_processing, ".txt", sep = "_", collapse = NULL), row.names = FALSE)

