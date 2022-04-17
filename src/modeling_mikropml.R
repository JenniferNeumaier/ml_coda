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
file_path <- args[1]

# Declare file path for metadata
metadata_file <- args[2]

# Declare Predictor
predictor <- args[3] 

# Number of repeats per data set
counter <- args[4]

if (counter == 0) {
  stop("Counter has to be set to at least 1!", call. = FALSE)
}

# Declare model type
model_name <- args[5]  #glmnet or xgboost

# Declare type of loss function
loss_function <- args[6] # AUC or RMSE

# Declare number of kfolds
# kfolds <- args[7]

# Declare number if CV times

# cv_times <- args[8]

# Declare amount of training fraction
training_frac <- args[7]


# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("stringr")
library("dyplr")
library("tidymodels")
library("mikropml")
library("codacore")

# Load src documents
source("./convenience.R")
source("./data_analysis.R")


# Import data
data_set <- read.table(file_path, header = TRUE)
metadata <- read.table(metadata_file, header = TRUE)

print("Imported data!")


# Data perparation ----
#------------------------------------------------#
#                                                #
#              DATA PREPARATION                  # 
#                                                #
#------------------------------------------------#

# set variables
merge_id <- names(metadata[1])
leftover_factors <- outersect(c(predictor, id), names(metadata))


# merge data used for mikrop models

data_set <- data_set %>% 
  left_join(metadata, by = merge_id) %>% 
  relocate(c(predictor, leftover_factors), .after = merge_id) %>% 
  dplyr::select(-c(merge_id, leftover_factors))



# Mikropml ----
#------------------------------------------------#
#                                                #
#                MIKROPML                         # 
#                                                #
#------------------------------------------------#

# set variables
mikrop_training <- 
mikrop_test <- list()


# go through each data set in data_list
for (i in 1:length(data_transformed)) { 
  
  # initiate temp lists
  training_res <- list()
  test_res <- list()
  
  # repeat experiment n times
  for (j in counter) { 
    
    ml_result <- run_ml(data_transformed[[i]],
                        model_name,
                        outcome_colname = predictor, 
                        kfold = 5,
                        cv_times = 10,
                        training_frac = training_frac,
                        seed = NA)
    
    # save results in temp lists
    if(loss_function == "AUC"){
      
      training_res[j] <- ml_result$performance[1]
      test_res[j] <- ml_result$performance[3]
      
    } else if (loss_function == "RMSE"){
      
      training_res[j] <- ml_result$performance[1]
      test_res[j] <- ml_result$performance[2]
      
    }
    
  }
  
  # save lists in final list
  mikrop_training[[i]] <- training_res
  mikrop_test[[i]] <- test_res
}


# output ----
#------------------------------------------------#
#                                                #
#               OUTPUT                           # 
#                                                #
#------------------------------------------------#

## mikropml
set_names <- names(data_transformed)

# curating training set
names(mikrop_training) <- set_names
training_df <- stack(mikrop_training)

training_df <- training_df %>% 
  rename(loss = values, dataset = ind) %>% 
  add_column(type = "training") %>% 
  separate(dataset, c("dataset", "filtering", "filler", "transformation"), "_") %>% 
  dplyr::select(-c(filler))

# curating test set
names(mikrop_test) <- set_names
test_df <- stack(mikrop_test)

test_df <- test_df %>% 
  rename(loss = values, dataset = ind) %>% 
  add_column(type = "test") %>% 
  separate(dataset, c("dataset", "filtering", "filler", "transformation"), "_") %>% 
  dplyr::select(-c(filler))

# combine and outout
final <- rbind(training_df, test_df)

write.table(final, paste("../out/table_codacore", predictor, model_name, loss_function, ".txt", sep = "_", collapse = NULL), row.names = FALSE)

