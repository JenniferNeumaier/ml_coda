# TODO:
# - NB! rbind() in fusion of test and train set will probably not work if two different denominators for test and train 
#       set has been chosen -> rewrite my_ALR function to also re-name columns and removes both denominators from final data set
# - put renaming of column names potentially in another script

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
file_path <- "../data/processed/filtered/EstMB_abundances_95.txt"

# Declare file path for metadata
metadata_file <- "../data/EstMB/EstMB_metadata.txt"

# Declare Predictor
predictor <- "I11" 

# Transformation
pre_processing <- "ALRr"

# Number of repeats per data set
n_repeats <- 5

# Declare model type
model_name <- "glmnet"  #glmnet or xgboost

# Declare type of loss function
loss_function <- "AUC" # AUC or RMSE

# Declare number of kfolds
kfold <- 5

# Declare number if CV times
cv_times <- 5

# Declare amount of training fraction
training_frac <- as.numeric(0.8)

# Declare seed 
seed <- 2022

# Method for Imputation
method <- "p-counts"


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
library("zCompositions")

# Load src documents
source("./convenience.R")
source("./data_analysis.R")


# Import data
data_set <- read.table(file_path, header = TRUE)
metadata <- read.table(metadata_file, header = TRUE)


# Data preparation ----
#------------------------------------------------#
#                                                #
#              DATA PREPARATION                  # 
#                                                #
#------------------------------------------------#

# rename column names
# colnames(data_set) <- gsub("\\.", "_", colnames(data_set))
metadata <- metadata %>%
  mutate(E11 = ifelse(E11 == "1", "DT2", "healthy")) %>%
  mutate(F41 = ifelse(F41 == "1", "PD", "healthy")) %>%
  mutate(I11 = ifelse(I11 == "1", "HHD", "healthy"))

# set variables
merge_id <- names(metadata[1])
# leftover_factors <- outersect(c(merge_id, predictor, "Country"), names(metadata))
leftover_factors <- outersect(c(merge_id, predictor), names(metadata))

# First merge
data_set <- data_set %>% 
  left_join(metadata, by = merge_id) %>% 
  relocate(c(names(metadata))) %>% 
  dplyr::select(-c(merge_id, leftover_factors))


# Split the data
set.seed(2022)
data_split <- initial_split(data_set, prop = training_frac, strata = predictor) 
train_data <- training(data_split) %>% as.data.frame()
test_data <- testing(data_split) %>% as.data.frame()

# save integers that control training set
idx <- data_split$in_id

# Data imputation ----
#------------------------------------------------#
#                                                #
#              DATA IMPUTATION                   # 
#                                                #
#------------------------------------------------#


train_data <- cbind(train_data[1], cmultRepl(train_data[,2:ncol(train_data)], output = method))

train_data <- train_data %>%
  mutate_if(is.numeric, round, digits=3) %>%
  dplyr::select(-X.1)


test_data <- cbind(test_data[1], cmultRepl(test_data[,2:ncol(test_data)], output = method))

test_data <- test_data %>%
  mutate_if(is.numeric, round, digits=3) %>%
  dplyr::select(-X.1)


print("Imputed data!")

# Data transformation ----
#------------------------------------------------#
#                                                #
#              DATA TRANSFORMATION               # 
#                                                #
#------------------------------------------------#

if(pre_processing == "ALR_optimal"){
  
  set.seed(seed)

  # split train data into IDs and features
  id <- train_data[1]
  features <- train_data[2:ncol(train_data)]
  
  # # find ALR reference
  alr.refs <- FINDALR(features)

  names(alr.refs$procrust.cor) <- colnames(features)
  res_procrustes <- alr.refs$procrust.cor[order(alr.refs$procrust.cor, decreasing = TRUE)][1:20]

  names(alr.refs$var.log) <- colnames(features)
  res_variance <- alr.refs$var.log[order(alr.refs$var.log)][1:20]

  diff <- intersect(names(res_procrustes), names(res_variance))
  denominator <- diff[1]
  
  ALR_list <- ALR(features, denom = which(colnames(features) == denominator))
  train_data <- cbind(id, as.data.frame(ALR_list$LR))
  colnames(train_data) <- gsub("\\/.*", "", colnames(train_data))
  
  # split test data into IDs and features 
  id <- test_data[1]
  features <- test_data[2:ncol(test_data)]
  
  # use ALR reference from train data
  ALR_list <- ALR(features, denom = which(colnames(features) == denominator))
  test_data <- cbind(id, as.data.frame(ALR_list$LR))
  colnames(test_data) <- gsub("\\/.*", "", colnames(test_data))

  } else if(pre_processing == "ALR_worst") {
    set.seed(seed)
    
    # split train data into IDs and features
    id <- train_data[1]
    features <- train_data[2:ncol(train_data)]
    
    # # find ALR reference
    alr.refs <- FINDALR(features)
    
    names(alr.refs$procrust.cor) <- colnames(features)
    res_procrustes <- alr.refs$procrust.cor[order(alr.refs$procrust.cor)][1:20]
    
    names(alr.refs$var.log) <- colnames(features)
    res_variance <- alr.refs$var.log[order(alr.refs$var.log)][1:20]
    
    diff <- intersect(names(res_procrustes), names(res_variance, decreasing = TRUE))
    denominator <- diff[1]
    
    ALR_list <- ALR(features, denom = which(colnames(features) == denominator))
    train_data <- cbind(id, as.data.frame(ALR_list$LR))
    colnames(train_data) <- gsub("\\/.*", "", colnames(train_data))
    
    # split test data into IDs and features 
    id <- test_data[1]
    features <- test_data[2:ncol(test_data)]
    
    # use ALR reference from train data
    ALR_list <- ALR(features, denom = which(colnames(features) == denominator))
    test_data <- cbind(id, as.data.frame(ALR_list$LR))
    colnames(test_data) <- gsub("\\/.*", "", colnames(test_data))
    
  } else if(pre_processing == "ALR_random") {
    set.seed(seed)
    
    # split train data into IDs and features
    id <- train_data[1]
    features <- train_data[2:ncol(train_data)]
    
    denominator <- sample(1:length(features), 1)
    
    ALR_list <- ALR(features, denom = denominator)
    train_data <- cbind(id, as.data.frame(ALR_list$LR))
    colnames(train_data) <- gsub("\\/.*", "", colnames(train_data))
    
    # split test data into IDs and features 
    id <- test_data[1]
    features <- test_data[2:ncol(test_data)]
    
    # use ALR reference from train data
    ALR_list <- ALR(features, denom = denominator)
    test_data <- cbind(id, as.data.frame(ALR_list$LR))
    colnames(test_data) <- gsub("\\/.*", "", colnames(test_data))
    
  } else {
  
  set.seed(seed)
  train_data <- transformation(train_data, pre_processing)
  test_data <- transformation(test_data, pre_processing)
  
}

# combine train and test set 
data_set <- rbind(train_data, test_data)



# Mikropml ----
#------------------------------------------------#
#                                                #
#                MIKROPML                       # 
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

write.table(final, paste("../out/simputation_stransformation", predictor, model_name, loss_function, pre_processing, ".txt", sep = "_", collapse = NULL), row.names = FALSE)

