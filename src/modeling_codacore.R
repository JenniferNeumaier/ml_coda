# TODO:
# - rewrite for a single data set
# - test if output works!


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


# merge PCOS_data with metadata
data_set <- data_set %>% 
  left_join(metadata, by = merge_id) %>% 
  relocate(c(predictor, leftover_factors), .after = merge_id) %>% 
  dplyr::select(-c(merge_id, leftover_factors))

# Split the data
set.seed(2022)
data_split <- initial_split(data_set, prop = training_frac, strata = predictor) 
train_data <- training(data_split) %>% as.data.frame()
test_data <- testing(data_split) %>% as.data.frame()


# codacore ----
#------------------------------------------------#
#                                                #
#               CODACORE                         # 
#                                                #
#------------------------------------------------#

# set up variables
coda_training <- list()
coda_test <- list()


# set objective for codacore based on chosen loss_function
if(loss_function == "AUC"){
  objective <- "binary classification"
} else if(loss_function == "RMSE"){
  objective <- "regression"
}


# generate features (x) and variable (y) for train
x_train <- train_data[, 2:ncol(train_data)]
y_train <- train_data[,1] 

# generate features (x) and variable (y) for test
x_test <- test_data[, 2:ncol(test_data)]
y_test <- test_data[,1] 



for (j in counter) {
  
  # initiate temp lists
  training_res <- list()
  test_res <- list()
  
  model_bal_0 <- codacore(
    x_train, 
    y_train, 
    objective = objective,
    logRatioType = "balances",
    cvParams = c(numFolds = 5),
    lambda = 0,
    verbose = FALSE, 
    fast = TRUE
  )
  
  model_bal_1 <- codacore(
    x_train, 
    y_train,
    objective = objective,
    logRatioType = "balances",  
    cvParams = c(numFolds = 5),
    lambda = 1,
    verbose = FALSE, 
    fast = TRUE    
  )
  
  if(loss_function == "AUC"){
    
    tryCatch(
      expr = {
        
        # save training AUC in lists
        training_res[[j]] <- c(as.numeric(model_bal_0$ensemble[[1]]$AUC), 
                               as.numeric(model_bal_1$ensemble[[1]]$AUC))
        
        # get test AUC by predicting manually
        pred_bal_0_lr1 <- predict(model_bal_0, x_test, logits = FALSE, numLogRatios = 1)
        pred_bal_0 <- predict(model_bal_0, x_test, logits = FALSE)
        
        pred_bal_1_lr1 <- predict(model_bal_1, x_test, logits = FALSE, numLogRatios = 1)
        pred_bal_1 <- predict(model_bal_1, x_test, logits = FALSE)
        
        test_res[[j]] <- c(as.numeric(pROC::auc(y_test, pred_bal_0_lr1, quiet = T)), 
                           as.numeric(pROC::auc(y_test, pred_bal_0, quiet = T)), 
                           as.numeric(pROC::auc(y_test, pred_bal_1_lr1, quiet = T)),
                           as.numeric(pROC::auc(y_test, pred_bal_1, quiet = T)))
      }
    )
    
  } else if (loss_function == "RMSE"){ 
    
    tryCatch(
      expr = { 
        
        # save training RMSE in lists
        training_res[[j]] <- c(as.numeric(model_bal_0$ensemble[[1]]$RMSE), 
                               as.numeric(model_bal_1$ensemble[[1]]$RMSE))
        
        # get test RMSE by predicting manually
        pred_bal_0_lr1 <- predict(model_bal_0, x_test, logits = FALSE, numLogRatios = 1)
        pred_bal_0 <- predict(model_bal_0, x_test, logits = FALSE)
        
        pred_bal_1_lr1 <- predict(model_bal_1, x_test, logits = FALSE, numLogRatios = 1)
        pred_bal_1 <- predict(model_bal_1, x_test, logits = FALSE)
        
        
        # test_res[[j]] <- c(as.numeric(sqrt(mean((y_test - pred_bal_0_lr1)^2))), 
        #                    as.numeric(sqrt(mean((y_test - pred_bal_0)^2))), 
        #                    as.numeric(sqrt(mean((y_test - pred_bal_1_lr1)^2))),
        #                    as.numeric(sqrt(mean((y_test - pred_bal_1)^2))))
        
        
        test_res[[j]] <- c(as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_0_lr1, quiet = T)), 
                           as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_0, quiet = T)), 
                           as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_1_lr1, quiet = T)),
                           as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_1, quiet = T)))
      }
    )
  }
  
  # save lists in final list
  coda_training[[i]] <- training_res
  coda_test[[i]] <- test_res
}


## codacore

# set names
set_names <- names(data_imputed)

# curating training set
names(coda_training) <- set_names

training_df <- data.frame()
for(i in 1:length(coda_training)){
  res <- do.call(rbind, coda_training[[i]])
  res <- res %>% 
    as.data.frame %>% 
    rename(model_bal_0 = V1, model_bal_1 = V2) %>% 
    stack() %>% 
    add_column(dataset = names(data_imputed)[i], type = "training") %>% 
    rename(loss = values, model = ind) %>% 
    separate(dataset, c("dataset", "filler", "filtering"), "_") %>% 
    dplyr::select(-c(filler))
  
  training_df <- rbind(training_df, res)
  
}

# curating test set
names(coda_test) <- set_names

test_df <- data.frame()
for(i in 1:length(coda_test)){
  res <- do.call(rbind, coda_test[[i]])
  res <- res %>% 
    as.data.frame %>% 
    rename(pred_bal_0_lr1 = V1, pred_bal_0 = V2, pred_bal_1_lr1 = V3, pred_bal_1= V4) %>% 
    stack() %>% 
    add_column(dataset = names(data_imputed)[i], type = "test") %>% 
    rename(loss = values, model = ind) %>% 
    separate(dataset, c("dataset", "filler", "filtering"), "_") %>% 
    dplyr::select(-c(filler))
  
  test_df <- rbind(test_df, res)
  
}

# combine and output
final <- rbind(training_df, test_df)

write.table(final, paste("../out/table_codacore", predictor, model_name, loss_function, ".txt", sep = "_", collapse = NULL), row.names = FALSE)

