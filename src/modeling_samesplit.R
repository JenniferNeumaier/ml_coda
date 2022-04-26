# TODO:
# - rename codacore models in something more distinct
# - check codacore for RMSE and add to for-loop (done)
# - add message to tryCatch()
# - add args for cvtimes and kfold (not as important)



# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)             
# args = c("RData/Abundance_filtered_species_p10a001.rds", "RData/Data_master.rds", "LASSO", "CLR", "1", "5", "SET2", "Rdata/ML", "F41", "I10") 
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
n_repeats <- args[4]

if (n_repeats == 0) {
  stop("Counter has to be set to at least 1!", call. = FALSE)
}

# Declare model type
model_name <- args[5] #glmnet or xgbTree (xgboost)

# Declare type of loss function
loss_function <- args[6] # AUC or RMSE

# Declare training fraction
training_frac <- as.numeric(args[7])

# Declare number of kfolds
# kfolds <- args[7]

# Declare number if CV times

# cv_times <- args[8]


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
library("codacore")


# Load src documents
source("./convenience.R")
source("./data_analysis.R")


# Import data
data_transformed <- import(path = file_path, pattern = ".*(transform_*).*", header = TRUE)
data_imputed <- import(path = file_path, pattern=".*(abundances_.*).*", header = TRUE)
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
leftover_factors <- outersect(c(merge_id, predictor), names(metadata))


# merge data used for mikrop models
for (i in 1:length(data_transformed)){
  data_transformed[[i]] <- data_transformed[[i]] %>% 
    left_join(metadata, by = merge_id) %>% 
    relocate(c(predictor, leftover_factors), .after = merge_id) %>% 
    dplyr::select(-c(merge_id, leftover_factors))
}


# merge data used for codacore and also split
train_data <- list()
test_data <- list()

for (i in 1:length(data_imputed)){
  
  # merge with metadata
  data_imputed[[i]] <- data_imputed[[i]] %>% 
    left_join(metadata, by = merge_id) %>% 
    relocate(c(predictor, leftover_factors), .after = merge_id) %>% 
    dplyr::select(-c(merge_id, leftover_factors))
  
  # Split the data
  set.seed(2022)
  data_split <- initial_split(data_imputed[[i]], prop = training_frac, strata = predictor) 
  train_data[[i]] <- training(data_split) %>% as.data.frame()
  test_data[[i]] <- testing(data_split) %>% as.data.frame()
  
}

# name train and test data
names(train_data) <- names(data_imputed)
names(test_data) <- names(data_imputed)


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
                        training_frac = idx, # training fraction taken from split for codacore
                        seed = NA)
    
    # save results in temp lists
    if(loss_function == "AUC"){
      
      training_res[j] <- ml_result$performance[1]
      test_res[j] <- ml_result$performance[3]
      
    } else if (loss_function == "RMSE"){
      
      training_res[j] <- ml_result$performance[1]
      test_res[j] <- ml_result$performance[2]
      
    }
 
    paste("Round", j, sep = " ", collapse = NULL)
    
  }
  
  # save lists in final list
  mikrop_training[[i]] <- training_res
  mikrop_test[[i]] <- test_res
  
  
  paste("Finished", data_transformed[[i]], sep = " ", collapse = NULL)
}


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


for (i in 1:length(train_data)) { 
  
  # generate features (x) and variable (y) for train
  x_train <- train_data[[i]]
  x_train <- x_train[, 2:ncol(x_train)]
  
  y_train <- train_data[[i]]
  y_train <- y_train[,1] 
  
  # generate features (x) and variable (y) for test
  x_test <- test_data[[i]]
  x_test <- x_test[, 2:ncol(x_test)]
  
  y_test <- test_data[[i]]
  y_test <- y_test[,1] 
  
  # initiate temp lists
  training_res <- list()
  test_res <- list()
  
  for (j in counter) {
    
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
  
  paste("Finished", train_data[[i]], sep = " ", collapse = NULL)
  
  }
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

write.table(final, paste("../out/table_mikrop", predictor, model_name, loss_function, ".txt", sep = "_", collapse = NULL), row.names = FALSE)


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

