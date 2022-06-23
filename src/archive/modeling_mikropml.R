# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)             

# Change parameters here
args = c("../data/EstMB/EstMB_90_transform_ALRw_.txt", 
         "../data/EstMB/EstMB_abundances_90.txt", 
         "../data/EstMB/EstMB_metadata.txt", 
         "F41", 
         "5", 
         "glmnet", 
         "AUC", 
         "0.8", 
         "5",
         "5") 

cat(args, sep = "\n")


if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 


# Declare transformed file
transformed_file <- args[1]

# Declare imputed file
imputed_file <- args[2]

# Declare metadata file
metadata_file <- args[3]

# Declare predictor
predictor <- args[4]

# Number of repeats per data set
n_repeats <- as.numeric(args[5])
counter <- seq(1,n_repeats,1)

if (n_repeats == 0) {
  stop("Counter has to be set to at least 1!", call. = FALSE)
}

# Declare model type
model_name <- args[6] #glmnet or xgbTree (xgboost)

# Declare type of loss function
loss_function <- args[7] # AUC or RMSE

# Declare training fraction
training_frac <- as.numeric(args[8])

# Declare number of kfolds
kfolds <- as.numeric(args[9])

# Declare number if CV times
cv_times <- as.numeric(args[10])


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
data_transformed <- read.table(transformed_file, header = TRUE)
data_imputed <- read.table(imputed_file, header = TRUE)
metadata <- read.table(metadata_file, header = TRUE)

print("Imported data!")


# Data perparation ----
#------------------------------------------------#
#                                                #
#              DATA PREPARATION                  # 
#                                                #
#------------------------------------------------#

# change types of predictor columns 
# DT2: Diabetes tyoe 2
# PD: Panic disorcer
# HHD: Hypertensive heart disease
metadata <- metadata %>% 
  mutate(E11 = ifelse(E11 == "1", "DT2", "healthy")) %>% 
  mutate(F41 = ifelse(F41 == "1", "PD", "healthy")) %>% 
  mutate(I11 = ifelse(I11 == "1", "HHD", "healthy"))

# set variables
merge_id <- names(metadata[1])
leftover_factors <- outersect(c(merge_id, predictor), names(metadata))

# merge data used for mikrop models
data_transformed <- data_transformed %>% 
  left_join(metadata, by = merge_id) %>% 
  relocate(c(predictor, leftover_factors), .after = merge_id) %>% 
  dplyr::select(-c(merge_id, leftover_factors)) 

# merge with metadata for codacore models
data_imputed <- data_imputed %>% 
  left_join(metadata, by = merge_id) %>% 
  relocate(c(predictor, leftover_factors), .after = merge_id) %>% 
  dplyr::select(-c(merge_id, leftover_factors))

# Split the data
set.seed(2022)
data_split <- initial_split(data_imputed, prop = training_frac, strata = predictor) 
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

# repeat experiment n times
for (j in counter) { 
  
  ml_result <- run_ml(data_transformed,
                      model_name,
                      outcome_colname = predictor, 
                      kfold = kfolds,
                      cv_times = cv_times,
                      training_frac = idx, # training fraction taken from split for codacore
                      seed = NA)
  

  # save results in temp lists
  if(loss_function == "AUC"){
    
    mikrop_training[[j]] <- ml_result$performance[1]
    mikrop_test[[j]] <- ml_result$performance[3]
    
  } else if (loss_function == "RMSE"){
    
    mikrop_training[j] <- ml_result$performance[1]
    mikrop_test[j] <- ml_result$performance[2]
    
  }
  
}

# only save last model 
mikrop_models <- ml_result


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

  model_bal_0 <- codacore(
    x_train, 
    y_train, 
    objective = objective,
    logRatioType = "balances",
    cvParams = c(numFolds = kfolds),
    lambda = 0,
    verbose = FALSE, 
    fast = TRUE
  )
  
  model_bal_1 <- codacore(
    x_train, 
    y_train,
    objective = objective,
    logRatioType = "balances",  
    cvParams = c(numFolds = kfolds),
    lambda = 1,
    verbose = FALSE, 
    fast = TRUE    
  )
  
  if(loss_function == "AUC"){
    
    tryCatch(
      expr = {
        
        # save training AUC in lists
        coda_training[[j]] <- c(as.numeric(model_bal_0$ensemble[[1]]$AUC), 
                               as.numeric(model_bal_1$ensemble[[1]]$AUC))
        
        # get test AUC by predicting manually
        pred_bal_0_lr1 <- predict(model_bal_0, x_test, logits = FALSE, numLogRatios = 1)
        pred_bal_0 <- predict(model_bal_0, x_test, logits = FALSE)
        
        pred_bal_1_lr1 <- predict(model_bal_1, x_test, logits = FALSE, numLogRatios = 1)
        pred_bal_1 <- predict(model_bal_1, x_test, logits = FALSE)
        
        coda_test[[j]] <- c(as.numeric(pROC::auc(y_test, pred_bal_0_lr1, quiet = T)), 
                           as.numeric(pROC::auc(y_test, pred_bal_0, quiet = T)), 
                           as.numeric(pROC::auc(y_test, pred_bal_1_lr1, quiet = T)),
                           as.numeric(pROC::auc(y_test, pred_bal_1, quiet = T)))
      }
    )
    
  } else if (loss_function == "RMSE"){ 
    
    tryCatch(
      expr = { 
        
        # save training RMSE in lists
        coda_training[[j]] <- c(as.numeric(model_bal_0$ensemble[[1]]$RMSE), 
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
        
        
        coda_test[[j]] <- c(as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_0_lr1, quiet = T)), 
                           as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_0, quiet = T)), 
                           as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_1_lr1, quiet = T)),
                           as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_1, quiet = T)))
      }
    )
  }
  
}

# # only save last model 
# coda_models <- c(model_bal_0, model_bal_1)

# output ----
#------------------------------------------------#
#                                                #
#               OUTPUT                           # 
#                                                #
#------------------------------------------------#

## mikropml
training_df <- mikrop_training %>% 
  as.data.frame() %>% 
  stack() %>% 
  rename(loss = values) %>% 
  add_column(dataset = "EstMB", filtering = "90", transformation = transformations, type = "training") %>% 
  dplyr::select(-ind)


test_df <- mikrop_test %>% 
  as.data.frame() %>% 
  stack() %>% 
  rename(loss = values) %>% 
  add_column(dataset = "EstMB", filtering = "90", transformation = transformations, type = "test") %>% 
  dplyr::select(-ind)


# combine and outout
final <- rbind(training_df, test_df)
write.table(final, paste("../out/table_mikrop", predictor, model_name, loss_function, transformations, "90.txt", sep = "_", collapse = NULL), row.names = FALSE)


## codacore
training_df <- coda_training %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rename("model_bal_0" = V1, "model_bal_1" = V2) %>% 
  stack() %>% 
  rename("loss" = values, "model" = ind) %>% 
  add_column(dataset = "EstMB", filtering = "90", type = "training")

test_df <- coda_test %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rename("pred_bal_0_lr1" = V1, "pred_bal_0" = V2, "pred_bal_1_lr1" = V3, "pred_bal_1" = V4) %>% 
  stack() %>% 
  rename("loss" = values, "model" = ind) %>% 
  add_column(dataset = "EstMB", filtering = "90", type = "test")

# combine and output
final <- rbind(training_df, test_df)
write.table(final, paste("../out/table_codacore", predictor, model_name, loss_function, transformations, "90.txt", sep = "_", collapse = NULL), row.names = FALSE)

# # save models
# saveRDS(mikrop_models, paste("../out/models_mikrop", model_name, predictor, transformations, "90.rds", sep = "_", collapse = NULL))
# saveRDS(coda_models, paste("../out/models_codacore", model_name, predictor, transformations, "90.rds", sep = "_", collapse = NULL))

