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
file_path <- "../data/processed/filtered/EstMB_abundances_90.txt"

# Declare file path for metadata
metadata_file <- "../data/EstMB/EstMB_metadata.txt"

# Declare Predictor
predictor <- "E11" 

# Number of repeats per data set
n_repeats <- 10

# Declare type of loss function
loss_function <- "AUC" # AUC or RMSE

# Declare number of kfolds
kfold <- 5

# Declare number if CV times
cv_times <- 10

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
library("codacore")

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


# codacore ----
#------------------------------------------------#
#                                                #
#               CODACORE                         # 
#                                                #
#------------------------------------------------#
counter <- seq(1,n_repeats,1)

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
        
        
        coda_test[[j]] <- c(as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_0_lr1, quiet = T)), 
                            as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_0, quiet = T)), 
                            as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_1_lr1, quiet = T)),
                            as.numeric(ModelMetrics::rmse(actual = y_test, predicted = pred_bal_1, quiet = T)))
      }
    )
  }
  
}

# output ----
#------------------------------------------------#
#                                                #
#               OUTPUT                           # 
#                                                #
#------------------------------------------------#


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
write.table(final, paste("../out/EstMB_90_codacore_simputation", predictor, loss_function, ".txt", sep = "_", collapse = NULL), row.names = FALSE)

