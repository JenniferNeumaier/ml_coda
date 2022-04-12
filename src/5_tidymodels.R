## Linear Regression (ENET)

## Preprocessing

# step_dummy: tell the script to not use SampleID as predictor
data_recipe <- recipe(PCOS_Riikka ~ ., data = train_data) %>% 
  step_dummy(SampleID)

# create workflow
wf <- workflow() %>% 
  add_recipe(data_recipe)


# Tuning ENET parameters
# how to pick penalty?
set.seed(123)

data_boots <- bootstraps(train_data, strata = PCOS_Riikka)

tune_spec <- linear_reg(penalty = tune(), mixture = 0.5) %>% 
  set_engine("glmnet")

lambda_grid <- grid_regular(
  penalty(),
  levels = 50)


doParallel::registerDoParallel()

set.seed(2020)

enet_grid <- tune_grid(
  wf %>% add_model(tune_spec),
  resamples = data_boots,
  grid = lambda_grid
)

# observe Parameter tuning

enet_grid %>% 
  collect_metrics %>% 
  ggplot(aes(penalty, mean, color = .metric)) +
  geom_errorbar(aes(ymin = mean - std_err,
                    ymax = mean + std_err),
                alpha = 0.5) +
  geom_line(size = 1.5, show.legend = FALSE) + 
  facet_wrap(~.metric, scales = "free", nrow = 2) + 
  scale_x_log10() + 
  theme(legend.position = "none")


# Use best parameters
lowest_RMSE <- enet_grid %>% 
  select_best("rmse", maximize = FALSE) 
# maximize = FALSE: lowest RMSE
# maximize = TRUE: highest RMSE

final_enet <- finalize_workflow(wf %>% add_model(tune_spec),
                                lowest_RMSE)

final_enet %>% 
  fit(train_data) %>% 
  pull_workflow_fit() %>% 
  vi(lambda = lowest_RMSE$penalty) %>% 
  mutate(Importance = abs(Importance),
         Variable = fct_reorder(Variable, Importance)) %>% 
  ggplot(aes(x = Importance, y = Variable, fill = Sign)) + 
  geom_col() + 
  scale_x_continuous(expand = c(0,0)) + 
  labs(y = NULL)


# Model specifications
# ENET = mixture = 0.5
# LASSO = mixture = 1
# Ridge = micture = 0

# set_enginge: which one?
enet_spec <- linear_reg(penalty = 0.1, mixture = 0.5) %>% 
  set_engine("glmnet")

wf <- workflow() %>% 
  add_recipe(data_recipe)

enet_fit <- wf %>% 
  add_model(enet_spec) %>% 
  fit(data = train_data)

enet_fit %>% 
  pull_workflow_fit() %>% 
  tidy()


## XGBoost
# https://www.r-bloggers.com/2020/05/using-xgboost-with-tidymodels/

# XGBoost model specification
xgboost_model <- boost_tree(
  trees = 1000,
  tree_depth = tune(), min_n = tune(), loss_reduction = tune(),
  sample_size = tune(), mtry = tune(),
  learn_rate = tune(),
) %>%
  set_engine("xgboost") %>% 
  set_mode("classification")

xgboost_model

# grid specification
xgboost_grid <- grid_latin_hypercube(
  min_n(),
  tree_depth(),
  learn_rate(),
  loss_reduction(),
  finalize(mtry(), train_data), #mtry() has unknown  
  sample_size = sample_prop(),
  size = 20
)

xgboost_grid


xgb_wf <- workflow() %>% 
  add_formula(data_set$y ~ .) %>% 
  add_model(xgboost_model)

set.seed(123)
vb_folds<- vfold_cv(train_data, strata = train_data$y)

doParallel::registerDoParallel()

set.seed(123)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = vb_folds,
  grid = xgboost_grid,
  control = control_grid(save_pred = TRUE)
)


## LDA

# Creating recipe and preprocessing

lda_recipe <- recipe(y ~ ., data = train_data) %>% 
  step_normalize(all_numeric(), -all_outcomes())

# Create Model
# LDA = frac_common_cov = 1
# QDA = frac_common_cov = 0

lda_model <- discrim_regularized(frac_common_cov = 1) %>% 
  set_engine("klaR") %>% 
  set_mode("classification")

lda_wf <- workflow() %>% 
  add_model(lda_model) %>% 
  add_recipe(lda_recipe)

last_fit_lda <- lda_wf %>% 
  last_fit(split = data_split)

last_fit_lda %>% collect_metrics()

lda_predictions <- last_fit_lda %>% 
  collect_predictions()

lda_predictions
