
# Import libraries
library("dplyr")
library("tidyverse")
library("ggplot2")

# Load src documents
source("./convenience.R")
source("./data_analysis.R")


# Boxplots for models ----
#------------------------------------------------#
#                                                #
#              BOXPLOTS FOR MODELS               # 
#                                                #
#------------------------------------------------#

# Import files
table_models <- read.table("../out/EstMB/EstMB_filtering.txt", header = TRUE) %>% 
  # filter(!set == "non_leaky_special", dataset == "I11") %>% 
  filter(dataset == "I11") %>% 
  # filter(filtering == "10", dataset == "glmnet")
  mutate(dplyr::across(model, factor, levels=c("model_bal_0",
                                        "model_bal_1",
                                        "pred_bal_0",
                                        "pred_bal_0_lr1",
                                        "pred_bal_1",
                                        "pred_bal_1_lr1",
                                        #"ALRo",
                                        "ALRr",
                                        #"ALRw",
                                        "CLR",
                                        "TSS")))

# facet_labels <- c(
#    "leaky" = "Leaky",
#    "non_leaky" = "non-Leaky"
# )

facet_labels <- c(
  "90" = "90% prevalence filtering",
  "95" = "95% prevalence filtering"
)

# facet_labels <- c(
#   "codacore" = "CODACORE",
#   "glmnet" = "GLMNET",
#   "svmRadial" = "SVMRADIAL"
# )


# plot
graph <- ggplot(table_models, aes(x = model, y = loss, col = type)) + 
  geom_boxplot() + 
  facet_wrap(~filtering, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "EstMB: Influence of Feature Size (HHD)", x = "Models & Transformations", y = "AUROC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph


# Boxplots for holdout matching ----
#------------------------------------------------#
#                                                #
#              BOXPLOTS FOR HOLDOUTS             # 
#                                                #
#------------------------------------------------#

# prepare tables from holdout modeling
tables <- import(path = "../out/CRC_holdout", pattern = ".*(mikrop*).*", header = TRUE)

df <- data.frame()
for(i in 1:length(tables)){
  res <- tables[[i]] %>% 
    add_column(set = "Holdout", temp = names(tables)[[i]]) %>% 
    separate(temp, c("temp", "filler1", "country", "predictor", "model", "filler2"), "_") %>% 
    dplyr::select(-c("temp", "filler1", "filler2"))
  
  df <- rbind(df, res)
}

# prepare table from standard modeling
table_standard <- read.table("../out/CRC/table_mikrop_Group_glmnet_AUC_.txt", header = TRUE)

table_standard <- table_standard %>% 
  add_column(set = "Standard", country = "80/20", predictor = "Group", model = "glmnet")

df_combined <- rbind(df, table_standard)

facet_labels <- c(
  "Holdout" = "Holdout",
  "Standard" = "Standard 80/20"
)

# boxplot
graph <- ggplot(df_combined, aes(x = transformation, y = loss, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~set, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "Influence of Transformation on Holdout vs. 80/20", x = "Transformations", y = "AUROC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph

# jitterplot
graph <- ggplot(df_combined, aes(x = transformation, y = loss, col = set, shape = type)) + 
  geom_jitter() + 
  # facet_wrap(~set, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "Influence of Transformation on Holdout vs. 80/20", x = "Transformations", y = "AUROC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph


