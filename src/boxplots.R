# Boxplots for models ----
#------------------------------------------------#
#                                                #
#              BOXPLOTS FOR MODELS               # 
#                                                #
#------------------------------------------------#


# Import libraries
library("dplyr")
library("tidyverse")
library("ggplot2")

# Load src documents
source("./convenience.R")
source("./data_analysis.R")

# Import files
table_models <- read.table("../out/CRC_holdout/table_codacore_GER_Group_glmnet_AUC_.txt", header = TRUE)

facet_labels <- c(
   "10" = "10% abundance filter",
   "50" = "50% abundance filter"
)

# for mikropml
graph <- ggplot(table_models, aes(x = transformation, y = loss, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~filtering, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "GLM for CRC data", x = "Transformations", y = "AUC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph

# for codacore
graph <- ggplot(table_models, aes(x = model, y = loss, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~filtering, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "Codacore for CRC data", x = "Codacore", y = "AUC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph


ggsave(plot = graph, filename = "./figs/figure.jpeg", device = "jpeg", height = 6, width = 7)


# Boxplots for ALR comparison ----
#------------------------------------------------#
#                                                #
#              BOXPLOTS FOR ALR                  # 
#                                                #
#------------------------------------------------#

tables <- import(path = "../out/leaky_nonleaky", pattern = ".*(PCOS).*", header = TRUE)

df <- data.frame()
for(i in 1:length(tables)){
  res <- tables[[i]] %>% 
    add_column(dataset = names(tables)[[i]]) %>% 
    separate(dataset, c("dataset", "preprocessing", "country"), "_") %>% 
    dplyr::select(-c("dataset"))
  
  df <- rbind(df, res)
}


graph <- ggplot(df, aes(x = transformation, y = loss, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~preprocessing, scale="free_x") + 
  labs(title = "Influence of Transformation on Data Leakage", x = "Transformations", y = "AUROC", ylim = c(0.40, 0.60)) +
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

graph <- ggplot(df_combined, aes(x = transformation, y = loss, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~set, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "Influence of Transformation on Holdout vs. 80/20", x = "Transformations", y = "AUROC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph


# subset 10% abundance filter and plot again

df_subset <- subset(df_combined, df_combined$filtering == "10")

facet_labels <- c(
  "Holdout" = "Holdout",
  "Standard" = "Standard 80/20"
)

graph <- ggplot(df_subset, aes(x = transformation, y = loss, col = type)) + 
  geom_boxplot() + 
  facet_wrap(~set, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "Influence of Transformation on Holdout vs. 80/20", x = "Transformations", y = "AUROC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph




# get mean of all data that is not 80/20
df_mean <- df_combined %>% 
  filter(!country == "80/20") %>% 
  group_by(filtering, transformation, type) %>% 
  summarise(mean = mean(loss))

facet_labels <- c(
  "10" = "10% abundance filter",
  "50" = "50% abundance filter"
)

graph <- ggplot(df_mean, aes(x = transformation, y = mean, fill = country)) + 
  geom_boxplot() + 
  facet_wrap(~type, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "Influence of Transformation on Holdout vs. 80/20", x = "Transformations", y = "AUROC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph

