# Boxplots for models ----
#------------------------------------------------#
#                                                #
#              BOXPLOTS FOR MODELS               # 
#                                                #
#------------------------------------------------#


# Import libraries
library("dplyr")
library("ggplot2")

# Load src documents
source("./src/convenience.R")
source("./src/data_analysis.R")

# Import files
table_models <- read.table("../out/table_codacore_BMI.txt", header = TRUE)

facet_labels <- c(
   "10" = "10% abundance filter",
   "50" = "50% abundance filter"
)

# for mikropml
graph <- ggplot(table_models, aes(x = transformation, y = loss, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~filtering, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "GLM for PCOS data", x = "Transformationen", y = "AUC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph

# for codacore
graph <- ggplot(table_models, aes(x = model, y = loss, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~filtering, labeller = as_labeller(facet_labels), scale="free_x") + 
  labs(title = "Codacore for PCOS data", x = "Codacore models", y = "RMSE") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph


ggsave(plot = graph, filename = "./figs/figure.jpeg", device = "jpeg", height = 6, width = 7)


# Boxplots for ALR comparison ----
#------------------------------------------------#
#                                                #
#              BOXPLOTS FOR ALR                  # 
#                                                #
#------------------------------------------------#

tables <- import(path = "./out", pattern = ".*(PCOS*).*", header = TRUE)

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
  labs(title = "Influence of Transformation on Preprocessing", x = "Transformations", y = "AUC", ylim = c(0.40, 0.60)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph

