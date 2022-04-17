# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#


# Import libraries
library("dplyr")
library("ggplot2")

# Import files

table_models <- read.table("../out/table_codacore_BMI.txt", header = TRUE)

facet_labels <- c(
   "10" = "10% abundance filter",
   "50" = "50% abundance filter"
)

# for mikropml
graph <- ggplot(table_models, aes(x = transformation, y = loss, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~filtering, labeller = as_labeller(facet_labels), scale="free") + 
  labs(title = "GLM for PCOS data", x = "Transformationen", y = "AUC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph

# for codacore
graph <- ggplot(table_models, aes(x = model, y = loss, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~filtering, labeller = as_labeller(facet_labels), scale="free") + 
  labs(title = "Codacore for PCOS data", x = "Codacore models", y = "RMSE") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
graph


ggsave(plot = graph, filename = "./figs/figure.jpeg", device = "jpeg", height = 6, width = 7)

