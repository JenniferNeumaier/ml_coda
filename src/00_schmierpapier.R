# combination of output tables

tables_codacore <- import(path = "../out/EstMB/I11/svmRadial/", pattern = ".*(codacore).*", header = TRUE)
tables_codacore <- do.call(rbind, tables_codacore) %>% 
  add_column(set = "codacore")
rownames(tables_codacore) <- NULL

tables_mikrop <- import(path = "../out/EstMB/I11/svmRadial/", pattern = ".*(mikrop).*", header = TRUE)
tables_mikrop <- do.call(rbind, tables_mikrop) %>% 
  add_column(set = "svmRadial") %>% 
  rename("model" = transformation)
rownames(tables_mikrop) <- NULL

table_final <- rbind(tables_codacore, tables_mikrop)

write.table(table_final, "../out/EstMB/EstMB_I11_model_summary.txt", row.names = FALSE)


# illustrative example compositional data

animals <- c("bee", "swallow", "rabbit", "wolf")
counts <- c(8, 7, 4, 1)

abundance_table_A <- data.frame(animals, counts, field = "A")
abundance_table_B <- data.frame(animals, counts = ceiling(counts/2), field = "B")

abundance_table <- rbind(abundance_table_A, abundance_table_B)

ggplot(abundance_table, aes(field, counts, fill = animals)) + 
  geom_bar(position = "stack", stat = "identity") + 
  labs(title = "Absolute Counts", x = "Samples", y = "Counts")

ggplot(abundance_table, aes(field, counts, fill = animals)) + 
  geom_bar(position = "fill", stat = "identity") + 
  labs(title = "Normalized Counts", x = "Samples", y = "Percentages")

# Ternary plots
library("Ternary")
library("dplyr")
library("ggplot2")

# Load src documents
source("./convenience.R")
source("./data_analysis.R")

# Import data
EstMB_data <- read.table("../data/EstMB/EstMB_abundances_90.txt", header = TRUE)
metadata <- read.table("../data/EstMB/EstMB_metadata.txt", header = TRUE)

metadata <- metadata %>% 
  mutate(E11 = ifelse(E11 == "1", "DT2", "healthy")) %>% 
  mutate(F41 = ifelse(F41 == "1", "PD", "healthy")) %>% 
  mutate(I11 = ifelse(I11 == "1", "HHD", "healthy"))

# set variables
merge_id <- names(metadata[1])
leftover_factors <- outersect(c(merge_id, "E11"), names(metadata))

# merge data used for mikrop models

EstMB_data <- EstMB_data %>% 
  select(0:4) %>% 
  left_join(metadata, by = merge_id) %>% 
  dplyr::relocate(c("E11", leftover_factors), .after = merge_id) %>% 
  dplyr::select(-c(merge_id, leftover_factors)) 

EstMB_DT2 <- subset(EstMB_data, E11 == "healthy")

# TernaryPlot
TernaryPlot(point = "up", 
            atip = "Heimdallarchaeota.archaeon", 
            btip = "Saccharolobus.solfataricus", 
            ctip = "Halosimplex.carlsbadense",
            axis.labels = seq(0, 500, by = 1),
            isometric = TRUE)
TernaryPoints(EstMB_DT2[2:4], pch = 16)

# ALR transformation
ALR_list <- ALR(EstMB_DT2[2:4], denom = 3)
EstMB_DT2_ALR <- cbind(EstMB_DT2[1], ALR_list$LR)

# Linear Relationship after Transformation
EstMB_DT2_ALR %>% 
  rename("Heimdallarchaeota.archaeon" = colnames(EstMB_DT2_ALR)[2],
         "Halosimplex.carlsbadense" = colnames(EstMB_DT2_ALR)[3]) %>% 
  ggplot(aes(x = Heimdallarchaeota.archaeon, y = Halosimplex.carlsbadense)) + 
  geom_point() + 
  labs(x = "Heimdallarchaeota.archaeon/Saccharolobus.solfataricus", y = "Halosimplex.carlsbadense/Saccharolobus.solfataricus")


# Checking on Correlation Coefficients between CRC_abundances_10 and CRC_abundances_50
CRC_all <- read.table("../data/processed/filtered/CRC_abundances_10.txt", header = TRUE) %>% 
  dplyr::select(-SampleID) %>% 
  as.matrix()

# Load the sub-compositional data 
CRC_sub_comp <- read.table("../data/CRC/CRC_abundances_10.txt", header = TRUE) %>% 
  dplyr::select(-SampleID) %>% 
  as.matrix()

# Get the correlation 
cor.all <- as.vector(cor(t(CRC_all[1:567,])))
cor.sub.comp <- as.vector(cor(t(CRC_sub_comp[1:567,])))

tmp <- as.data.frame(cbind(cor.all,cor.sub.comp))
names(tmp) <- c('correlation_all', 'correlation_sub_comp')
tmp$abs.diff <- as.factor(ifelse(abs(tmp$correlation_all - tmp$correlation_sub_comp)>0.5,1,0))

# Plot
ggplot(tmp,aes(correlation_all,correlation_sub_comp, color=abs.diff)) + 
  geom_point(size=2) + 
  labs(title = "CRC: Subcompositional Coherence", x = "Corr_before_Imputation", y = "Corr_after_Imputation") + 
  scale_colour_manual(values = c("1" = "Red", "0" = "Blue")) + 
  theme(legend.position = "none")

# Tables; calculate mean + sd for all groups
temp <- read.table("../out/EstMB/EstMB_data_leakage.txt", header = TRUE) %>% 
  # filter(filtering == "90", set == "svmRadial") %>% 
  mutate(dplyr::across(model, factor, levels=c("model_bal_0",
                                               "model_bal_1",
                                               "pred_bal_0",
                                               "pred_bal_0_lr1",
                                               "pred_bal_1",
                                               "pred_bal_1_lr1",
                                               "ALRo",
                                               "ALRr",
                                               "ALRw",
                                               "CLR",
                                               "TSS")))


ag <- aggregate(loss ~ ., temp, function(x) c(mean = mean(x), sd = sd(x)))
