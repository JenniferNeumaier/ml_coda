## Importing imputed data sets
setwd("D:/Master/Master Biologie/Sommersemester 2022/Master thesis/Github")

# Importing metadata, taxonomy, and otu_table (!very big) from PCOS data set
PCOS_abundances_10_pc <- read.table("./data/processed/PCOS_abundances_10_pc.txt", header = TRUE)
PCOS_abundances_50_pc <- read.table("./data/processed/PCOS_abundances_50_pc.txt", header = TRUE)

PCOS_metadata <- read.table("./data/processed/PCOS_metadata.txt", header = TRUE)

# colorectal cancer (CRC) data set (large file, therefore ignored until needed)
CRC_abundances_10_pc <- read.table("./data/processed/CRC_abundances_10_pc.txt", header = TRUE)
CRC_abundances_50_pc <- read.table("./data/processed/CRC_abundances_50_pc.txt", header = TRUE)

CRC_metadata <- read.table("./data/processed/CRC_metadata.txt", header = TRUE)

# EstMB data set (large file, therefore ignored until needed)
EstMB_abundances_10_pc <- read.table("./data/processed/EstMB_abundances_10_pc.txt", header = TRUE)
EstMB_abundances_50_pc <- read.table("./data/processed/EstMB_abundances_50_pc.txt", header = TRUE)

EstMB_metadata <- read.table("./data/processed/EstMB_metadata.txt", header = TRUE)