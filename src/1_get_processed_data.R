## Importing processed data sets

# Importing metadata, taxonomy, and otu_table (!very big) from PCOS data set
PCOS_abundances_10 <- read.csv("./data/processed/PCOS_abundances_10.csv")
PCOS_abundances_50 <- read.csv("./data/processed/PCOS_abundances_50.csv")

PCOS_metadata <- read.csv("./data/processed/PCOS_metadata.csv")

# colorectal cancer (CRC) data set (large file, therefore ignored until needed)
# CRC_abundances <- read.csv("./data/processed/CRC_abundances.csv")


# EstMB data set (large file, therefore ignored until needed)
EstMB_abundances_10 <- read.csv("./data/processed/EstMB_abundances_10.csv")
EstMB_abundances_50 <- read.csv("./data/processed/EstMB_abundances_50.csv")


EstMB_metadata <- read.csv("./data/raw/EstMB_metadata.csv")
