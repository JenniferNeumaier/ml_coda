## Importing raw data sets

# Importing metadata, taxonomy, and otu_table (!very big) from PCOS data set
PCOS_abundances <- read.delim("./data/raw/PCOS_abundances.txt")
PCOS_metadata <- read.csv("./data/raw/PCOS_metadata.csv")
PCOS_taxonomy <- read.csv("./data/raw/PCOS_taxonomy.csv")

# Keeping only phenotype from metadata
PCOS_metadata <- read.csv("./data/raw/PCOS_metadata.csv")

# colorectal cancer (CRC) data set (large file, therefore ignored until needed)
CRC_abundances <- read.csv("./data/raw/CRC_abundances.csv")


# EstMB data set (large file, therefore ignored until needed)
EstMB_abundances <- read.csv("./data/raw/EstMB_abundances.csv")
EstMB_metadata <- read.csv("./data/raw/EstMB_metadata.csv")

