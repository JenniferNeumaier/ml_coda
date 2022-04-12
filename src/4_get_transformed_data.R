## Importing transformed data sets

# PCOS
PCOS_10_TSS <- read.table("./data/transformed/pseudocount/PCOS_10_pc_TSS.txt", header = TRUE)
PCOS_10_CLR <- read.table("./data/transformed/pseudocount/PCOS_10_pc_CLR.txt", header = TRUE)
PCOS_10_ALR_optimal <- read.table("./data/transformed/pseudocount/PCOS_10_pc_ALR_optimal.txt", header = TRUE)
PCOS_10_ALR_worst <- read.table("./data/transformed/pseudocount/PCOS_10_pc_ALR_worst.txt", header = TRUE)
PCOS_10_ALR_random <- read.table("./data/transformed/pseudocount/PCOS_10_pc_ALR_random.txt", header = TRUE)

PCOS_50_TSS <- read.table("./data/transformed/pseudocount/PCOS_50_pc_TSS.txt", header = TRUE)
PCOS_50_CLR <- read.table("./data/transformed/pseudocount/PCOS_50_pc_CLR.txt", header = TRUE)
PCOS_50_ALR_optimal <- read.table("./data/transformed/pseudocount/PCOS_50_pc_ALR_optimal.txt", header = TRUE)
PCOS_50_ALR_worst <- read.table("./data/transformed/pseudocount/PCOS_50_pc_ALR_worst.txt", header = TRUE)
PCOS_50_ALR_random <- read.table("./data/transformed/pseudocount/PCOS_50_pc_ALR_random.txt", header = TRUE)

PCOS_metadata <- read.table("./data/processed/PCOS_metadata.txt", header = TRUE)

CRC_metadata <- read.table("./data/processed/CRC_metadata.txt", header = TRUE)