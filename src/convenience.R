sigmoid <- function(x) {
  return(1/(1 + exp(-x)))
}
 

logit <- function(x, w) {
  random_number <- runif(length(x), min  = 0, max = 1)
  return(sample(random_number, 1) <= sigmoid((x - 0.5) %*% w))
}


data_sampler <- function(n, k, f){ 
  
  #     Data generator that generates n x k feature matrix and a target vector
  #     
  #     Returns a data frame with columns x1, ..., xk from randomised uniform distribution.
  #     y is computed using labelling function f.
  dataframe <- data.frame(matrix(runif(n*k, min = 0, max = 1), n, k))
  dataframe <- sapply(dataframe, function(x) {x >= 0.5})
  # dataframe[ ,ncol(dataframe) + 1] <- apply(dataframe, 1, function(x) f)
  
  return(dataframe)
}


microbiome_sampler <- function(n,k){
  columns <- paste(rep("OTU", k), seq(1,k,1), sep = "_")
  dataframe <- data.frame(round(matrix(runif(n*k, min = 1, max = 2000), n, k)))
  colnames(dataframe) <- columns
  dataframe$y <- rnorm(n)
  dataframe$y <- ifelse(dataframe$y > 0, 1, 0)
  dataframe <- dataframe %>% 
    relocate(y)
  return(dataframe)
} 

# Imports several data objects into a list object
import <- function(path = path, pattern = pattern, header = header) {
  # import cluster data
  temp <- list.files(path = path, pattern = pattern, full.names = TRUE)
  data_list <- lapply(temp, function(x) read.table(x, header = header))
  
  # extract cluster names and name data_list
  data_names <- lapply(temp, function(x) sub("\\.[[:alnum:]]+$", "", basename(as.character(x))))
  names(data_list) <- data_names
  return(data_list)
}

outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}





