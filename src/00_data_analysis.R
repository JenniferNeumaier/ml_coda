library("easyCODA")


# modified from https://github.com/nhanhocu/metamicrobiomeR/blob/master/R/taxa.filter.R
taxa.filter<-function(taxtab, percent.filter=percent.filter, relabund.filter=relabund.filter){
  # filter (remove) taxa with relative abundance <relabund.filter and available in <percent.filter of number of samples
  taxdat <- as.data.frame(apply(taxtab[,2:ncol(taxtab)], 2, function(x) as.numeric(as.character(x))))
  # get assigned taxa only
  taxlist<-colnames(taxdat[,2:ncol(taxdat)])
  # filter using percent.filter
  taxtest <- apply(taxdat[,taxlist],2,function(x){length(x[!is.na(x)&x>0])})
  taxget<-taxtest[taxtest>= as.numeric(percent.filter)*(nrow(taxdat))]
  #filter using relabund.filter
  taxtestm<-apply(taxdat[,taxlist],2,mean,na.rm=T)
  taxgetm<-taxtestm[taxtestm>as.numeric(relabund.filter)]
  taxlistf<-c(names(taxget)[names(taxget) %in% names(taxgetm)])
  # make new processed data frame
  dataframe <- cbind(taxtab[1], subset(taxdat, select = c(taxlistf)))
  return(dataframe)
}


my_TSS <- function(data){
  res <- data %>% 
    mutate_if(is.numeric, function(x) (x-min(x))/(max(x)-min(x)))
  return(res)
}


my_CLR <- function(data){
  idx <- data[1]
  features <- data[2:ncol(data)]
  CLR_list <- CLR(features)
  res <- cbind(idx, as.data.frame(CLR_list$LR))
  return(res)
}


my_random_ALR <- function(data){
  idx <- data[1]
  features <- data[2:ncol(data)]
  denominator <- sample(1:length(features), 1)
  ALR_list <- ALR(features, denom = denominator)
  res <- cbind(idx, as.data.frame(ALR_list$LR))
  return(res)
}

my_optimal_ALR <- function(data){
  # split data into IDs and features
  idx <- data[1]
  features <- data[2:ncol(data)]
  
  # find good alr references
  alr.refs <- FINDALR(features)
  
  # top 20 for correlation
  names(alr.refs$procrust.cor) <- colnames(features)
  res_procrustes <- alr.refs$procrust.cor[order(alr.refs$procrust.cor, decreasing=TRUE)][1:20]
  
  # top 20 for log variance
  names(alr.refs$var.log) <- colnames(features)
  res_variance <- alr.refs$var.log[order(alr.refs$var.log)][1:20]
  
  # find best overlap
  diff <- intersect(names(res_procrustes), names(res_variance))
  denominator <- diff[1]
  
  # if no overlap can be found, output res_procrustes and res_variance so that the user can decide manually
  if (is_empty(denominator) == TRUE){
    
    return(res_procrustes, res_variance)
    
  }
  else {
    
    ALR_list <- ALR(features, denom = which(colnames(features) == denominator))
    res <- cbind(idx, as.data.frame(ALR_list$LR))
    
    return(res)
  }
}


my_worst_ALR <- function(data){
  # split data into IDs and features
  idx <- data[1]
  features <- data[2:ncol(data)]
  
  # find good alr references
  alr.refs <- FINDALR(features)
  
  # last 20 for correlation (features with lowes correlation)
  names(alr.refs$procrust.cor) <- colnames(features)
  res_procrustes <- alr.refs$procrust.cor[order(alr.refs$procrust.cor)][1:20]
  
  # last 20 for log variance (features with highest variance)
  names(alr.refs$var.log) <- colnames(features)
  res_variance <- alr.refs$var.log[order(alr.refs$var.log, decreasing = TRUE)][1:20]
  
  # find best overlap
  diff <- intersect(names(res_procrustes), names(res_variance))
  denominator <- diff[1]
  
  # if no overlap can be found, output res_procrustes and res_variance so that the user can decide manually
  if (is_empty(denominator) == TRUE){
    
    return(res_procrustes, res_variance)
    
  }
  else {
    
    ALR_list <- ALR(features, denom = which(colnames(features) == denominator))
    res <- cbind(idx, as.data.frame(ALR_list$LR))
    
    return(res)
  }
}


FINDALR <- function(data, weight = FALSE) { 
  
  ### -------------------------------------------------------------------
  ### https://github.com/michaelgreenacre/CODAinPractice
  ###
  ### function to identify the best reference for a set of ALRs
  ### various statistics are computed for each reference to assist
  ### the choice of best reference
  ### data is a normalized data matrix
  ### equal weighting is default for the logratio geometry here
  ### row (sample) weighting not implemented in this version
  
  if(sum(data[1,]!=1)) data <- data/rowSums(data)
  
  ### first compute the exact logratio geometry
  data.lra <- LRA(data, weight=weight)
  data.lra.rpc <- data.lra$rowpcoord 
  tot.var <- sum(data.lra$sv^2)
  
  ### loop on all the potential references, computing Procrustes correlation
  ### of each set of ALRs with the exact geometry
  procrust.cor <- rep(0, ncol(data))
  dim <- min(nrow(data), ncol(data)) - 1
  for(j in 1:ncol(data)) {
    # ALR transformation
    alr <- ALR(data, denom=j, weight=weight)
    # ALR geometry using PCA 'by hand' using SVD, without or with weighting
    if(!weight) {
      alr.svd <- svd(sqrt(1/nrow(alr$LR)) * sweep(alr$LR, 2, colMeans(alr$LR)) * sqrt(1/ncol(alr$LR)))
      alr.rpc <- sqrt(nrow(alr$LR)) * alr.svd$u %*% diag(alr.svd$d)
    }
    if(weight) {
      c <- colMeans(data)
      cc <- c*c[j]
      cc <- cc[-j]
      alr.svd <- svd(sqrt(1/nrow(alr$LR)) * sweep(alr$LR, 2, colMeans(alr$LR)) %*% diag(sqrt(cc)))
      alr.rpc <- sqrt(nrow(alr$LR)) * alr.svd$u %*% diag(alr.svd$d)
    }
    procrust.cor[j] <- protest(alr.rpc[,1:dim],data.lra.rpc, permutations=0)$t0
  }
  
  ### the variances of the log-transformed parts
  var.log <- as.numeric(apply(log(data), 2, var))
  
  ### highest Procrustes correlation
  procrust.max <- max(procrust.cor)
  
  ### which reference gives maximum Procrustes
  procrust.ref <- which(procrust.cor==procrust.max)
  
  ### lowest log variance
  var.min <- min(var.log)
  
  ### which reference gives lowest log variance
  var.ref <- which(var.log==var.min)
  
  return(list(tot.var=tot.var, procrust.cor=procrust.cor, 
              procrust.max=procrust.max, procrust.ref=procrust.ref,
              var.log=var.log, var.min=var.min, var.ref=var.ref))
}


transformation <- function(data = data, method = c("ALR_optimal", "ALR_worst", "ALR_random", "CLR", "TSS")) {
  
  
  ifelse(method == "ALR_optimal", res <- my_optimal_ALR(data), 
         ifelse(method == "ALR_worst", res <- my_worst_ALR(data), 
                ifelse(method == "ALR_random",  res <- my_random_ALR(data),
                       ifelse(method == "CLR",  res <- my_CLR(data), 
                              ifelse(method == "TSS", res <- my_TSS(data))
                       )
                )
         )
  )
  return(res)
}