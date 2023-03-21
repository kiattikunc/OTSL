ebic <- function(node, parents, data, args) {
  
  n = nrow(data)
  N = length(data)
  gamma_penalty=args$gamma_penalty
  gamma_penalty=as.numeric(gamma_penalty)
  if (length(parents) == 0) {
    
    counts = table(data[, node])
    nparams = dim(counts) - 1
    a1 = sum(counts * log(counts / n),na.rm=TRUE) - nparams * log(n) / 2 - gamma_penalty *nparams* log(N)
    return(a1)
    # 
  }else if (length(unique(data[, node])) == 1) {
    
    counts = table(data[, node])
    nparams = dim(counts) - 1
    a3 = sum(counts * log(counts / n),na.rm=TRUE) - nparams * log(n) / 2 - gamma_penalty *nparams* log(N)
    return(a3)
    # 
  }
  else {
    
    counts = table(data[, node], configs(data[, parents, drop = FALSE]))
    nparams = ncol(counts) * (nrow(counts) - 1)
    a2 = sum(counts * log(prop.table(counts, 2)),na.rm=TRUE)  - nparams * log(n) / 2 - gamma_penalty *nparams* log(N) 
    return(a2)
    
    
  }#ELSE
}