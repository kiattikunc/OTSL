average_model_MI <- function(learnt_cpdag,data_train){
  library(bnlearn)
  library(BNSL)
  MI_val <- c()
  MI_bnsl<- c()
  MI_sh<-c()
  #learnt_cpdag<- hc_runner 
  #data_train <- dataset_test
  all_edges <- learnt_cpdag$arcs
  data_train_factor <- data_train
  
  if (nrow(all_edges)>0)
  {
    cols <- colnames(data_train)
    for(j in cols){
      data_train_factor[[j]]<- as.numeric(data_train[[j]])
      data_train_factor[[j]]<- data_train_factor[[j]]-1
    }
  
    for (edge_index in 1: nrow(all_edges))
    {
      MI_val <- rbind(MI_val,ci.test(all_edges[edge_index,1],all_edges[edge_index,2], data = data_train, test = "mi")$statistic)
      MI_bnsl <- rbind(MI_bnsl, mi(sapply(data_train_factor[all_edges[edge_index,1]], as.numeric), sapply(data_train_factor[all_edges[edge_index,2]], as.numeric), proc=2))
      MI_sh <- rbind(MI_sh,ci.test(all_edges[edge_index,1],all_edges[edge_index,2], data = data_train, test = "mi-sh")$statistic)
      
       }
  }
  else{
    MI_val = NA
    MI_bnsl =NA
    MI_sh =NA
  }
  return_list <- list('avg_MI'=mean(MI_val,na.rm=TRUE),'avg_MI_sh'=mean(MI_sh,na.rm=TRUE),'total_MI'=Reduce('+', MI_val),'total_bnslMI'=Reduce('+', MI_bnsl),'total_MI_sh'=Reduce('+', MI_sh))
}