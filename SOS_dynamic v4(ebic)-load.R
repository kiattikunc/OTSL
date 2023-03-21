# SOS_dynamic v3: fix reading configuration issue  as.numeric(as.character(all_config)) and compare the grpahical accuracy given pure obs data
# SOS_dynamic v3.1 add tuning score AIC and change format name 
# SOS_dynamic v3.2 add penalises and prior in to the consideration 
# SOS_dynamic v3.8.6 change to EBIC
# SOS_dynamic v3.8.6 penalised for both BDeu and BIC
# SOS_dynamic v3.8.6 rescale BIC to penalty [0,1]
# SOS_dynamic v3.8.6.1. use BIC from bnlearn
# SOS_dynamic v3.8.6.5. change to mcmc in MCMC
# SOS_dynamic v4 clean up all comments
# with network score name and resampling for the test and training data 


library(pcalg)
library(tictoc)
library(igraph)
library(bnlearn)
library(DOT)
library(BiDAG)
library(datasets)
library(caret)
source("precisionrecall_PAG.R")
source("average_KL_model.R")
source("graphutils.R")
source("average_MI_model.R")
source("ebic.R")
############### change configurations ################
iteration = 10 # number of bootstraping 

Rp = TRUE # bootstraping 
round = seq(1,iteration, by=1)
total_fold = iteration
Pd =seq(1,20, by=1)
Prior =seq(1,20, by=1)
alpha_list = c(0.05, 0.01, 0.1)
latent <- c("")
set_list = c('ASIA','ALARM','SPORTS','PROPERTY','hepar2','hailfinder')
none_value = rep(c('none'), times=20)
samplesize_list = c(1000,10000)
alpha_stat <- 0.05
s = 100

## select BDeu,BIC,AIC,LL,MI,MI-sh,avg_MI-sh, avg_MI tuning_score ## BDeu refers to BDeu-iss and BIC refers to EBIC-iss
tuning_score_list = c('BIC','BDeu')

#######################################################



for (tuning_score in tuning_score_list)
{
  
  
  ##setting for mcmc
  mcmc_config <- data.frame(algoname=rep('mcmc', times=20),
                            indtest=none_value,
                            score=rep(c('bdeu'), times=20),alpha=none_value,penalty=none_value,prior=Prior)
  ##setting for hc
  hc_config <- data.frame(algoname=rep('hc', times=40),
                          indtest=rep(c('none'), times=40),
                          score=c(rep(c('bdeu'), times=20),rep(c('bic'), times=20)),alpha=c(none_value,none_value),penalty=c(none_value,Pd),prior=c(Prior,none_value))
  ##setting for FGS_tetrad
  FGS_config <- data.frame(algoname=rep('FGS_tetrad', times=40),
                           indtest=rep(c('none'), times=40),
                           score=c(rep(c('bdeu'), times=20),rep(c('bic'), times=20)),alpha=c(none_value,none_value),penalty=c(none_value,Pd),prior=c(Prior,none_value))
  
  ##setting for pc
  pcstable_config <- data.frame(algoname=rep('pc', times=9),
                                indtest=c(rep(c('x2'), times=3),rep(c('mi'), times=3),rep(c('mi-sh'), times=3)),
                                score=rep(c('none'), times=9),alpha=alpha_list,penalty=rep(c('none'), times=9),prior=rep(c('none'), times=9))
  
  ##setting for mmhck=
  mmhc_config <- data.frame(algoname=rep('mmhc', times=40),
                            indtest=rep(c('x2'), times=40),
                            score=c(rep(c('bdeu'), times=20),rep(c('bic'), times=20)),alpha=c(rep(0.05, times=10),rep(0.01, times=10),rep(0.05, times=10),rep(0.01, times=10)),penalty=c(none_value,seq(1,10, by=1),seq(1,10, by=1)),prior=c(seq(1,10, by=1),seq(1,10, by=1),none_value))
  
  # write.table(mmhc_config, paste0("output/mmhc_config.csv"), sep=",", col.names=FALSE, row.names=FALSE)
  # write.table(pcstable_config, paste0("output/pc_config.csv"), sep=",", col.names=FALSE, row.names=FALSE)
  # write.table(FGS_config, paste0("output/fges_config.csv"), sep=",", col.names=FALSE, row.names=FALSE)
  # write.table(hc_config, paste0("output/hc_config.csv"), sep=",", col.names=FALSE, row.names=FALSE)
  # write.table(mcmc_config, paste0("output/mcmc_config.csv"), sep=",", col.names=FALSE, row.names=FALSE)
  
  
  all_config <- rbind(hc_config,mcmc_config,pcstable_config,mmhc_config,FGS_config)
  
  
  
  ################function to convert tetrad format to Matrix format ###############
  TetradGetAdjmat <- function(tetradrunner) {
    p <- length(tetradrunner$nodes)
    adjmat <- matrix(0, p, p)
    for (e in edges) {
      edgevec <- unlist(strsplit(e, " "))
      i <- match(edgevec[1], tetradrunner$nodes)
      j <- match(edgevec[3], tetradrunner$nodes)
      if (edgevec[2] == '---') {
        edge <- c(edgevec[1], edgevec[3])
        adjmat[i, j] <- 3
        adjmat[j, i] <- 3
      } else if (edgevec[2] == "-->") {
        adjmat[i, j] <- 2
        adjmat[j, i] <- 3
      }
      else if (edgevec[2] == "<--") {
        adjmat[i, j] <- 3
        adjmat[j, i] <- 2
      }
      else if (edgevec[2] == "o-o") {
        edge <- c(edgevec[1], edgevec[3])
        adjmat[i, j] <- 1
        adjmat[j, i] <- 1
      }
      else if (edgevec[2] == "o->") {
        adjmat[i, j] <- 2
        adjmat[j, i] <- 1
      }
      else if (edgevec[2] == "<-o") {
        adjmat[i, j] <- 1
        adjmat[j, i] <- 2
      }
      else if (edgevec[2] == "<->") {
        edge <- c(edgevec[1], edgevec[3])
        adjmat[i, j] <- 2
        adjmat[j, i] <- 2
      }
      
    }
    return(adjmat)
  }
  
  ############################################################################
  
  for (set in set_list)
  { 
    
    
    
    
    if (set=='ALARM' | set== 'SPORTS' |set== 'ASIA' |set=='PROPERTY')
    {
      DAG_set <-paste0('Input/DAGtrue_',set,'.csv')
      
      DAG_RAW <- read.csv(DAG_set,header = TRUE,na.strings=c(""),check.name=FALSE)
      
      ########################################################
      
      
      data_matrix_bidirect <-as.matrix(DAG_RAW)
      
      direct_tmp <- matrix(0, nrow = 0, ncol = 2)
      
      
      
      ############## no bi-directed edge to true DAG (if possible)
      for (row_index in 1:nrow(data_matrix_bidirect))
      {
        if (data_matrix_bidirect [row_index,3]!="<->")
        {
          direct <- c(data_matrix_bidirect [row_index,2],data_matrix_bidirect [row_index,4])
          direct <- rbind(direct_tmp , direct)
          direct_tmp <- direct 
        }
      }
      # Convert friends matrix to an igraph object
      g <- graph.edgelist(direct, directed = TRUE)
      NEL <- igraph.to.graphNEL(g)
      true_dag <- as.bn(NEL)
      true_dag_amat <- amat(true_dag)
      true_dag_amat <- true_dag_amat[order(rownames(true_dag_amat)),order(colnames(true_dag_amat))]
      true_dag_amat_temp <- true_dag_amat
      true_dag_amat_temp[true_dag_amat !=0] <-0
      true_dag_amat_temp[true_dag_amat ==0 & t(true_dag_amat)==1] <- 3
      true_dag_amat_temp[true_dag_amat ==1 & t(true_dag_amat)==0] <- 2
      true_pag_amat <- true_dag_amat_temp
    }
    if (set=='hepar2' | set=='hailfinder')
    {
      load(paste0("Input/",set,".rda"))
      
      true_dag_amat <- amat(bn)
      true_dag_amat <- true_dag_amat[order(rownames(true_dag_amat)),order(colnames(true_dag_amat))]
      true_dag_amat_temp <- true_dag_amat
      true_dag_amat_temp[true_dag_amat !=0] <-0
      true_dag_amat_temp[true_dag_amat ==0 & t(true_dag_amat)==1] <- 3
      true_dag_amat_temp[true_dag_amat ==1 & t(true_dag_amat)==0] <- 2
      true_pag_amat <- true_dag_amat_temp
    }
    
    
    for (n in samplesize_list)
    {
      
      log2.txt <-c()
      tic.clearlog()
      tic("SOS runtime")
      
      perf_round <-list()
      avg_perf<- c()
      avg_df<-c()
      if (n==1000)
      {
        training_set = paste0("Input/trainingData_",set,"_N_1k.csv")
      }
      if (n==10000)
      {
        training_set = paste0("Input/trainingData_",set,"_N_10k.csv")
      }
      
      ########## read data sets and convenrt to factors  ######
      dataset_truedag <- read.csv(training_set,header = TRUE,na.strings=c(""),check.names = FALSE)
      dataset_truedag <- dataset_truedag[,order(colnames(dataset_truedag))]
      for(i in 1:ncol(dataset_truedag)){
        
        dataset_truedag[,i] <- as.factor(dataset_truedag[,i])
        
      }
      
      
      ########## check the number of factors for each feather and remove the level with one ######
      
      feather <- sapply(dataset_truedag[,sapply(dataset_truedag, is.factor)], nlevels)
      for (feather_index in colnames(dataset_truedag))
      {
        if (feather[feather_index] ==1)
        {
          dataset_truedag <- dataset_truedag[,!names(dataset_truedag) %in% feather_index]
          true_dag_amat <- true_dag_amat[, colnames(true_dag_amat) != feather_index] 
          true_dag_amat <- true_dag_amat[rownames(true_dag_amat) != feather_index, ] 
        }
      }
      
      ####################################################################################### 
      
      
      for (config_index in 1:nrow(all_config))
      {
        algoname = as.character(all_config[config_index,1])
        testname = as.character(all_config[config_index,2])
        scorename = as.character(all_config[config_index,3])
        alphavalue = (all_config[config_index,4])
        Penalty = (all_config[config_index,5])
        structureprior = (all_config[config_index,6])
        
        if (scorename=='bdeu')  structureprior =  as.numeric(as.character(all_config[config_index,6]))
        if (scorename=='bic')   Penalty =  as.numeric(as.character(all_config[config_index,5]))
        if (algoname=='pc')     alphavalue = as.numeric(as.character((all_config[config_index,4])))
        if (algoname=='mmhc')     alphavalue = as.numeric(as.character((all_config[config_index,4])))
        
        ########## reset list of objective score for each round ##########
        perf <-c()
        df <-c() 
        accuracy <-c()
        
        ########## read learnt DAG from memory  ##########
        list_config_G1 <- readRDS(file=paste0('outputG/',set,'/performance_SOS_dyn_v3.5.1_test_avg_',length(round),'_',set,'_',n,'.rds'))
        list_config_G2 <- readRDS(file=paste0('outputG/',set,'/performance_SOS_dyn_v3.5.1_test(FGS_tetrad)_avg_',length(round),'_',set,'_',n,'.rds'))
        list_config_G<- append(list_config_G1,list_config_G2)
        ##################################################
        
        for (r in round)
        {
          print(paste0('Case =',set,' Sample size =',n,' algorithm=',algoname,' score=',scorename,' test=',testname,' alpha=',alphavalue,' Penalty=',Penalty,' Structure=',structureprior,' round=',r))
          
          train_number = n*(total_fold-1)/total_fold
          set.seed(r)
          dataset <- dataset_truedag[sample(nrow(dataset_truedag),train_number, replace= Rp), ]
          
          #### Do not support latent variable in this project yet ######
          #  dataset <- dataset[ , !(names(dataset) %in% latent)]
          
          dataset_test <- dataset_truedag[!rownames(dataset_truedag) %in%  rownames(dataset), ]
          test_number = n/total_fold
          dataset_test <- dataset_test[sample(nrow(dataset_test),test_number, replace= TRUE), ]
          V = ncol(dataset_test)
          
          # # export train and test data 
          # dataset_train_pcalg <- dataset[,order(colnames(dataset))]
          # dataset_test_pcalg <- dataset_test[,order(colnames(dataset_test))]
          # cols <- colnames(dataset_train_pcalg)
          # for(j in cols){
          #   dataset_train_pcalg[[j]]<- as.numeric(dataset[[j]])
          #   dataset_train_pcalg[[j]]<- dataset_train_pcalg[[j]]-1
          # }
          # cols <- colnames(dataset_test_pcalg)
          # for(j in cols){
          #   dataset_test_pcalg[[j]]<- as.numeric(dataset_test[[j]])
          #   dataset_test_pcalg[[j]]<- dataset_test_pcalg[[j]]-1
          # }
          #write.table(dataset_train_pcalg, paste0("output/",set,"/train_",r,"_",set,"_",n,".csv"), sep=",", col.names=FALSE, row.names=FALSE)
          #write.table(dataset_test_pcalg, paste0("output/",set,"/test_",r,"_",set,"_",n,".csv"), sep=",", col.names=FALSE, row.names=FALSE)
          
          
          if (algoname =='FGS_tetrad')
          {
            tetradrunner <- list_config_G[[config_index]][[r]]
            edges <- tetradrunner$edges
            fci.est <-  TetradGetAdjmat(tetradrunner)
            
            # Convert pag matrix to bnlearn
            learnt_cpdag_amat <- fci.est 
            learnt_cpdag_amat[fci.est!=0]<-0
            learnt_cpdag_amat[fci.est==2 & t(fci.est)==3]<-1
            learnt_cpdag_amat[fci.est==3 & t(fci.est)==2]<-0
            learnt_cpdag_amat[fci.est==3 & t(fci.est)==3]<-1
            
            learnt_cpdag = empty.graph(colnames(dataset))
            
            amat(learnt_cpdag)= learnt_cpdag_amat
            sample_dag_from_cpdag <- cextend(learnt_cpdag, strict = TRUE, debug = FALSE)
            
            
            
            {
              performance_forest_list <- list()
              
              free_param <- nparams(sample_dag_from_cpdag,dataset_test)
              
              if (scorename=='bdeu')
              {
                performance_bde <- score(sample_dag_from_cpdag, dataset_test, type = "bde", iss= structureprior)
                performance_bic <- score(sample_dag_from_cpdag, dataset_test, type = "bic",  k= 1) -  free_param* log(V)
                
                #performance_bic <- score(hc_runner, dataset_test, type = "bic",  k= 1)
                #performance_bic <- score(hc_runner,dataset_test,type = "custom",fun=ebic,args=list(gamma_penalty=0))
                performance_aic <- score(sample_dag_from_cpdag, dataset_test, type = "aic",  k= 1)
                #performance_bic <- score(hc_runner, dataset_test, type = "bic",  k= 1) -  (1-(1/20)*(structureprior-1)) *free_param* log(V)
                
              }
              
              if (scorename=='bic')
              {
                performance_bde <- score(sample_dag_from_cpdag, dataset_test, type = "bde", iss= 1)
                
                performance_bic <- score(sample_dag_from_cpdag, dataset_test, type = "bic",  k= 1) -  (1/20)*(Penalty-1) *free_param* log(V)
                #performance_bic <- score(hc_runner,dataset_test,type = "custom",fun=ebic,args=list(gamma_penalty=(1/20)*(Penalty-1)))
                performance_aic <- score(sample_dag_from_cpdag, dataset_test, type = "aic",  k= 1)
              }
              
              performance_ll <- score(sample_dag_from_cpdag, dataset_test, type = "loglik")
              performance_mi <- average_model_MI(sample_dag_from_cpdag, dataset_test)$total_MI
              performance_mi_sh <- average_model_MI(sample_dag_from_cpdag, dataset_test)$total_MI_sh
              performance_avg_mi <- average_model_MI(sample_dag_from_cpdag, dataset_test)$avg_MI
              performance_avg_mi_sh <- average_model_MI(sample_dag_from_cpdag, dataset_test)$avg_MI_sh
              
              performance <- c(performance_bde,performance_bic,performance_ll,performance_aic,free_param, performance_avg_mi,performance_avg_mi_sh,performance_mi,performance_mi_sh)
              
              
            }
            
            
            
            perf <- rbind(perf, performance) 
            value = data.frame(perf)
            colnames(value) <- c("Avg_BDeu","Avg_bic","Avg_ll","Avg_aic","avg_free_param","Avg_mi","Avg_mi_sh","total_mi","total_mi_sh")
            
            
          }
          
          if (algoname =='hc')
          {
            
            
            hc_runner <- list_config_G[[config_index]][[r]]
            
            edges <- hc_runner$arcs
            data_matrix <-amat(hc_runner)
            
            
            learnt_dag_amat <- true_pag_amat 
            learnt_dag_amat[learnt_dag_amat!=0]<-0
            learnt_dag_amat[data_matrix==1 & t(data_matrix)==0]<-2
            learnt_dag_amat[data_matrix==0 & t(data_matrix)==1]<-3
            
            # # export train and test data 
            # dataset_train_pcalg <- dataset[,order(colnames(dataset))]
            # dataset_test_pcalg <- dataset_test[,order(colnames(dataset_test))]
            # cols <- colnames(dataset_train_pcalg)
            # for(j in cols){
            #   dataset_train_pcalg[[j]]<- as.numeric(dataset[[j]])
            #   dataset_train_pcalg[[j]]<- dataset_train_pcalg[[j]]-1
            # }
            # cols <- colnames(dataset_test_pcalg)
            # for(j in cols){
            #   dataset_test_pcalg[[j]]<- as.numeric(dataset_test[[j]])
            #   dataset_test_pcalg[[j]]<- dataset_test_pcalg[[j]]-1
            # }
            # 
            # if (scorename=='bdeu')
            # {
            #   write.table(learnt_dag_amat, paste0("output/",set,"/dag_",algoname,"_",scorename,"_",structureprior,"_",r,"_",set,"_",n,".csv"), sep=",", col.names=FALSE, row.names=FALSE)
            #   
            # }
            # if (scorename=='bic')
            # {
            #   write.table(learnt_dag_amat, paste0("output/",set,"/dag_",algoname,"_",scorename,"_",Penalty,"_",r,"_",set,"_",n,".csv"), sep=",", col.names=FALSE, row.names=FALSE)
            #   
            # }
            # 
            
            
            {
              
              performance_forest_list <- list()
              
              
              free_param <- nparams(hc_runner,dataset_test)
              if (scorename=='bdeu')
              {
                performance_bde <- score(hc_runner, dataset_test, type = "bde", iss= structureprior)
                performance_bic <- score(hc_runner, dataset_test, type = "bic",  k= 1) -  free_param* log(V)
                
                #performance_bic <- score(hc_runner, dataset_test, type = "bic",  k= 1)
                #performance_bic <- score(hc_runner,dataset_test,type = "custom",fun=ebic,args=list(gamma_penalty=0))
                performance_aic <- score(hc_runner, dataset_test, type = "aic",  k= 1)
                #performance_bic <- score(hc_runner, dataset_test, type = "bic",  k= 1) -  (1-(1/20)*(structureprior-1)) *free_param* log(V)
                
              }
              
              if (scorename=='bic')
              {
                performance_bde <- score(hc_runner, dataset_test, type = "bde", iss= 1)
                
                
                performance_bic <- score(hc_runner, dataset_test, type = "bic",  k= 1) -  (1/20)*(Penalty-1) *free_param* log(V)
                #performance_bic <- score(hc_runner,dataset_test,type = "custom",fun=ebic,args=list(gamma_penalty=(1/20)*(Penalty-1)))
                performance_aic <- score(hc_runner, dataset_test, type = "aic",  k= 1)
              }
              
              
              performance_ll <- score(hc_runner, dataset_test,type = "loglik")
              performance_mi <- average_model_MI(hc_runner, dataset_test)$total_MI
              performance_mi_sh <- average_model_MI(hc_runner, dataset_test)$total_MI_sh
              performance_avg_mi <- average_model_MI(hc_runner, dataset_test)$avg_MI
              performance_avg_mi_sh <- average_model_MI(hc_runner, dataset_test)$avg_MI_sh
              
              
              performance <- c(performance_bde,performance_bic,performance_ll,performance_aic,free_param, performance_avg_mi,performance_avg_mi_sh,performance_mi,performance_mi_sh)
              
              
            }
            
            
            perf <- rbind(perf, performance) ######## average KL of hyperparameter
            value = data.frame(perf)
            colnames(value) <- c("Avg_BDeu","Avg_bic","Avg_ll","Avg_aic","avg_free_param","Avg_mi","Avg_mi_sh","total_mi","total_mi_sh") 
            
          }
          
          if (algoname =='pc')
          {
            
            pc_runner <- list_config_G[[config_index]][[r]] 
            ERROR <- FALSE
            tryCatch( { sample_dag_from_cpdag <- cextend(pc_runner, strict = TRUE, debug = FALSE); ERROR <<- FALSE }
                      , error = function(e) {ERROR <<- TRUE})
            
            
            if(!ERROR)
            {
              
              data_matrix <-amat(pc_runner)
              
              learnt_dag_amat <- true_pag_amat 
              learnt_dag_amat[learnt_dag_amat!=0]<-0
              learnt_dag_amat[data_matrix==1 & t(data_matrix)==0]<-2
              learnt_dag_amat[data_matrix==0 & t(data_matrix)==1]<-3
              learnt_dag_amat[data_matrix==1 & t(data_matrix)==1]<-3
              
              # # export train and test data 
              # dataset_train_pcalg <- dataset[,order(colnames(dataset))]
              # dataset_test_pcalg <- dataset_test[,order(colnames(dataset_test))]
              # cols <- colnames(dataset_train_pcalg)
              # for(j in cols){
              #   dataset_train_pcalg[[j]]<- as.numeric(dataset[[j]])
              #   dataset_train_pcalg[[j]]<- dataset_train_pcalg[[j]]-1
              # }
              # cols <- colnames(dataset_test_pcalg)
              # for(j in cols){
              #   dataset_test_pcalg[[j]]<- as.numeric(dataset_test[[j]])
              #   dataset_test_pcalg[[j]]<- dataset_test_pcalg[[j]]-1
              # }
              
              {
                
                performance_forest_list <- list()
                
                
                free_param <- nparams(sample_dag_from_cpdag,dataset_test)
                performance_bic <- score(sample_dag_from_cpdag, dataset_test, type = "bic",  k= 1)
                #performance_bic <- score(sample_dag_from_cpdag,dataset_test, type = "custom",fun=ebic,args=list(gamma_penalty=0))
                performance_bde <- score(sample_dag_from_cpdag, dataset_test, type = "bde", iss= 1)
                
                #performance_bic <- score(sample_dag_from_cpdag,dataset,type = "custom",fun=ebic,args=list(gamma_penalty=0))
                
                performance_aic <- score(sample_dag_from_cpdag, dataset_test, type = "aic",  k= 1)
                performance_ll <- score(sample_dag_from_cpdag, dataset_test, type = "loglik")
                performance_mi <- average_model_MI(sample_dag_from_cpdag, dataset_test)$total_MI
                performance_mi_sh <- average_model_MI(sample_dag_from_cpdag, dataset_test)$total_MI_sh
                performance_avg_mi <- average_model_MI(sample_dag_from_cpdag, dataset_test)$avg_MI
                performance_avg_mi_sh <- average_model_MI(sample_dag_from_cpdag, dataset_test)$avg_MI_sh
                
                performance <- c(performance_bde,performance_bic,performance_ll,performance_aic,free_param, performance_avg_mi,performance_avg_mi_sh,performance_mi,performance_mi_sh)
                
              }
              
            }
            else
            {
              performance <- c(NA,NA,NA,NA,NA, NA,NA,NA,NA)
              
            }
            
            perf <- rbind(perf, performance) ######## average KL of hyperparameter
            value = data.frame(perf)
            colnames(value) <- c("Avg_BDeu","Avg_bic","Avg_ll","Avg_aic","avg_free_param","Avg_mi","Avg_mi_sh","total_mi","total_mi_sh")
            perf_round[[r]] <- perf
            
            
          }
          if (algoname =='mcmc')
          {
            
            edge_amat <- list_config_G[[config_index]][[r]] 
            
            
            learnt_dag_amat <- true_pag_amat 
            learnt_dag_amat[learnt_dag_amat!=0]<-0
            learnt_dag_amat[edge_amat==1 & t(edge_amat)==0]<-2
            learnt_dag_amat[edge_amat==0 & t(edge_amat)==1]<-3
            
            
            learnt_dag = empty.graph(colnames(edge_amat))
            
            amat(learnt_dag)= edge_amat
            
            # # export train and test data 
            # dataset_train_pcalg <- dataset[,order(colnames(dataset))]
            # dataset_test_pcalg <- dataset_test[,order(colnames(dataset_test))]
            # cols <- colnames(dataset_train_pcalg)
            # for(j in cols){
            #   dataset_train_pcalg[[j]]<- as.numeric(dataset[[j]])
            #   dataset_train_pcalg[[j]]<- dataset_train_pcalg[[j]]-1
            # }
            # cols <- colnames(dataset_test_pcalg)
            # for(j in cols){
            #   dataset_test_pcalg[[j]]<- as.numeric(dataset_test[[j]])
            #   dataset_test_pcalg[[j]]<- dataset_test_pcalg[[j]]-1
            # }
            # 
            # 
            #write.table(learnt_dag_amat, paste0("output/",set,"/dag_",algoname,"_",scorename,"_",structureprior,"_",r,"_",set,"_",n,".csv"), sep=",", col.names=FALSE, row.names=FALSE)
            
            
            {
              
              performance_forest_list <- list()
              
              free_param <- nparams(learnt_dag,dataset_test)
              if (scorename=='bdeu')
              {
                performance_bde <- score(learnt_dag, dataset_test, type = "bde", iss= structureprior)
                
                #performance_bic <- score(learnt_dag, dataset_test, type = "bic",  k= 1) - (1-(1/10)*(structureprior) *free_param* log(V))
                
                performance_bic <- score(learnt_dag, dataset_test, type = "bic",  k= 1)
                #performance_bic <- score(learnt_dag,dataset_test,type = "custom",fun=ebic,args=list(gamma_penalty=0))
                performance_aic <- score(learnt_dag, dataset_test, type = "aic",  k= 1)
              }
              
              
              performance_ll <- score(learnt_dag, dataset_test, type = "loglik")
              performance_mi <- average_model_MI(learnt_dag, dataset_test)$total_MI
              performance_mi_sh <- average_model_MI(learnt_dag, dataset_test)$total_MI_sh
              performance_avg_mi <- average_model_MI(learnt_dag, dataset_test)$avg_MI
              performance_avg_mi_sh <- average_model_MI(learnt_dag, dataset_test)$avg_MI_sh
              
              performance <- c(performance_bde,performance_bic,performance_ll,performance_aic,free_param, performance_avg_mi,performance_avg_mi_sh,performance_mi,performance_mi_sh)
              
              
              
            }
            
            
            perf <- rbind(perf, performance)
            value = data.frame(perf)
            colnames(value) <- c("Avg_BDeu","Avg_bic","Avg_ll","Avg_aic","avg_free_param","Avg_mi","Avg_mi_sh","total_mi","total_mi_sh")
            
            
          }
          if (algoname =='mmhc')
          {
            
            mmhc_runner <- list_config_G[[config_index]][[r]] 
            
            edges <- mmhc_runner$arcs
            data_matrix <-amat(mmhc_runner)
            
            learnt_dag_amat <- true_pag_amat 
            learnt_dag_amat[learnt_dag_amat!=0]<-0
            learnt_dag_amat[data_matrix==1 & t(data_matrix)==0]<-2
            learnt_dag_amat[data_matrix==0 & t(data_matrix)==1]<-3
            
            # # export train and test data 
            # dataset_train_pcalg <- dataset[,order(colnames(dataset))]
            # dataset_test_pcalg <- dataset_test[,order(colnames(dataset_test))]
            # cols <- colnames(dataset_train_pcalg)
            # for(j in cols){
            #   dataset_train_pcalg[[j]]<- as.numeric(dataset[[j]])
            #   dataset_train_pcalg[[j]]<- dataset_train_pcalg[[j]]-1
            # }
            # cols <- colnames(dataset_test_pcalg)
            # for(j in cols){
            #   dataset_test_pcalg[[j]]<- as.numeric(dataset_test[[j]])
            #   dataset_test_pcalg[[j]]<- dataset_test_pcalg[[j]]-1
            # }
            # 
            
            
            
            {
              
              performance_forest_list <- list()
              
              
              
              free_param <- nparams(mmhc_runner,dataset_test)
              
              
              if (scorename=='bdeu')
              {
                performance_bde <- score(mmhc_runner, dataset_test, type = "bde", iss= structureprior)
                performance_bic <- score(mmhc_runner, dataset_test, type = "bic",  k= 1) -  free_param* log(V)
                
                #performance_bic <- score(hc_runner, dataset_test, type = "bic",  k= 1)
                #performance_bic <- score(hc_runner,dataset_test,type = "custom",fun=ebic,args=list(gamma_penalty=0))
                performance_aic <- score(mmhc_runner, dataset_test, type = "aic",  k= 1)
                #performance_bic <- score(hc_runner, dataset_test, type = "bic",  k= 1) -  (1-(1/20)*(structureprior-1)) *free_param* log(V)
                
              }
              
              if (scorename=='bic')
              {
                performance_bde <- score(mmhc_runner, dataset_test, type = "bde", iss= 1)
                
                
                performance_bic <- score(mmhc_runner, dataset_test, type = "bic",  k= 1) -  (1/10)*(Penalty-1) *free_param* log(V)
                #performance_bic <- score(hc_runner,dataset_test,type = "custom",fun=ebic,args=list(gamma_penalty=(1/20)*(Penalty-1)))
                performance_aic <- score(mmhc_runner, dataset_test, type = "aic",  k= 1)
              }
              
              
              performance_ll <- score(mmhc_runner, dataset_test,  type = "loglik")
              performance_mi <- average_model_MI(mmhc_runner, dataset_test)$total_MI
              performance_mi_sh <- average_model_MI(mmhc_runner, dataset_test)$total_MI_sh
              performance_avg_mi <- average_model_MI(mmhc_runner, dataset_test)$avg_MI
              performance_avg_mi_sh <- average_model_MI(mmhc_runner, dataset_test)$avg_MI_sh
              performance <- c(performance_bde,performance_bic,performance_ll,performance_aic,free_param, performance_avg_mi,performance_avg_mi_sh,performance_mi,performance_mi_sh)
              
              
            }
            
            
            
            perf <- rbind(perf, performance) 
            value = data.frame(perf)
            colnames(value) <- c("Avg_BDeu","Avg_bic","Avg_ll","Avg_aic","avg_free_param","Avg_mi","Avg_mi_sh","total_mi","total_mi_sh") 
            
            
          }
          
          
        }
        
        ################### graphical accuracy of observational data ############################
        
        {
          
          
          if (algoname =='FGS_tetrad')
          {
            
            
            tetradrunner <- readRDS(file=paste0('outputG/',set,'/performance_SOS_dyn_v3.5.1_dataset_truedag_(',algoname,')_avg_',config_index-(nrow(all_config)-40),'_',set,'_',n,'.rds'))
            
            
            edges <- tetradrunner$edges
            fci.est <-  TetradGetAdjmat(tetradrunner)
            
            # Convert pag matrix to bnlearn
            learnt_cpdag_amat <- fci.est 
            learnt_cpdag_amat[fci.est!=0]<-0
            learnt_cpdag_amat[fci.est==2 & t(fci.est)==3]<-1
            learnt_cpdag_amat[fci.est==3 & t(fci.est)==2]<-0
            learnt_cpdag_amat[fci.est==3 & t(fci.est)==3]<-1
            
            learnt_cpdag = empty.graph(colnames(dataset))
            
            amat(learnt_cpdag)= learnt_cpdag_amat
            sample_dag_from_cpdag <- cextend(learnt_cpdag, strict = TRUE, debug = FALSE)
            
            
            results <- precisionrecall(fci.est,true_dag_amat_temp)
            
            
            df_tempt<- data.frame(matrix(unlist(results), ncol = length(results), byrow=TRUE))
            
            dfFitList <- lapply(results, '[', 'fitted')
            
            colnames(df_tempt) <- names(dfFitList)
            df <- df_tempt
            colnames(df) <- names(dfFitList)
            
          }
          
          if (algoname =='hc')
          {
            
            
            hc_runner <- readRDS(file=paste0('outputG/',set,'/performance_SOS_dyn_v3.5.1_dataset_truedag_(',algoname,')_avg_',config_index,'_',set,'_',n,'.rds'))
            
            edges <- hc_runner$arcs
            data_matrix <-amat(hc_runner)
            
            learnt_dag_amat <- true_pag_amat 
            learnt_dag_amat[learnt_dag_amat!=0]<-0
            learnt_dag_amat[data_matrix==1 & t(data_matrix)==0]<-2
            learnt_dag_amat[data_matrix==0 & t(data_matrix)==1]<-3
            
            results <- precisionrecall(learnt_dag_amat,true_dag_amat_temp)
            df_tempt<- data.frame(matrix(unlist(results), ncol = length(results), byrow=TRUE))
            
            dfFitList <- lapply(results, '[', 'fitted')
            
            df <- df_tempt
            
          }
          
          if (algoname =='pc')
          {
            
            pc_runner <- readRDS(file=paste0('outputG/',set,'/performance_SOS_dyn_v3.5.1_dataset_truedag_(',algoname,')_avg_',config_index,'_',set,'_',n,'.rds'))
            
            ERROR <- FALSE
            tryCatch( { sample_dag_from_cpdag <- cextend(pc_runner, strict = TRUE, debug = FALSE); ERROR <<- FALSE }
                      , error = function(e) {ERROR <<- TRUE})
            
            if(!ERROR)
            {
              
              data_matrix <-amat(pc_runner)
              
              learnt_dag_amat <- true_pag_amat 
              learnt_dag_amat[learnt_dag_amat!=0]<-0
              learnt_dag_amat[data_matrix==1 & t(data_matrix)==0]<-2
              learnt_dag_amat[data_matrix==0 & t(data_matrix)==1]<-3
              learnt_dag_amat[data_matrix==1 & t(data_matrix)==1]<-3
              
              results <- precisionrecall(learnt_dag_amat,true_dag_amat_temp)
              df_tempt<- data.frame(matrix(unlist(results), ncol = length(results), byrow=TRUE))
              
              dfFitList <- lapply(results, '[', 'fitted')
              
              df <- df_tempt
              
            }
            else
            {
              df <-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
              
            }
            
          }
          if (algoname =='mcmc')
          {
            
            edge_amat <- readRDS(file=paste0('outputG/',set,'/performance_SOS_dyn_v3.5.1_dataset_truedag_(',algoname,')_avg_',config_index,'_',set,'_',n,'.rds'))
            
            learnt_dag_amat <- true_pag_amat 
            learnt_dag_amat[learnt_dag_amat!=0]<-0
            learnt_dag_amat[edge_amat==1 & t(edge_amat)==0]<-2
            learnt_dag_amat[edge_amat==0 & t(edge_amat)==1]<-3
            
            learnt_dag = empty.graph(colnames(edge_amat))
            
            amat(learnt_dag)= edge_amat
            
            results <- precisionrecall(learnt_dag_amat,true_dag_amat_temp)
            df_tempt<- data.frame(matrix(unlist(results), ncol = length(results), byrow=TRUE))
            
            
            dfFitList <- lapply(results, '[', 'fitted')
            
            
            df <- df_tempt
            
          }
          if (algoname =='mmhc')
          {
            
            mmhc_runner <- readRDS(file=paste0('outputG/',set,'/performance_SOS_dyn_v3.5.1_dataset_truedag_(',algoname,')_avg_',config_index,'_',set,'_',n,'.rds'))
            
            edges <- mmhc_runner$arcs
            data_matrix <-amat(mmhc_runner)
            
            learnt_dag_amat <- true_pag_amat 
            learnt_dag_amat[learnt_dag_amat!=0]<-0
            learnt_dag_amat[data_matrix==1 & t(data_matrix)==0]<-2
            learnt_dag_amat[data_matrix==0 & t(data_matrix)==1]<-3
            
            results <- precisionrecall(learnt_dag_amat,true_dag_amat_temp)
            df_tempt<- data.frame(matrix(unlist(results), ncol = length(results), byrow=TRUE))
            
            dfFitList <- lapply(results, '[', 'fitted')
            
            df <- df_tempt
            
          }
          
        }
        
        avg_perf_temp <-colMeans(value,na.rm = TRUE)
        
        df_temp<- c(config_index,df)
        df_temp = data.frame(df_temp)
        colnames(df_temp) <- c("config id",	"precision",	"recall",	"BSF",	"f1","shd",	"tp",	"fp",	"fn",	"tn",	"nEdgesGtPag","a","i","nEdgesPag","tp_undirect","tp_bidirect","tp_o_direct")
        
        accuracy<- rbind(accuracy,df_temp)
        
        avg_perf_tempt<- c(accuracy,avg_perf_temp)
        avg_perf_tempt = data.frame(avg_perf_tempt)
        avg_perf<- rbind(avg_perf,avg_perf_tempt)
        
      }
      avg_perf <- cbind(all_config,avg_perf)
      
      ##################### select algorithms with dynamic score AIC/MIC ##########################
      {
        
        ##################### select algorithms and the hyper-paramaters with dynamic score BDeu/BIC/MIC ##########################
        avg_perf$tuning_F1 <- NA
        avg_perf$tuning_SHD <-  NA
        if (tuning_score =='BDeu') tuning_score_col = 24
        if (tuning_score =='BIC') tuning_score_col = 25
        if (tuning_score =='AIC') tuning_score_col = 27
        if (tuning_score =='LL') tuning_score_col = 26
        if (tuning_score =='MI') tuning_score_col = 31
        if (tuning_score =='MI-sh') tuning_score_col = 32
        if (tuning_score =='avg_MI') tuning_score_col = 29
        if (tuning_score =='avg_MI-sh') tuning_score_col = 30
        
        for (config_index3 in 1:nrow(all_config))
        {
          
          if (all_config[config_index3,1] == 'mmhc')
          {
            if (tuning_score=='BIC')
            {
              default_position = 21
            }
            if (tuning_score=='BDeu')
            {
              default_position = 1
            }
            
            range_compare=40
            max_range = config_index3+range_compare-1
            range = config_index3:max_range
            max_range2 = max_range+(range_compare)
            range2 = (max_range+1):max_range2
            
            
            optimal_config = which.max(avg_perf[range,tuning_score_col])
            #F1
            avg_perf[range[optimal_config],33] =   (as.numeric(avg_perf[range[optimal_config],11])  - as.numeric(avg_perf[range[default_position],11]))/as.numeric(avg_perf[range[default_position],11])
            
            #SHD
            avg_perf[range[optimal_config],34] =   (as.numeric(avg_perf[range[default_position],12])- as.numeric(avg_perf[range[optimal_config],12]))/as.numeric(avg_perf[range[default_position],12])
            write.csv(optimal_config, paste0('output/',set,'/bestconfig_SOS_v4_',tuning_score,"_",all_config[config_index3,1],'_',set,'_',n,'.csv'), row.names=F)
            
            
            break
            
            
          }
          
        } 
        for (config_index3 in 1:nrow(all_config))
        {
          if (all_config[config_index3,1] == 'hc')
          {
            if (tuning_score=='BIC')
            {
              default_position = 21
            }
            if (tuning_score=='BDeu')
            {
              default_position = 1
            }
            range_compare=40
            max_range = config_index3+range_compare-1
            range = config_index3:max_range
            max_range2 = max_range+(range_compare)
            range2 = (max_range+1):max_range2
            
            optimal_config = which.max(avg_perf[range,tuning_score_col])
            #F1
            avg_perf[range[optimal_config],33] =   (as.numeric(avg_perf[range[optimal_config],11])  - as.numeric(avg_perf[range[default_position],11]))/as.numeric(avg_perf[range[default_position],11])
            
            #SHD
            avg_perf[range[optimal_config],34] =   (as.numeric(avg_perf[range[default_position],12])- as.numeric(avg_perf[range[optimal_config],12]))/as.numeric(avg_perf[range[default_position],12])
            write.csv(optimal_config, paste0('output/',set,'/bestconfig_SOS_v4_',tuning_score,"_",all_config[config_index3,1],'_',set,'_',n,'.csv'), row.names=F)
            break
            
            
          }
          
        } 
        
        for (config_index3 in 1:nrow(all_config))
        {
          if (tuning_score=='BIC')
          {
            default_position = 21
          }
          if (tuning_score=='BDeu')
          {
            default_position = 1
          }
          if (all_config[config_index3,1] == 'FGS_tetrad')
          {
            range_compare=40
            max_range = config_index3+range_compare-1
            range = config_index3:max_range
            max_range2 = max_range+(range_compare)
            range2 = (max_range+1):max_range2
            
            optimal_config = which.max(avg_perf[range,tuning_score_col])
            #F1
            avg_perf[range[optimal_config],33] =   (as.numeric(avg_perf[range[optimal_config],11])  - as.numeric(avg_perf[range[default_position],11]))/as.numeric(avg_perf[range[default_position],11])
            
            #SHD
            avg_perf[range[optimal_config],34] =   (as.numeric(avg_perf[range[default_position],12])- as.numeric(avg_perf[range[optimal_config],12]))/as.numeric(avg_perf[range[default_position],12])
            write.csv(optimal_config, paste0('output/',set,'/bestconfig_SOS_v4_',tuning_score,"_",all_config[config_index3,1],'_',set,'_',n,'.csv'), row.names=F)
            
            
            break
          }
          
        } 
        
        for (config_index3 in 1:nrow(all_config))
        {
          
          if (all_config[config_index3,1] == 'mcmc')
          {
            range_compare=20
            max_range = config_index3+range_compare-1
            range = config_index3:max_range
            
            optimal_config = which.max(avg_perf[range,tuning_score_col])
            #F1
            avg_perf[range[optimal_config],33] =   (as.numeric(avg_perf[range[optimal_config],11])  - as.numeric(avg_perf[range[1],11]))/as.numeric(avg_perf[range[1],11])
            
            #SHD
            avg_perf[range[optimal_config],34] =   (as.numeric(avg_perf[range[1],12])- as.numeric(avg_perf[range[optimal_config],12]))/as.numeric(avg_perf[range[1],12])
            write.csv(optimal_config, paste0('output/',set,'/bestconfig_SOS_v4_',tuning_score,"_",all_config[config_index3,1],'_',set,'_',n,'.csv'), row.names=F)
            
            break
          }
        } 
        
        for (config_index3 in 1:nrow(all_config))
        {
          
          if (all_config[config_index3,1] == 'pc')
          {
            range_compare=9
            max_range = config_index3+range_compare-1
            range = config_index3:max_range
            max_range2 = max_range+(range_compare)
            range2 = (max_range+1):max_range2
            optimal_config = which.max(avg_perf[range,tuning_score_col])
            #F1
            avg_perf[range[optimal_config],33] =   (as.numeric(avg_perf[range[optimal_config],11])  - as.numeric(avg_perf[range[1],11]))/as.numeric(avg_perf[range[1],11])
            
            #SHD
            avg_perf[range[optimal_config],34] =   (as.numeric(avg_perf[range[1],12])- as.numeric(avg_perf[range[optimal_config],12]))/as.numeric(avg_perf[range[1],12])
            write.csv(optimal_config, paste0('output/',set,'/bestconfig_SOS_v4_',tuning_score,"_",all_config[config_index3,1],'_',set,'_',n,'.csv'), row.names=F)
            
            break
          }
        } 
        
        
      }
      
      #write.csv(avg_df, paste0("output/",set,"/graphicalaccuracy_SOS_dyn_",algoname,"_",scorename,"_avg",length(round),"_",set,"_",n,".csv"), row.names=F
      write.csv(avg_perf, paste0("output/",set,"/performance_SOS_dyn_v4_test_(",tuning_score,")_avg_",length(round),"_",set,"_",n,".csv"), row.names=F)
      
      
      
      ############## STOP #########################
      toc(log = TRUE, quiet = TRUE)
      log.txt <- tic.log(format = TRUE)
      log2.txt <- rbind(log2.txt, log.txt)
      #############################################
      
      write.csv(log2.txt, paste0("output/",set,"/runtime_SOS_dyn_v4_(",tuning_score,")_",set,"_",n,".csv"), row.names=F)
      
    }
  }
}