# SOS_dynamic v2.1.1: to average all metrics from BIC (no penalty) and MI/nE with input algos and their hyperparams 
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
source("exportamattoCSV.R")
############### change configurations ################
Rp = TRUE # bootstraping 
round = seq(1,10, by=1)
total_fold = 10
Pd =seq(1,20, by=1)
Prior =seq(1,20, by=1)
alpha_list = c(0.05, 0.01, 0.1)
latent <- c("")
row_graph=1
set_list = c('diarrhoea')
#set_list = c('ASIA')
algoname = 'hc'
none_value = rep(c('none'), times=20)
alphavalue = 0.05

##### input parameter ####
structureprior =1
Penalty =1
scorename='bdeu'
testname='x2'


###############################################################################################################

for (set in set_list)
{ 
  
  {
    
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
    
  }
  
  
  
  {
    
    log2.txt <-c()
    
    
    perf_round <-list()
    avg_perf<- c()
    avg_df<-c()
    
    training_set = paste0("Input/trainingData_",set,".csv")
    
    
    dataset_truedag <- read.csv(training_set,header = TRUE,na.strings=c(""),check.names = FALSE)
   # dataset_truedag <- dataset_truedag[sample(nrow(dataset_truedag),train_number, replace= Rp), ]
    
    
    dataset_truedag <- dataset_truedag[,order(colnames(dataset_truedag))]
    n= nrow(dataset_truedag) 
    for(i in 1:ncol(dataset_truedag)){
      
      dataset_truedag[,i] <- as.factor(dataset_truedag[,i])
      
    }
    
    
    ########## check the number of factors for each feather ######
    
    feather <- sapply(dataset_truedag[,sapply(dataset_truedag, is.factor)], nlevels)
    for (feather_index in colnames(dataset_truedag))
    {
      if (feather[feather_index] ==1)
      {
        dataset_truedag <- dataset_truedag[,!names(dataset_truedag) %in% feather_index]
       
      }
    }
    
    ############################################### 
    if (algoname =='mcmc')
    {
      
      tic.clearlog()
      tic("structure runtime")
      
      
      myScore<-scoreparameters("bdecat",dataset_truedag,bdecatpar = list(chi = 0.05, edgepf = structureprior))
      MCMCchains<-list()
      set.seed(1)
      MCMCchains[[1]]<-sampleBN(myScore,"order", hardlimit = 2000)
      edge_amat<-as.matrix(getDAG(MCMCchains[[1]],amat = TRUE, cp = FALSE))
      
      
      learnt_dag_amat <- matrix(data =0,nrow=ncol(dataset_truedag),ncol=ncol(dataset_truedag)) 
      learnt_dag_amat[learnt_dag_amat!=0]<-0
      learnt_dag_amat[edge_amat==1 & t(edge_amat)==0]<-2
      learnt_dag_amat[edge_amat==0 & t(edge_amat)==1]<-3
      
      
      learnt_dag = empty.graph(colnames(edge_amat))
      
      amat(learnt_dag)= edge_amat
      #graphviz.plot(learnt_dag, layout = "dot")
      
      
      edges <- learnt_dag$arcs
      
      
      data_matrix <-amat(learnt_dag)
      
      varNames = colnames(dataset_truedag)
      learnt_PAG<- adjtodot(learnt_dag_amat, varNames)
      dot(learnt_PAG)
      
      ############## STOP #########################
      toc(log = TRUE, quiet = TRUE)
      log.txt <- tic.log(format = TRUE)
      log2.txt <- rbind(log2.txt, log.txt)
      #############################################
       
    }

    if (algoname =='hc')
    {
      
      tic.clearlog()
      tic("structure runtime")
      
      
      if (scorename=='bdeu')  hc_runner <- hc(dataset_truedag,score = 'bde',iss= structureprior)
      if (scorename=='bic') hc_runner <- hc(dataset_truedag,score = 'bic',k= Penalty)
     # hc_runner <- h2pc(dataset_truedag)
      
      edges <- hc_runner$arcs
      data_matrix <-amat(hc_runner)
      
      
      
      learnt_dag_amat <-  matrix(0, length(dataset_truedag), length(dataset_truedag))
      learnt_dag_amat[learnt_dag_amat!=0]<-0
      learnt_dag_amat[data_matrix==1 & t(data_matrix)==0]<-2
      learnt_dag_amat[data_matrix==0 & t(data_matrix)==1]<-3

      varNames = colnames(dataset_truedag)
      learnt_PAG<- adjtodot(learnt_dag_amat, varNames)
      dot(learnt_PAG)
      
      ############## STOP #########################
      toc(log = TRUE, quiet = TRUE)
      log.txt <- tic.log(format = TRUE)
      log2.txt <- rbind(log2.txt, log.txt)
      #############################################
      
      
      write.csv(log2.txt, paste0("output/",set,"/runtime_SOS_dyn_",algoname,"_",scorename,"_",set,".csv"), row.names=F)
      
    }
    
    if (algoname =='mmhc')
    {
      
      
      tic.clearlog()
      tic("structure runtime")
      
      if (scorename=='bdeu')mmhc_runner <- mmhc(dataset_truedag,  restrict.args=list(test = testname, alpha = alphavalue),maximize.args = list(score = 'bde',iss= structureprior))
      
      if (scorename=='bic') mmhc_runner <- mmhc(dataset_truedag,  restrict.args=list(test = testname, alpha = alphavalue),maximize.args = list(score = 'bic',k=Penalty ))
     
      
      #gs_runner <- gs(dataset_truedag,  test = testname, alpha = alphavalue, undirected = FALSE)
      
      #mmhc_runner <- gs_runner 
      edges <- mmhc_runner$arcs
      data_matrix <-amat(mmhc_runner)
      
      
      # 
      learnt_dag_amat <- matrix(0,nrow=ncol(dataset_truedag),ncol=ncol(dataset_truedag))
      learnt_dag_amat[data_matrix==1 & t(data_matrix)==0]<-2
      learnt_dag_amat[data_matrix==0 & t(data_matrix)==1]<-3
     
     
      varNames = colnames(dataset_truedag)
      colnames(learnt_dag_amat) =colnames(dataset_truedag)
      rownames(learnt_dag_amat) =colnames(dataset_truedag)
      learnt_PAG<- adjtodot(learnt_dag_amat, varNames)
      dot(learnt_PAG)
      
      write.csv(log2.txt, paste0("output/",set,"/runtime_SOS_dyn_",algoname,"_",scorename,"_",set,".csv"), row.names=F)
      
      
    }
    
    if (algoname =='pc')
    {
      
      
      tic.clearlog()
      tic("structure runtime")
      
      
      pc_runner <- pc.stable(dataset_truedag,  test = testname, alpha = alphavalue, undirected = FALSE)
      #graphviz.plot(pc_runner, layout = "dot")
      
      
      data_matrix <-amat(pc_runner)
      
      learnt_dag_amat <- matrix(nrow=ncol(dataset),ncol=ncol(dataset))
      learnt_dag_amat[learnt_dag_amat!=0]<-0
      learnt_dag_amat[data_matrix==1 & t(data_matrix)==0]<-2
      learnt_dag_amat[data_matrix==0 & t(data_matrix)==1]<-3
      learnt_dag_amat[data_matrix==1 & t(data_matrix)==1]<-3
      
     
      
      ############## STOP #########################
      toc(log = TRUE, quiet = TRUE)
      log.txt <- tic.log(format = TRUE)
      log2.txt <- rbind(log2.txt, log.txt)
      #############################################
      
      
      write.csv(log2.txt, paste0("output/",set,"/runtime_SOS_dyn_",algoname,"_",alphavalue,"_",set,".csv"), row.names=F)
      
      
      
    }
  }
}


#exportamattoCSV(learnt_dag_amat,paste0("output/",set,"/DAG_SOS_dyn_",algoname,"_",scorename,"_",set,"_kmean.csv"))
