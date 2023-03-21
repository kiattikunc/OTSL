


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
#set_list = c('ASIA','ALARM','SPORTS','PROPERTY','hepar2','hailfinder')

none_value = rep(c('none'), times=20)
samplesize_list = c(1000,10000)
alpha_stat <- 0.05
s = 100

## select BDeu,EBIC,BIC,LL,MI,BDeus,avg_MI-sh, avg_MI tuning_score ## BDeu refers to BDeu-iss, BDeus; Standard BDeue  and EBIC refers to EBIC-iss
tuning_score = c('LL')

#default ='reverse'
default ='reverse'
##setting for mcmc
{
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
  # 
}

all_config <- rbind(hc_config,mcmc_config,pcstable_config,mmhc_config,FGS_config)
avg_perf_all <- read.csv(paste0("output/performance_SOS_dyn_v4_avg.csv"),header = TRUE,na.strings=c(""),check.name=FALSE)
range_data=149
next_rang_first=1
next_rang_last=149
avg_perf_combine <-c()
for (set in set_list)
{ 
  for (n in samplesize_list)
##################### select algorithms with dynamic score AIC/MIC ##########################
{

  {
    
    avg_perf <- avg_perf_all[next_rang_first:next_rang_last,]
    next_rang_first <- next_rang_first+range_data
    next_rang_last <- next_rang_last+range_data
  }
  
  ##################### select algorithms and the hyper-paramaters with dynamic score BDeu/BIC/MIC ##########################
  avg_perf$tuning_F1 <- NA
  avg_perf$tuning_SHD <-  NA
  if (tuning_score =='BDeu') tuning_score_col = 24+3
  if (tuning_score =='EBIC') tuning_score_col = 25+3
  if (tuning_score =='BIC') tuning_score_col = 27+3
  if (tuning_score =='LL') tuning_score_col = 26+3
  if (tuning_score =='MI') tuning_score_col = 31+3
  if (tuning_score =='BDeus') tuning_score_col = 33
  if (tuning_score =='avg_MI') tuning_score_col = 32
  if (tuning_score =='avg_MI-sh') tuning_score_col = 35
  
  for (config_index3 in 1:nrow(all_config))
  {
    
    if (all_config[config_index3,1] == 'mmhc')
    {
      if (default=='')
      {
        default_position = 1
      }
      if (default=='reverse')
      {
        default_position = 21
      }
      
      range_compare=40
      max_range = config_index3+range_compare-1
      range = config_index3:max_range
      max_range2 = max_range+(range_compare)
      range2 = (max_range+1):max_range2
      
      
      optimal_config = which.max(avg_perf[range,tuning_score_col])
      #F1
      avg_perf[range[optimal_config],36] =   (as.numeric(avg_perf[range[optimal_config],14])  - as.numeric(avg_perf[range[default_position],14]))/as.numeric(avg_perf[range[default_position],14])
      
      #SHD
      avg_perf[range[optimal_config],37] =   (as.numeric(avg_perf[range[default_position],15])- as.numeric(avg_perf[range[optimal_config],15]))/as.numeric(avg_perf[range[default_position],15])
      write.csv(optimal_config, paste0('output/',set,'/bestconfig_SOS_v4_',tuning_score,"_",all_config[config_index3,1],'_',set,'_',n,'.csv'), row.names=F)
      
      
      break
      
      
    }
    
  } 
  for (config_index3 in 1:nrow(all_config))
  {
    if (all_config[config_index3,1] == 'hc')
    {
      if (default=='')
      {
        default_position = 1
      }
      if (default=='reverse')
      {
        default_position = 21
      }
      range_compare=40
      max_range = config_index3+range_compare-1
      range = config_index3:max_range
      max_range2 = max_range+(range_compare)
      range2 = (max_range+1):max_range2
      
      optimal_config = which.max(avg_perf[range,tuning_score_col])
      #F1
      avg_perf[range[optimal_config],36] =   (as.numeric(avg_perf[range[optimal_config],14])  - as.numeric(avg_perf[range[default_position],14]))/as.numeric(avg_perf[range[default_position],14])
      
      #SHD
      avg_perf[range[optimal_config],37] =   (as.numeric(avg_perf[range[default_position],15])- as.numeric(avg_perf[range[optimal_config],15]))/as.numeric(avg_perf[range[default_position],15])
    
        write.csv(optimal_config, paste0('output/',set,'/bestconfig_SOS_v4_',tuning_score,"_",all_config[config_index3,1],'_',set,'_',n,'.csv'), row.names=F)
      break
      
      
    }
    
  } 
  
  for (config_index3 in 1:nrow(all_config))
  {
    if (default=='')
    {
      default_position = 1
    }
    if (default=='reverse')
    {
      default_position = 21
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
      avg_perf[range[optimal_config],36] =   (as.numeric(avg_perf[range[optimal_config],14])  - as.numeric(avg_perf[range[default_position],14]))/as.numeric(avg_perf[range[default_position],14])
      
      #SHD
      avg_perf[range[optimal_config],37] =   (as.numeric(avg_perf[range[default_position],15])- as.numeric(avg_perf[range[optimal_config],15]))/as.numeric(avg_perf[range[default_position],15])
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
      avg_perf[range[optimal_config],36] =   (as.numeric(avg_perf[range[optimal_config],14])  - as.numeric(avg_perf[range[1],14]))/as.numeric(avg_perf[range[1],14])
      
      #SHD
      avg_perf[range[optimal_config],37] =   (as.numeric(avg_perf[range[1],15])- as.numeric(avg_perf[range[optimal_config],15]))/as.numeric(avg_perf[range[1],15])
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
      avg_perf[range[optimal_config],36] =   (as.numeric(avg_perf[range[optimal_config],14])  - as.numeric(avg_perf[range[1],14]))/as.numeric(avg_perf[range[1],14])
      
      #SHD
      avg_perf[range[optimal_config],37] =   (as.numeric(avg_perf[range[1],15])- as.numeric(avg_perf[range[optimal_config],15]))/as.numeric(avg_perf[range[1],15])
      write.csv(optimal_config, paste0('output/',set,'/bestconfig_SOS_v4_',tuning_score,"_",all_config[config_index3,1],'_',set,'_',n,'.csv'), row.names=F)
      
      break
    }
  } 
  
  
  #write.csv(avg_df, paste0("output/",set,"/graphicalaccuracy_SOS_dyn_",algoname,"_",scorename,"_avg",length(round),"_",set,"_",n,".csv"), row.names=F)
  write.csv(avg_perf, paste0("output/",set,"/performance_SOS_dyn_v4_(",tuning_score," ",default,")_avg_",length(round),"_",set,"_",n,".csv"), row.names=F)
  avg_perf_combine <-rbind(avg_perf_combine,avg_perf)
  }
  
}

write.csv(avg_perf_combine,paste0("output/performance_SOS_dyn_v4_(",tuning_score," ",default,")_avg_",length(round),".csv"), row.names=F)