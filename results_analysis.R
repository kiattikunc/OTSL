

path="results_all.csv"
results_all <- read.csv(path,header = TRUE,na.strings=c(""),check.name=FALSE)
results_all$Best_config_F1 <- as.numeric(results_all$Best_config_F1)
results_all$Best_config_SHD <-as.numeric(results_all$Best_config_SHD)


#set_list = c('ALARM')
set_list = c('ASIA','ALARM','SPORTS','PROPERTY','hepar2','hailfinder')
algo_list = c('hc','mcmc','mmhc','FGS_tetrad','pc')
#algo_list = c('FGS_tetrad')
tuning_list = c('STAR','SCORE_AIC','SCORE_BIC','OCT')
#tuning_list = c('OCT')
samplesize_list = c(1000,10000)
for (set in set_list)
{ 
   for (algoname in algo_list)
   {
      for (samplesize in samplesize_list)
      {
         for (tuning in tuning_list)
         {
           
           
            if (tuning=='SCORE_AIC'){
               matlab_results_all <-paste0("C:/Users/kiattikun/Documents/GitHub/OCT/output/",set,'/bestconfig_TuningScore_',algoname,'_',set,'_',samplesize,'.csv')
               matlab_results_all <- read.csv(matlab_results_all,header = FALSE,na.strings=c(""),check.name=FALSE)
               best_config = matlab_results_all$V2
            }
            else if (tuning=='SCORE_BIC'){
               matlab_results_all <-paste0("C:/Users/kiattikun/Documents/GitHub/OCT/output/",set,'/bestconfig_TuningScore_',algoname,'_',set,'_',samplesize,'.csv')
               matlab_results_all <- read.csv(matlab_results_all,header = FALSE,na.strings=c(""),check.name=FALSE)
               best_config = matlab_results_all$V1
            }
            else if (tuning=='STAR'){
               matlab_results_all <-paste0("C:/Users/kiattikun/Documents/GitHub/OCT/output/",set,'/bestconfig_STAR_',algoname,'_',set,'_',samplesize,'.csv')
               matlab_results_all <- read.csv(matlab_results_all,header = FALSE,na.strings=c(""),check.name=FALSE)
               best_config = matlab_results_all$V1
            }
            else if (tuning=='SOS'){
               matlab_results_all <-paste0("C:/Users/kiattikun/Documents/GitHub/SamplingOutofSample/output/",set,'/bestconfig_SOS_',algoname,'_',set,'_',samplesize,'.csv')
               matlab_results_all <- read.csv(matlab_results_all,header = FALSE,na.strings=c(""),check.name=FALSE)
               best_config = matlab_results_all$V1
            }
            else{
               matlab_results_all <-paste0("C:/Users/kiattikun/Documents/GitHub/OCT/output/",set,'/bestconfig_',tuning,'_',algoname,'_',set,'_',samplesize,'.csv')
               
               matlab_results_all <- read.csv(matlab_results_all,header = FALSE,na.strings=c(""),check.name=FALSE)
               
               best_config = matlab_results_all$V1
            }
              
               
            for (results_all_index in 1:nrow(results_all))
            {
               if (results_all[results_all_index,1]==tuning && results_all[results_all_index,2]==samplesize  && results_all[results_all_index,5]==set && results_all[results_all_index,7]==algoname)
               { 
                  #F1
                  results_all[results_all_index+best_config-1,53]= (as.numeric(results_all[results_all_index+best_config-1,17])  - as.numeric(results_all[results_all_index,17]))/as.numeric(results_all[results_all_index,17])
                  #SHD
                  results_all[results_all_index+best_config-1,54]= (as.numeric(results_all[results_all_index,18])- as.numeric(results_all[results_all_index+best_config-1,18]))/as.numeric(results_all[results_all_index,18])
                  print('update')   
                  break 
                    
               }
            }
          
         }
     
      }
   }
  
     
}
      
write.csv(results_all,'results_all_output.csv')
