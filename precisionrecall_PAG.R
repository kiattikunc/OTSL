precisionrecall <- function(pag,gtpag){
  
  # modified adjacency matrix for PAG MAG DAG
  # o-o = (1,1)
  # o-> = (1,2)
  # o- = (1,3)
  # - = (3,3)
  # -> = (3,2)
  # <-> = (2,2)
  
  
  nVars<-nrow(gtpag)
  tps <- (((pag-gtpag) ==0) & (t(pag)-t(gtpag) ==0)) & (pag!=0)
  tp <- 0.5*length(tps[tps==TRUE])
  
  
  # o-o vs -> or o-o vs o-> 
  tps_variant <- (pag ==1) & (t(pag) ==1) & (t(gtpag) ==2) | (pag ==1) & (t(pag) ==1) & (gtpag ==2)  | (gtpag ==1) & (t(gtpag) ==1) &  (t(pag) ==2) | (gtpag ==1) & (t(gtpag) ==1) &  (pag ==2)
  
  # - vs -> or *- vs -* 
  tps_undirect <- (pag ==3) & (t(pag) ==3) & ((gtpag) !=3) & (t(gtpag) ==3) | (gtpag ==3) & (t(gtpag) ==3) & ((pag) !=3) & (t(pag) ==3)
  
  tps_bidirect <- ((pag ==2) & (t(pag) ==2) & ((gtpag ==2) & (t(gtpag) ==2)))
  # <-> vs o-> or <-> vs o-o
  tps_o_direct <- (pag ==1 & t(pag) ==2 & gtpag !=1 & t(gtpag) ==2) | (gtpag ==1 & t(gtpag) ==2 & pag !=1 & t(pag) ==2) 
  # o-> vs ->  or  o-> vs <->
  
  tp_undirect <- 0.5*0.5*length(tps_undirect[tps_undirect==TRUE])
  tp_bidirect <- 0.5*length(tps_bidirect[tps_bidirect==TRUE])
  tp_variant <- 0.5*0.5*length(tps_variant[tps_variant==TRUE])
  tp_o_direct <- 0.5*length(tps_o_direct[tps_o_direct==TRUE])
  
  
  
  
  tns <-  gtpag==0 & t(gtpag)==0 & pag==0 & t(pag)==0
  tn<- 0.5*(length(tns[tns==TRUE]) -length(diag(tns)))
  
  
  temp2 <- gtpag!=0
  temp3 <- (gtpag ==2)& t(gtpag== 2)
  nEdgesGtPag <- 0.5*length(temp2[temp2==TRUE])
  nBiGtEdges <- 0.5*length(temp3[temp3==TRUE])
  
  
  
  temp2 <- pag!=0
  nEdgesMag <- 0.5*length(temp2[temp2==TRUE])
  
  temp3 <- (pag ==2)& t(pag== 2)
  nBiEdges <- length(temp3[temp3==TRUE])
  
  temp3 <- (pag ==1)& t(pag== 1)
  nUnEdges <- 0.5*length(temp3[temp3==TRUE])
  
  
  
  
  temp4 <- pag!=0
  temp5 <- (pag ==2)& t(pag== 2)
  temp6 <- (pag ==3)& t(pag== 3)
  temp7 <- (pag ==1)& t(pag== 1)
  nEdgesPag <- 0.5*length(temp4[temp4==TRUE])
  nBiEdges <- 0.5*length(temp5[temp5==TRUE])
  nUnEdges <- 0.5*length(temp6[temp6==TRUE])+tp_undirect
  oEdges <- length(temp7[temp7==TRUE])
  
  
  
  tp<-  tp+tp_o_direct+tp_variant+2*tp_undirect
  
  
  a <- nEdgesGtPag
  
  i <-(nVars*(nVars-1)/2)-a
  
  fp <- i -tn
  fn <- nEdgesGtPag - tp 
  BSF <- 0.5*(tp/a+tn/i-fp/i-fn/a)
  #nCorrectEdges = nrow(((pag-gtpag)==0)& ((pag'-gtpag')==0) & ~~pag);
  
  
  precision <- tp/(nEdgesPag)
  recall <- tp/nEdgesGtPag
  f1 <- 2*(precision*recall)/(precision+recall)
  shd <- sum(pag !=gtpag &  t(pag) !=t(gtpag))/2 +tp_undirect
  return_list <- list('precision'=precision,'recall'=recall,'BSF'=BSF,'f1'=f1, 'shd'=shd ,'tp'=tp,'fp'=fp,'fn'=fn,'tn'=tn,'nEdgesGtPag'=nEdgesGtPag,'a'=a,'i'=i,'nEdgesPag'=nEdgesPag,'tp_undirect'=tp_undirect,'tp_bidirect'=tp_bidirect,'tp_o_direct'=tp_o_direct)
  
}












## SHD as defined in Tsamardinos et al. (2006)
##  True   Estimated   Penalty
##  *-*                  1  (a1)
##         *-*           1  (a1)
##  ->     *-           1   (a3)
##  -*      <-          1   (a3)
##  ->     <-          1    (a4)               
##  <->     *-          1  (a5)     
##   *-     <->          1  (a5)  