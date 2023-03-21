
exportamattoCSV <- function(amat,path){


# <- matrix(sample(1:9,9),ncol=3,nrow=3)
#path <-'output/MMPC/1.csv'
ID <-1

export.data <- data.frame()

  for (row in rownames(amat)) {

    for (col in colnames(amat)) {

      if(c(amat[row,col])==2 & c(amat[col,row])==3){

        from <-row
        direction <-"->"
        to <-col

        export.newdata <-data.frame(ID,from,direction,to)

        ID<-ID+1
        export.data <- rbind(export.data,export.newdata)
      }

      else if(c(amat[row,col])==2 & c(amat[col,row])==1){

        from <-row
        direction <-"o->"
        to <-col

        export.newdata <-data.frame(ID,from,direction,to)

        ID<-ID+1
        export.data <- rbind(export.data,export.newdata)

      }
      else if(c(amat[row,col])==2 & c(amat[col,row])==2){

        from <-row
        direction <-"<->"
        to <-col

        export.newdata <-data.frame(ID,from,direction,to)

        ID<-ID+1
        export.data <- rbind(export.data,export.newdata)

      }
      else if(c(amat[row,col])==3 & c(amat[col,row])==3){

        from <-row
        direction <-"-"
        to <-col

        export.newdata <-data.frame(ID,from,direction,to)

        ID<-ID+1
        export.data <- rbind(export.data,export.newdata)

      }
      else if(c(amat[row,col])==1 & c(amat[col,row])==1){

        from <-row
        direction <-"o-o"
        to <-col

        export.newdata <-data.frame(ID,from,direction,to)

        ID<-ID+1
        export.data <- rbind(export.data,export.newdata)

      }
      else if(c(amat[row,col])==1 & c(amat[col,row])==0){

        from <-row
        direction <-"->"
        to <-col

        export.newdata <-data.frame(ID,from,direction,to)

        ID<-ID+1
        export.data <- rbind(export.data,export.newdata)

      }




    }

}
  #sum( amat_pcalG==3)
  colnames(export.data) <- c("ID","Variable 1","Dependency","Variable 2")

  write.csv(export.data,path,row.names = FALSE,quote=FALSE)


}

