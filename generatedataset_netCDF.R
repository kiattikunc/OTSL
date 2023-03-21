library(ncdf4)
library(lattice)
library(RColorBrewer)
library(arules)
ncfname <- paste("Input/air.mon.mean", ".nc", sep="")
colname_lon= c()
#ncfname <- paste("Input/air.sig995.mon.1981-2010.ltm", ".nc", sep="")
 # note: tmp means temperature (not temporary)
ncin <- nc_open(ncfname)
print(ncin)

lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
head(lon)
m <- 1
tmp_slice <- tmp_array[,,m]
levelplot(tmp_array[,,10] ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
          col.regions=(rev(brewer.pal(10,"RdBu"))), main="Mean July Temperature (C)")

lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)
head(lat)

print(c(nlon,nlat))

time <- ncvar_get(ncin,"time")
time

tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

tunits


tmp_array <- ncvar_get(ncin,"air")
dlname <- ncatt_get(ncin,"air","long_name")
dunits <- ncatt_get(ncin,"air","units")
fillvalue <- ncatt_get(ncin,"air","_FillValue")
dim(tmp_array)
tmp_array[1,1,1]
tmp_array2 <- array(0, c(36,18,900))
array
ncvar_get(ncin,"air")
for (a in 1:36)
{
   for (b in 1:18)
   {
          {
            tmp_array2[a,b,] =tmp_array[4*(a-1)+1,4*(b-1)+1,]
          }
    }
}
tmp_array2[2,2,101]
tmp_array[5,5,101]
df_vector <- matrix(tmp_array2,ncol=dim(tmp_array2)[1]*dim(tmp_array2)[2], nrow=dim(tmp_array2)[3], byrow=TRUE)
df <- data.frame(matrix(tmp_array2,ncol=dim(tmp_array2)[1]*dim(tmp_array2)[2], nrow=dim(tmp_array2)[3], byrow=TRUE))
i=1
for (a in 1:36)
{
 
if (a==1) colname_lon = paste0('lon',lon[i])
if (a!=1)
  { 
  i=i+4
  colname_lon = c(colname_lon,paste0('lon',lon[i]))

  }

}


j=1
for (b in 1:18)
{
  
  if (b==1) colname_lat = paste0('lat',lat[j])
  if (b!=1)
  { 
    j=j+4
    colname_lat = c(colname_lat,paste0('lat',lat[j]))
    
  }
  
}
c=1
colname=c()
for (a in 1:36)
{
  for (b in 1:18)
  {
    colname[c]= paste0(colname_lon[a],colname_lat[b])
     c=c+1
  }
}
colnames(df)<-colname
library(bnlearn)
weather = bnlearn::discretize(df, method = 'hartemink', breaks = 4, ibreaks = 10)
write.table(weather, paste0("input/trainingData_weather_kmean.csv"), sep=",", col.names=TRUE, row.names=FALSE)

weather = arules::discretize(df_vector, method="frequency",categories=4)
weather2 = matrix(weather,ncol=dim(tmp_array2)[1]*dim(tmp_array2)[2], nrow=dim(tmp_array2)[3], byrow=TRUE)
weather2 = data.frame(weather2)
colnames(weather2)<-colname
write.table(weather2, paste0("input/trainingData_weather_kmean.csv"), sep=",", col.names=TRUE, row.names=FALSE)


library(raster)
r1 <- raster(ncfname, level=1)
r2 <- raster(ncfname, level=2)

plot(r1)
image(r1) 
s <- stack(r1, r2)
plot(s)


