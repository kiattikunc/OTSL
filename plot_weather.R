library(plotly)

library(dplyr)

# airport locations

air <- read.csv('https://raw.githubusercontent.com/plotly/datasets/master/2011_february_us_airport_traffic.csv')

# flights between airports

flights <- read.csv('https://raw.githubusercontent.com/plotly/datasets/master/2011_february_aa_flight_paths.csv')

flights$id <- seq_len(nrow(flights))


# map projection

geo <- list(
  
  scope = 'world',
  
  projection = list(type = 'azimuthal equal area'),
  
  showland = TRUE,
  
  landcolor = toRGB("gray95"),
  
  countrycolor = toRGB("gray80")
  
)


fig <- plot_geo(locationmode = 'USA-states', color = I("red"))

fig <- fig %>% add_markers(
  
  data = air, x = ~long, y = ~lat, text = ~airport,
  
  size = ~cnt, hoverinfo = "text", alpha = 0.5
  
)

fig <- fig %>% add_segments(
  
  data = group_by(flights, id),
  
  x = ~start_lon, xend = ~end_lon,
  
  y = ~start_lat, yend = ~end_lat,
  
  alpha = 0.3, size = I(1), hoverinfo = "none"
  
)

fig <- fig %>% layout(
  
  title = 'Feb. 2011 American Airline flight paths<br>(Hover for airport names)',
  
  geo = geo, showlegend = FALSE
  
)

fig

NAME <- c("A", "A", "B", "B", "C", "C")
YEAR <- c(2016, 2011, 2016, 2011, 2016, 2011)
YEAR <- as.factor(YEAR)
VALUE <- c(1, 4, 1, 5, 2, 8)
DATA <- data.frame(NAME, YEAR, VALUE)

ggplot(DATA, aes(x=VALUE, y=NAME)) + 
  geom_point(size=5, aes(colour=YEAR)) +
  geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))




library(igraph)
library(ggraph)
library(maptools)
data(wrld_simpl)

# create network
# nodes
#DAG_SOS_dyn_pc_bdeu_weather.csv
actors <- read.csv("Input/node_weather.csv",header = TRUE,na.strings=c(""),check.names = FALSE)
relations <- read.csv(paste0("Input/DAG_SOS_dyn_mmhc_bdeu_weather.csv"),header = TRUE,na.strings=c(""),check.names = FALSE)
              
relations2 <- read.csv(paste0("Input/DAG_SOS_dyn_mmhc_bdeu_weather_kmean.csv"),header = TRUE,na.strings=c(""),check.names = FALSE)

# Convert friends matrix to an igraph object
#g <- graph.edgelist(direct, directed = TRUE)

g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
g2 <- graph_from_data_frame(relations2, directed=TRUE, vertices=actors)

# get lat long coordinates for the layout
lo <- layout.norm(as.matrix(actors[, c("long","lat")]),xmin = -2.06,xmax =2.06,ymin = -1,ymax = 1)
#layout.norm(as.matrix(actors[, c("long","lat")]))
#plot
#plot.igraph(g, layout=lo, color="green",rescale=T, vertex.label= NA)

#
dev.new(width=20, height=12, unit="in")
plot(g2, layout=lo,vertex.size=3,rescale=F,
     vertex.label.dist=1, edge.color="orange", vertex.color="orange", edge.arrow.size=0.5,vertex.label= NA)
#par(new=TRUE)
#plot(g2, layout=lo,vertex.size=3,rescale=F,
#     vertex.label.dist=1, edge.color="red", vertex.color="orange", edge.arrow.size=1,vertex.label= NA)

par(new=TRUE)
plot(wrld_simpl,axes=FALSE,ann=TRUE)






          
          
