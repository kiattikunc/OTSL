
library(ggplot2)
training_set = paste0("Input/raw data BDeu 10k.csv")
dataset_truedag <- read.csv(training_set,header = TRUE,na.strings=c(""),check.names = FALSE)
colnames(dataset_truedag) = c('F1','Cases','Score')

# ggplot(dataset_truedag, aes(x=Cases, y=Score, fill=F1)) +
#   geom_boxplot() +
#   scale_y_continuous(
#                      breaks = pretty(c(0,1), n = 10),
#                      limits = c(0.1,0.8)) +
#   scale_fill_brewer(palette="Spectral")
     

ggplot(data = dataset_truedag,
       aes(x = Cases, y = Score, fill = F1)) +
  geom_boxplot() + 
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = 4, size = 2,
               show.legend = FALSE)+
  scale_y_continuous(
    breaks = pretty(c(0,1), n = 10),
    limits = c(0.1,0.95)) +
  scale_fill_brewer(palette="Spectral") +
  theme(legend.position="none")
  #theme(legend.position="right")

#theme(legend.position="none")