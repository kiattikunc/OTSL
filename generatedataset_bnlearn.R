library(bnlearn)

load("Input/hailfinder.rda")
dataset = rbn(bn, 10000)


write.csv(dataset, paste0("Input/trainingData_hailfinder_N_10k.csv"), row.names=F)

load("Input/munin1.rda")
dataset = rbn(new, 10000)


write.csv(dataset, paste0("Input/trainingData_munin_N_10k.csv"), row.names=F)
