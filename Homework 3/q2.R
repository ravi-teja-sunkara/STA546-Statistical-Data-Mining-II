rm(list = ls())

library(gRain)
library(Rgraphviz)
library(gRbase)
library(ggm)

## model 
q2 <- list(~Burglary, ~Earthquake, ~TV, ~Nap,~Johncall|Burglary:Earthquake:TV, ~Marycall|Burglary:Earthquake:Nap:TV)
dag_q2 <- dagList(q2)
plot(dag_q2)
