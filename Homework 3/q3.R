rm(list = ls())

library(gRain)
library(Rgraphviz)
library(gRbase)
library(ggm)

q3 <- list(~A,~B,~C|A,~D|A:B,~E|B,~F|C:A:E,~G|D:E,~H|F:G)
dag_q3 <- dagList(q3)
plot(dag_q3)

### Inquiring D-separation
dSep(as(dag_q3,"matrix"),"C","G", NULL)
dSep(as(dag_q3,"matrix"),"C","E", NULL)
dSep(as(dag_q3,"matrix"),"C","E",c("G"))
dSep(as(dag_q3,"matrix"),"A","G",c("D","E"))
dSep(as(dag_q3,"matrix"),"A","G",c("D"))