rm(list = ls())

# Installing and loading packages
install.packages("bnlearn")
source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "Rgraphviz", "RBGL"))
install.packages("gRain")
install.packages('ggm')

library(gRain)
library(Rgraphviz)
library(gRbase)
library(ggm)
library(caret)

#################################################################################################
#                               Question - 1                                                    #
#################################################################################################
data('cad1', package = 'gRbase')
names(cad1)

#####################   PART - A  ######################
cad_optimal <- list(~Sex, ~SuffHeartF, ~Smoker|Sex, ~Inherit|Smoker, ~Hyperchol|Smoker:SuffHeartF, ~CAD|Inherit:Hyperchol)
cad_optimal_dag <- dagList(cad_optimal)

##### Inquire about D-separation (all of the below are TRUE)
dSep(as(cad_optimal_dag, 'matrix'),'Sex','Hyperchol', 'Smoker') 
dSep(as(cad_optimal_dag, 'matrix'), 'Sex', 'CAD', c('Inherit', 'Hyperchol'))
dSep(as(cad_optimal_dag, 'matrix'), 'Sex', 'CAD', c('Smoker', 'Hyperchol')) # given Inherit/Smoker and Hyperchol
dSep(as(cad_optimal_dag, 'matrix'), 'Sex', 'Inherit', 'Smoker')
dSep(as(cad_optimal_dag, 'matrix'), 'Sex', 'Inherit', c('Smoker', 'CAD')) # Sex indp of Inher given Smoker or Hyper
dSep(as(cad_optimal_dag, 'matrix'), 'Sex', 'SuffHeartF', 'Smoker')

dSep(as(cad_optimal_dag, 'matrix'), 'Smoker', 'CAD', c('Inherit', 'Hyperchol'))
dSep(as(cad_optimal_dag, 'matrix'), 'Smoker', 'SuffHeartF', NULL)


dSep(as(cad_optimal_dag, 'matrix'), 'SuffHeartF', 'CAD', c('Hyperchol', 'Inherit'))
dSep(as(cad_optimal_dag, 'matrix'), 'SuffHeartF', 'CAD', c('Hyperchol', 'Smoker')) #given Inherit/Smoker and Hyperchol
dSep(as(cad_optimal_dag, 'matrix'), 'SuffHeartF', 'Inherit', c('Hyperchol', 'Smoker'))

dSep(as(cad_optimal_dag, 'matrix'), 'Inherit', 'SuffHeartF', c('CAD', 'Hyperchol', 'Smoker'))
dSep(as(cad_optimal_dag, 'matrix'), 'Inherit', 'Hyperchol', 'Smoker')

### Building the network
pp <- extractCPT(cad1, cad_optimal_dag) # A list of conditional probability tables
cpp <- compileCPT(pp)
pn <- grain(cpp)

querygrain(pn) # to get the CP Tables
table(cad1$Sex) # Female = 47 (No), Male = 189 (Yes) # verifying using table fucntions. The probabilities assigned are same using extractCPT as well

plot(pn) # plot of the graph

#######################   PART - B  ###########################
cad_compile <- compile(pn)
summary(cad_compile)

cad_compile_prop <- propagate(cad_compile)
cad_compile_prop.ev <- setFinding(cad_compile_prop, nodes = c('Sex', 'Hyperchol'), states = c('Female', 'Yes'))

## absorbing
querygrain(cad_compile_prop.ev, nodes = c('SuffHeartF', 'CAD'), type = 'marginal')
querygrain(cad_compile_prop.ev, nodes = c('SuffHeartF', 'CAD'), type = 'joint')
querygrain(cad_compile_prop.ev, nodes = c('SuffHeartF', 'CAD'), type = 'conditional')

## not absorbed
querygrain(cad_compile_prop, nodes = c('SuffHeartF', 'CAD'), type = 'marginal')
querygrain(cad_compile_prop, nodes = c('SuffHeartF', 'CAD'), type = 'joint')
querygrain(cad_compile_prop, nodes = c('SuffHeartF', 'CAD'), type = 'conditional')

# 
print(abs) #  values after absorbing evidence
print(not_abs)

#######################   PART - C  ###########################
sim_c <- simulate(cad_compile_prop.ev, nsim = 5, seed= 121)
View(sim_c)
sim_c

# prediction of class
response = c('Smoker', 'CAD')
predict.grain(cad_compile_prop, response, predictors = setdiff(names(sim_c), response), newdata = sim_c, type = 'class' )
# probability values of predictions
predict.grain(cad_compile_prop, response = c('Smoker', 'CAD'), predictors = setdiff(names(sim_c), response), newdata = sim_c, type = 'distribution' ) 

#######################   PART - D  ###########################
sim_d <- simulate(cad_compile_prop.ev, nsim = 500)
write.table(sim_d, file = 'q1_simulate_d.txt', sep = ' ', col.names = T, row.names = F)

predictions_d <- predict(cad_compile_prop, response = c('Smoker', 'CAD'), predictors = c("Sex","SuffHeartF","Hyperchol","Inherit"), newdata = sim_d, type = 'class' )
predict_smoker <- factor(predictions_d$pred$Smoker)
levels(predict_smoker) <- c(levels(predict_smoker), 'No')

predict_cad <- factor(predictions_d$pred$CAD)
levels(predict_cad) <- c(levels(predict_cad), 'No')

cf_smoker <- confusionMatrix(predict_smoker, sim_d$Smoker)
cf_smoker
cf_cad <- confusionMatrix(predict_cad, sim_d$CAD)
cf_cad


#################################################################################################
#                               Question - 2                                                    #
#################################################################################################
q2 <- list(~Burglary, ~Earthquake, ~TV, ~Nap,~Johncall|Burglary:Earthquake:TV, ~Marycall|Burglary:Earthquake:Nap:TV)
dag_q2 <- dagList(q2)
plot(dag_q2)

#################################################################################################
#                               Question - 3                                                    #
#################################################################################################
q3 <- list(~A,~B,~C|A,~D|A:B,~E|B,~F|C:A:E,~G|D:E,~H|F:G)
dag_q3 <- dagList(q3)
plot(dag_q3)

### Inquiring D-separation
dSep(as(dag_q3,"matrix"),"C","G", NULL)
dSep(as(dag_q3,"matrix"),"C","E", NULL)
dSep(as(dag_q3,"matrix"),"C","E",c("G"))
dSep(as(dag_q3,"matrix"),"A","G",c("D","E"))
dSep(as(dag_q3,"matrix"),"A","G",c("D"))
