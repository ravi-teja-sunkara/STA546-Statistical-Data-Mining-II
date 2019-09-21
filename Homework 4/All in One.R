rm(list=ls())

#############################################################################
#                   QUESTION 1                                              #
#############################################################################
# Installing and loading packages
pckgs <- c("dplyr", "ggplot2", 'factoextra', 'cluster', 'NbClust', 'corrplot', 'fpc', 'SnowballC', 'lsa', 'proxy', 'kmed', 'gRain', 'ggm', 'gRim', 
           'bnlearn', 'igraph')

install.packages(pckgs, dependencies = TRUE)
install.packages(c('BiocInstaller', 'yaml'), dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rgraphviz")
BiocManager::install("RBGL")

pckgs <- append(pckgs, c('RBGL', 'Rgraphviz'))
lapply(pckgs, library, character.only = TRUE)

parkinsons <- read.delim("./parkinsons_updrs.data", header = TRUE, sep = ",")
View(parkinsons)

parkinsons <- parkinsons %>% select(-`subject.`)

##########################
#       Part 1           #
##########################

parkinsons2 <- parkinsons %>% select(-motor_UPDRS, -total_UPDRS)

clust <- hclust(dist(parkinsons2 %>% as.matrix))

clust %>% plot(labels = FALSE, hang = 0, main = 'Cluster')
clustcut <- cutree(clust, k = 10)

total_UPDRS <- parkinsons$`total_UPDRS` %>% cut(breaks=c(0,10,20,30,40,50,60))
motor_UPDRS <- parkinsons$`motor_UPDRS` %>% cut(breaks=c(0, 10, 20, 30, 40))

table(clustcut, total_UPDRS)
table(clustcut, motor_UPDRS)

kmean6 <- kmeans(parkinsons2,6)
table(kmean6$cluster, total_UPDRS)

kmean4 <- kmeans(parkinsons2,4)
table(kmean4$cluster, motor_UPDRS)


#####################
#     Part 2        #
#####################

parkinsons <- read.delim("./parkinsons_updrs.data", header = TRUE, sep = ",")
#View(parkinsons)
str(parkinsons)
colnames(parkinsons)

parkinsons <- parkinsons %>% select(-`subject.`)

# Performing cuts
factoredparkinsons <- parkinsons
factoredparkinsons$age <- cut(parkinsons$age,breaks=c(30,50,70,100))
factoredparkinsons$sex <-  as.factor(parkinsons$sex)
factoredparkinsons$test_time <- cut(parkinsons$test_time,breaks=c(-5, 50, 100, 150, 220))
factoredparkinsons$total_UPDRS <- cut(parkinsons$total_UPDRS,breaks=c(0,10,20,30,40,50,60))
factoredparkinsons$motor_UPDRS <-cut(parkinsons$motor_UPDRS,breaks=c(0, 10, 20, 30, 40))
factoredparkinsons$`Jitter...` <-  cut(parkinsons$`Jitter...`,breaks=c(0, 0.01, 1))
factoredparkinsons$`Jitter.Abs.` <- parkinsons$`Jitter.Abs.` %>% cut(breaks=c(0, 0.00008, 1))
factoredparkinsons$`Jitter.RAP` <- parkinsons$`Jitter.RAP` %>% cut(breaks=c(0, 0.006, 1))
factoredparkinsons$`Jitter.PPQ5` <- parkinsons$`Jitter.PPQ5` %>% cut(breaks=c(0, 0.007, 1))
factoredparkinsons$`Jitter.DDP` <-  cut(parkinsons$`Jitter.DDP`,breaks=c(0,0.02,1))
factoredparkinsons$`Shimmer` <-  cut(parkinsons$`Shimmer`,breaks=c(0, 0.08, 1))
factoredparkinsons$`Shimmer.dB.` <-  cut(parkinsons$`Shimmer.dB.`,breaks=c(0, 0.7, 2.2))
factoredparkinsons$`Shimmer.APQ3` <-  cut(parkinsons$`Shimmer.APQ3`,breaks=c(0,0.035, 1))
factoredparkinsons$`Shimmer.APQ5` <-  cut(parkinsons$`Shimmer.APQ5`,breaks=c(0,0.04,1))
factoredparkinsons$`Shimmer.APQ11` <-  cut(parkinsons$`Shimmer.APQ11`,breaks=c(0, 0.05,1))
factoredparkinsons$`Shimmer.DDA` <-  cut(parkinsons$`Shimmer.DDA`,breaks=c(0, 0.1,1))
factoredparkinsons$`NHR` <-  cut(parkinsons$`NHR`,breaks=c(0, 0.06, 1))
factoredparkinsons$`HNR` <- cut(parkinsons$`HNR`,breaks=c(0, 14, 40))
factoredparkinsons$`RPDE` <-  cut(parkinsons$`RPDE`,breaks=c(0, 0.55, 1))
factoredparkinsons$`DFA` <-  cut(parkinsons$`DFA`,breaks=c(0.5,0.6,0.7,0.87))
factoredparkinsons$`PPE` <-  cut(parkinsons$`PPE`,breaks=c(0, 0.25,1))


factoredparkinsons <- factoredparkinsons %>% select(-motor_UPDRS)
factoredparkinsons <- factoredparkinsons %>% as.data.frame()
factoredparkinsons<-na.omit(factoredparkinsons)
bn <- mclust(factoredparkinsons)
network <- bn %>% amat %>% as("graphNEL")
plot(network, cex = 25)


blocks <- rep(1,20)
blocks[4] <- 2
blm <- matrix(0, 20, 20)
colnames(blm) <- factoredparkinsons %>% names
rownames(blm) <- factoredparkinsons %>% names
blm[blocks == 2, blocks < 2] <- 1
blist <- data.frame(get.edgelist(as(blm, "igraph")))
names(blist) <- c("from", "to")
bn2 <- factoredparkinsons %>% hc(blacklist=blist)
network.constrained <- as(amat(bn2), "graphNEL")
plot(network.constrained, cex = 2)

bn2_to_gRain <- extractCPT(factoredparkinsons, network.constrained, smooth=0.5)
bn2_to_gRain <- compileCPT(bn2_to_gRain)
bn2_to_gRain <- grain(bn2_to_gRain)

bn2_to_gRain_evidence <- bn2_to_gRain %>% setFinding(nodes = c("total_UPDRS"), states = c("(50,60]"))

with_ev <- querygrain(bn2_to_gRain_evidence, nodes = c("Jitter...", "Jitter.Abs.", "Jitter.RAP", "Jitter.PPQ5", "Jitter.DDP"), type = "marginal")
no_ev <- querygrain(bn2_to_gRain, nodes = c("Jitter...", "Jitter.Abs.", "Jitter.RAP", "Jitter.PPQ5", "Jitter.DDP"), type = "marginal")

with_ev
no_ev


#############################################################################
#                   QUESTION 2                                              #
#############################################################################
# Installing and loading packages
pckgs <- c('magrittr', 'dplyr', 'arules', 'gRain', 'Rgraphviz', 'gRbase', 'ggm')
install.packages(pckgs, dependencies = TRUE)
lapply(pckgs, library, character.only = TRUE)

trainingdata = read.csv("./train.csv", header = TRUE)
testingdata = read.csv("./test.csv", header = TRUE)
genderdata = read.csv("./gender_submission.csv", header = TRUE)

testingdata <- cbind(testingdata, genderdata)
testingdata <- testingdata[,-12]
titanic_data <- rbind(trainingdata,testingdata)
rm(testingdata)
rm(trainingdata)

titanic_data <- na.omit(titanic_data)
titanic_data <- as.data.frame(titanic_data)

titanic_data[["Age"]] = ordered(cut(titanic_data[["Age"]], c(0,13,59,80)), labels = c("child", "adult","senior"))
titanic_data[["Survived"]] = factor(titanic_data$Survived, labels = c("yes","no"))
titanic_data[["Pclass"]] = as.factor(titanic_data$Pclass)
titanic_data$Parch <- as.factor(titanic_data$Parch)
titanic_data[["SibSp"]] <- as.factor(titanic_data$SibSp)



titanic_data[["PassengerId"]] = NULL
titanic_data[["Ticket"]] = NULL
titanic_data[["Name"]] = NULL
titanic_data[["Cabin"]] = NULL
titanic_data[["Embarked"]] = NULL
titanic_data[["Fare"]] = NULL
titanic_data = titanic_data[,-7]

titanic_data_matrix <- as(titanic_data, "transactions")

summary(titanic_data_matrix)

###  Binary Incidence Matrix  ###

x11()
itemFrequencyPlot(titanic_data_matrix, support = 0, cex.names = 0.8)
itemFrequencyPlot(titanic_data_matrix, support = 0.1, cex.names = 0.8)

rules <- apriori(titanic_data_matrix, parameter = list(support = 0.01, confidence = 0.6))
summary(rules)
rules
inspect(rules[1:10])


### Survived Passengers  ####

rules<-apriori(titanic_data_matrix, parameter=list(support=0.01,confidence=0.6))

rules_survived <- subset(rules, subset = rhs %in% "Survived=yes" & lift>1)
summary(rules_survived)
survived <- inspect(rules_survived)
survived$lhs
inspect(head(sort(rules_survived,by="confidence"),n=10))
inspect((sort(rules_survived,by="lift")))

#### Dead passengers #####

rules<-apriori(titanic_data_matrix, parameter=list(support=0.01,confidence=0.6))

rules_dead <- subset(rules, subset = rhs %in% "Survived=no" & lift>1)
summary(rules_dead)
dead <- inspect(rules_dead)

inspect(head(sort(rules_dead,by="confidence"),n=10))
inspect(head(sort(rules_dead,by="lift"),n=20))


#############################################################################
#                   QUESTION 3                                             #
#############################################################################
### Loading Required Libraries ###
library(igraph)
library(Rgraphviz)

### Webgraph A ###
set.seed(123)
nodes <- data.frame( names = c("A","B","C","D","E","F"))
relations <- data.frame(from = c("B","B","C","D","D","E","F"),  
                        to = c("C","E","A","B","E","D","C"))
g <- graph.data.frame(relations, directed = TRUE, vertices = nodes)
plot(g,  col = 'blue')


### Page Rank Algorithm ###
damping = c(0.05, 0.25, 0.5, 0.7, 0.99)
for (i in damping){
  i
  pr <- page.rank(g, damping = i)
  print(pr$vector)
}

### Webgraph B ###
nodes <- data.frame( names = c("A","B","C","D","E","F","G","H"))
relations <- data.frame(
  from = c("B","C","D","E","F","G","H"),  
  to = c("A","A","B","B","C","C","C")  
)
g <- graph.data.frame(relations, directed = TRUE, vertices = nodes)
plot(g)


#### Page Rank Algorithm ####
pr <- page.rank(g, damping = 0.15)
print(pr$vector)


#############################################################################
#                   QUESTION 4                                              #
#############################################################################
# Loading libraries
install.packages('kohonen')
library(glasso)
library(gRbase)
library(gRim)
library(gRain)
library(graph)
library(kohonen)

data("state")

# pre process the data
data_us<-data.frame(state.x77)
ncol(data_us)
colnames(data_us)

par(mfrow=c(3,3))
hist(data_us[,1],main='Population')
hist(data_us[,2],main='Income')
hist(data_us[,3],main='Illiteracy')
hist(data_us[,4],main='Life Expectation')
hist(data_us[,5],main='Murder')
hist(data_us[,6],main='HS.Grad')
hist(data_us[,7],main='Frost')
hist(data_us[,8],main='Area')

x11()
plot(data_us)

#Look at partial correlation
?cov.wt
state_us_cov<- cov.wt(data_us,method='ML')
state_us_cov$cov

partial_us_cor<-cov2pcor(state_us_cov$cov)
heatmap(partial_us_cor)

# Use the graphical lasso package to "learn GGM's"
ls("package:glasso")
?glasso

state_cov<- state_us_cov$cov
#Estimate a single graph
s1.lasso <- glasso(state_cov, rho=0.1)

names(s1.lasso)

?glasso

s1.lasso$wi

edges_1<- s1.lasso$wi != 0
edges_1

s1.lasso$wi[1:4,1:4]

diag(edges_1)<-FALSE
names(data_us)
edges_1

#Convert for plotting

g1.lasso<-as(edges_1,"graphNEL") 

g1.lasso

nodes(g1.lasso)<-names(data_us)

glasso1.net<- cmod(g1.lasso,data=data_us)
x11()
plot(glasso1.net)

# Estimate over a range of rho's


rhos <- c(1, 5, 7,10,15)
?glassopath

s2.lasso<- glassopath(state_cov, rho=rhos)

s2.lasso$wi
s2.lasso$rholist

graphics.off()
for(i in 1:length(rhos)){
  edges_2<- s2.lasso$wi[,,i] != 0
  diag(edges_2)<-FALSE
  names(data_us)
  #Convert for plotting
  g2.lasso<-as(edges_2,"graphNEL") 
  nodes(g2.lasso)<-names(data_us)
  
  glasso2.net<- cmod(g2.lasso,data=data_us)
  x11()
  plot(glasso2.net)
}

# SOM
data_us_train <- as.matrix(scale(data_us))
state.scaled <- scale(state.x77)
state.som1 <- som(state.scaled, grid=somgrid(5,4,"hexagonal"))
plot(state.som1, main = "State Data")
som_grid <- somgrid(xdim = 5, ydim=10, topo="hexagonal")
plot(som_grid)
som_model <- supersom(data_us_train, grid=som_grid, rlen=100, alpha=c(0.05,0.01), keep.data = TRUE)
x11()
plot(som_model, type = "property", property = getCodes(som_model, 1)[,1], main=names(som_model$data)[7])
x11()
plot(som_model$changes)

#############################################################################
#                   QUESTION 5                                              #
#############################################################################
# to compute a full Euclidean distance matrix for
dis = function(x)
{
  x = as.matrix(x)
  u = apply(x*x,1,sum) %*% matrix(1.0,1,nrow(x))
  sqrt(abs(u + t(u) - 2 * x %*% t(x)))
}

iorder = function(m)
{
  N = nrow(m) + 1
  iorder = rep(0,N)
  iorder[1] = m[N-1,1]
  iorder[2] = m[N-1,2]
  loc = 2
  for(i in seq(N-2,1))
  {
    for(j in seq(1,loc))
    {
      if(iorder[j] == i)
      {
        iorder[j] = m[i,1]
        if(j==loc)
        {
          loc = loc + 1
          iorder[loc] = m[i,2]
        } else
        {
          loc = loc + 1
          for(k in seq(loc, j+2)) iorder[k] = iorder[k-1]
          iorder[j+1] = m[i,2]
        }
      }
    }
  }
  -iorder
}

hc = function(d, method=c("single","complete","average"))
{
  if(!is.matrix(d)) d = as.matrix(d)
  # Pick a clustering function:
  method_fn = switch(match.arg(method),
                     single   = min,
                     complete = max,
                     average  = mean)
  N = nrow(d)
  diag(d)=Inf
  n = -(1:N)                       # Tracks group membership
  m = matrix(0,nrow=N-1, ncol=2)   # hclust merge output
  h = rep(0,N-1)                   # hclust height output
  for(j in seq(1,N-1))
  {
    # Find smallest distance and corresponding indices
    h[j] = min(d)
    
    i = which(d - h[j] == 0, arr.ind=TRUE)
    i = i[1,,drop=FALSE]
    p = n[i]
    p = p[order(p)]
    m[j,] = p
    grp = c(i, which(n %in% n[i[1,n[i]>0]]))
    n[grp] = j
    r = apply(d[i,],2,method_fn)
    
    d[min(i),] = d[,min(i)] = r
    d[min(i),min(i)]        = Inf
    d[max(i),] = d[,max(i)] = Inf
  }
  
  structure(list(merge = m, height = h, order = iorder(m),
                 labels = rownames(d), method = method, 
                 call = match.call(), dist.method = "euclidean"), 
            class = "hclust")
}

# Comparing ours and original hclust on single 
actual_hclust_singl = hclust(dist(USArrests),method="single")
defined_hclust_singl = hc(dis(USArrests), method="single")
plot(actual_hclust_singl, main = 'In-built function Single Linkage')
plot(defined_hclust_singl, main = 'Our defined function Single Linkage') # our defined function output

actual_hclust_avg = hclust(dist(USArrests),method="average")
plot(actual_hclust_avg, main = 'In-built function Average Linkage')
defined_hclust_avg = hc(dis(USArrests), method="average")
plot(defined_hclust_avg, main = 'Our defined function Average Linkage')

actual_hclust_compl = hclust(dist(USArrests),method="average")
plot(actual_hclust_compl)
defined_hclust_compl = hc(dis(USArrests), method="complete")
plot(actual_hclust_compl)