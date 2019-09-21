rm(list=ls())

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
