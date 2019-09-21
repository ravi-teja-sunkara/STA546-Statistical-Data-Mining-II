rm(list=ls())

# installing and loading packages
install.packages('cluster')
install.packages('caret')
install.packages('stringi')
install.packages('factoextra')
install.packages('NbClust')

library(cluster)
library(caret)
library(stringi)
library(factoextra)
library(NbClust)

# data
load('./primate.scapulae.rdata')
View(primate.scapulae)
summary(primate.scapulae)
str(primate.scapulae)

# Finding optimal number of clusters
fviz_nbclust(primate.scapulae[, 1:9], hcut, method='silhouette') #finding the best number of clusters
fviz_nbclust(primate.scapulae[, 1:9], hcut, method='wss') + geom_vline(xintercept = 3, linetype=2)

# Gap statistic
set.seed(123)
fviz_nbclust(na.omit(primate.scapulae[, 1:9]), hcut, nstart = 25, method='gap_stat', nboot=500)

#################################################
###################### Clustering ################
##################################################
# from the data we can see that we have 5 classes
# calculating the distance for HC
d <- dist(primate.scapulae[, 1:9]) # not taking the factor variables
dim(as.matrix(d))

# Single-linkage Hierarchical clustering
sl <- hclust(d, method='single')
plot(sl, hang=-1)
rect.hclust(sl, k=5, border=2:7)

ct <- cutree(sl, k=5)
si <- silhouette(ct, dist = d)
plot(si)

ct <- cutree(sl, k=3)
si <- silhouette(ct, dist = d)
plot(si)

confusionMatrix(factor(ct), primate.scapulae$classdigit) #accuracy is just 12.38%

# Average Linking Hierarchical clustering
al <- hclust(d, method='average')
plot(al, hang=-1)
rect.hclust(sl, k=5, border=2:7)

ct <- cutree(al, k=5)
si <- silhouette(ct, dist = d)
plot(si)

ct <- cutree(al, k=3)
si <- silhouette(ct, dist = d)
plot(si)

confusionMatrix(factor(ct), primate.scapulae$classdigit)

# Complete Linkage complete clustering
cl <- hclust(d, method='complete')
plot(cl, hang=-1)
rect.hclust(sl, k=5, border=2:7)

ct <- cutree(cl, k=5)
si <- silhouette(ct, dist = d)
plot(si)

ct <- cutree(cl, k=3)
si <- silhouette(ct, dist = d)
plot(si)

confusionMatrix(factor(ct), primate.scapulae$classdigit)


### K-means
primate.scapulae_k <- na.omit(primate.scapulae[]) #removing missing values as kmeans is giving an error for missing values
dats <- primate.scapulae_k[,1:9]
head(dats)

# silhouette method for determining optimal number of clusters
fviz_nbclust(dats, kmeans, method='silhouette')

km3 <- kmeans(dats, centers = 3, nstart=10) # for 3 clusters
km5 <- kmeans(dats, centers = 5, nstart=10) #for 5 clusters

# how well does the clustering match the labels
confusionMatrix(factor(km3$cluster), primate.scapulae_k$classdigit)
confusionMatrix(as.factor(km5$cluster), primate.scapulae_k$classdigit)
