################################################################
#######       Quesiton 1                      ##############
################################################################
rm(list=ls())

# installing and loading packages
library(cluster)

# data
data("USArrests")
View(USArrests)
summary(USArrests)
str(USArrests)

# a. Performing hierarchical clustering with complete linkage and Euclidean distance
d <- dist(USArrests, method='euclidean')
dim(as.matrix(d))

hclust_cmpl <- hclust(d, method = 'complete')
plot(hclust_cmpl, hang = -1)


# b. cut at height that results in 3 distinct clusters.
rect.hclust(hclust_cmpl, k=3, border=2:6)
abline(h=136, col='red')

ct <- cutree(hclust_cmpl, h=125) 
si <- silhouette(ct, dist = d)
plot(si)

View(ct) #to view the classifications

# c. Scaling the variables to SD=1 and clustering
USArrests_scaled <- as.data.frame(scale(USArrests))

d_scaled <- dist(USArrests_scaled, method = 'euclidean')
hclust_cmpl_scaled <- hclust(d_scaled, method = 'complete')
plot(hclust_cmpl_scaled, hang = -1)

# d. obtaining 3 clusters on scaled data
rect.hclust(hclust_cmpl_scaled, k=3, border=2:6)
abline(h=4.4, col='red')

ct_scaled <- cutree(hclust_cmpl_scaled, k=3)
si_scaled <- silhouette(ct_scaled, dist = d)
plot(si_scaled)

View(ct_scaled)

# viewing the classifications before and after scaling
table(ct, ct_scaled)

# For 3 clusters, scaled data is performing poorly, so finding optimum cluster value for scaled data
# and plotting the same 
# silhouette method for determining optimal number of clusters
fviz_nbclust(USArrests_scaled, hcut, method='silhouette') #2 is ideal from this and also from dendrogram

plot(hclust_cmpl_scaled, hang = -1)
rect.hclust(hclust_cmpl_scaled, k=2, border=2:3)
abline(h=5.25, col='cyan')



################################################################
#######       Quesiton 2                      ##############
################################################################
rm(list=ls())

# installing and loading packages
library(cluster)

# data
df <- read.csv('./Ch10Ex11.csv', header = FALSE)
df_scaled <- as.data.frame(scale(df))
View(df_scaled)
summary(df_scaled)
str(df_scaled)

# hierarchical clustering to the samples using correlation based distance

# The function cor() compute pairwise correlation coefficients between the columns of the data. We want correlations between
# the gene samples so applying cor() directly on the scaled data

# should we apply transpose or not.. should the correlation for observations or columns??? (The genes(1000) separate the samples?)
# the samples should be separated in space i.e., into clusters

corl <- cor(df_scaled, method = 'pearson') #computing correlation matrix
dist_cor <- as.dist(1 - corl)
dim(as.matrix(dist_cor))

# clustering
hclust_sngl <- hclust(dist_cor, method = 'single')
plot(hclust_sngl, main = 'Single Linkage')
rect.hclust(hclust_sngl, k=2, border=2:3)

hclust_avg <- hclust(dist_cor, method = 'average')
plot(hclust_avg, main = 'Average Linkage')
rect.hclust(hclust_avg, k=3, border=2:4)

hclust_cmplt <- hclust(dist_cor, method = 'complete')
plot(hclust_cmplt, main = 'Complete Linkage')
rect.hclust(hclust_cmplt, k=2, border=2:3)

hclust_w <- hclust(dist_cor, method = 'ward.D')
plot(hclust_w, main = 'Ward.D Method')
rect.hclust(hclust_w, k=2, border=2:3)

# C) which genes differ the most across the two groups - using pca
pca_gene <- prcomp(t(df))
summary(pca_gene)
head(pca_gene$rotation)

total_load <- apply(pca_gene$rotation, 1, sum)
index <- order(abs(total_load), decreasing = TRUE)
index[1:10]



###############################################################
#######       Quesiton 3                     ##############
################################################################
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


###############################################################
#######       Quesiton 4                     ##############
################################################################
rm(list=ls())
install.packages('kohonen')
install.packages('mlbench')

library(kohonen)
library(caret)
library(mlbench)

# Data
data('BreastCancer')
View(BreastCancer)
summary(BreastCancer)
str(BreastCancer)
head(BreastCancer)

BreastCancer <- na.omit(data.frame(BreastCancer))
Class <- BreastCancer$Class
BreastCancer$Class <- NULL
BreastCancer$Id <- NULL

BreastCancer <- data.frame(apply(BreastCancer, 2, function(x) as.numeric(as.character(x))))

BreastCancer.scaled <- scale(BreastCancer)

# fit an SOM
set.seed(121)
som_grid <- somgrid(xdim = 5, ydim = 5, topo = 'hexagonal')
BreastCancer.som <- som(BreastCancer.scaled, grid = som_grid, rlen=1000, mode='batch')  #batch som

codes <- BreastCancer.som$codes[[1]]

plot(BreastCancer.som, main = 'Default SOM plot') #plot of som

plot(BreastCancer.som, type = "changes", main = "Changes Plot")
plot(BreastCancer.som, type = "count", main = "Counts Plot")
plot(BreastCancer.som, type = "mapping", main = "Mapping Type SOM", pchs=20)

coolBlueHotRed <- function(n, alpha = 1){rainbow(n, end=4/6, alpha = alpha)[n:1]}
plot(BreastCancer.som, type = "dist.neighbours", palette.name = terrain.colors)

# component plane plots
for (i in 1:33){
  x11()
  plot(BreastCancer.som, type = "property", property=codes[,i], main = colnames(codes)[i])
}

# clustering
d <- dist(codes)
hc <- hclust(d)

plot(hc)
rect.hclust(hc, k=2, border=2:3)

som_cluster <- cutree(hc, h = 7)

my_pal <- c("green", "blue")
my_bhcol <- my_pal[som_cluster]

graphics.off()

plot(BreastCancer.som, type = "mapping", col = "black", bgcol = my_bhcol)
add.cluster.boundaries(BreastCancer.som, som_cluster)



########## SUPERVISED SOM ###############
data('BreastCancer')
str(BreastCancer)
BreastCancer <- na.omit(data.frame(BreastCancer))
BreastCancer$Id <- NULL
Class <- BreastCancer$Class
BreastCancer$Class <- NULL
BreastCancer <- data.frame(apply(BreastCancer, 2, function(x) as.numeric(as.character(x))))
BreastCancer$Class <- Class


set.seed(121)
training_indices <- sample(nrow(BreastCancer), 600)
BreastCancer.train <- scale(BreastCancer[training_indices, 1:9])

BreastCancer.test <- scale(BreastCancer[-training_indices, 1:9], center = attr(BreastCancer.train, "scaled:center"), scale = attr(BreastCancer.train, "scaled:scale"))

BreastCancer.som2 <- xyf(BreastCancer.train, classvec2classmat(Class[training_indices]), grid=somgrid(5, 5, 'hexagonal'), rlen = 100)

# prediction
Class.prediction <- predict(BreastCancer.som2, newdata = BreastCancer.test, trainX = BreastCancer.train, trainY= BreastCancer$Class[training_indices], whatmap = 1)
table(as.integer(BreastCancer$Class[-training_indices]), Class.prediction$predictions[[2]])
confusionMatrix(BreastCancer$Class[-training_indices], Class.prediction$predictions[[2]])


