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

