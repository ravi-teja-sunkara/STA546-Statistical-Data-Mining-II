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
