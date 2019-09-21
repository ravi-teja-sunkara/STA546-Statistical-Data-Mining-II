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
