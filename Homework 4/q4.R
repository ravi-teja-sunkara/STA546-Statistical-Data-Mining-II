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
