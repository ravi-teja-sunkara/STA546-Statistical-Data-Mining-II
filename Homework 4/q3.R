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

