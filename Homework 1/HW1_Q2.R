rm(list=ls())
options(warning = -1)
install.packages('arules')
install.packages('MASS')
install.packages('reshape2')
install.packages('ggplot2')

library(arules)
library(MASS)
library(ggplot2)
library(reshape2)
library(ggplot2)

# Data
boston = Boston
head(boston)
summary(boston)

melt_boston <- melt(boston)

# Visualizing
ggplot(melt_boston, aes(x=value, fill=variable)) + geom_histogram(bins=30) + facet_wrap(variable~., scales = 'free')
# hist(boston$crim)

# Dealing with variables
boston$crim <- ordered(cut(boston$crim, c(0, 4, max(boston$crim)), labels=c('Safe', 'Dangerous'))) #based on average
boston$zn <- ordered(cut(boston$zn, c(0,11.36, max(boston$zn)), labels=c('Low', 'High'))) #based on mean/average
boston$indus <- ordered(cut(boston$indus,  c(0,5.19, 18.10, max(boston$indus)), labels=c('Low non-retail area', 'Medium non-retail area', 'High non-retail area'))) # based on quartile values
boston$chas <- ordered(cut(boston$chas, c(0, 0.5, max(boston$chas)), labels=c('Unbound', 'Bounds'))) #only takes 2 values
boston$nox <- ordered(cut(boston$nox, c(0, .4490, 0.6240, max(boston$nox)), labels=c('Low', 'Medium', 'High')))
boston$rm <- ordered(cut(boston$rm, c(0, 5.886, 6.623, max(boston$rm)), labels=c('Few', 'Sufficient', 'Too-many')))
boston$age<- ordered(cut(boston$age, c(0, 45.02, 94.08, max(boston$age)), labels=c('Low', 'Medium', 'High')))
boston$dis <- ordered(cut(boston$dis, c(0, 2.1, 5.1, max(boston$dis)), labels=c('Near-by', 'Average-Distance', 'Far')))
boston$rad <- ordered(cut(boston$rad, c(0, 5.1, 24.01, max(boston$rad)), labels=c('Low', 'Medium', 'High')))
boston$tax <- ordered(cut(boston$tax, c(0, 279, 666, max(boston$tax)), labels=c('Low Tax', 'Average Tax', 'High Tax')))
boston$ptratio <- ordered(cut(boston$ptratio, c(0, 17.40, 20.20, max(boston$ptratio)), labels=c('Less students with more teachers', 'Adequate number of teachers','Teacher Shortage')))
boston$black <- ordered(cut(boston$black, c(0, 375.36, 396.23, max(boston$black)), labels=c('Few Blacks', 'Average number of blacks', 'Many Black people')))
boston$lstat <- ordered(cut(boston$lstat, c(0, 6.95, 16.95, max(boston$lstat)), labels=c('Less Percent', 'Average Percentage', 'High Percentage')))
boston$medv <- ordered(cut(boston$medv, c(0, 17, 25,  max(boston$medv)), labels=c('Cheap', 'Average Price', 'Costly')))
                       
# binary incidence matrix
boston_matrix <- as(boston, 'transactions')
summary(boston_matrix)
itemFrequencyPlot(boston_matrix, support=0.02, cex.names=0.7)

# applying the apriori algorithm
rules <- apriori(boston_matrix, parameter = list(support=0.01, confidence=0.7, minlen=2))
summary(rules)
sample(labels(rules), size=5)


# Part c - Low Crime Area and as close to city as possible as measured by dis
ruleslowCrimeNearbyDis <- subset(rules, subset = lhs %ain% 'crim=Safe' &rhs %in% 'dis=Near-by' & lift>1.25)
summary(ruleslowCrimeNearbyDis)
inspect(head(sort(ruleslowCrimeNearbyDis, by='confidence'), n=10))
# ruleslowCrimeNearbyDis <- subset(rules, subset = rhs %ain% 'crim=Safe' &lhs %in% 'dis=Near-by' & lift>1.25)
ruleslowCrimeNearbyDis <- subset(rules, subset = lhs %ain% c('crim=Safe', 'dis=Near-by') & lift>1.25)
summary(ruleslowCrimeNearbyDis)
inspect(head(sort(ruleslowCrimeNearbyDis, by='confidence'), n=10))

# Part D - schools with low pupil-teacher ratio
ruleslowptratio <- subset(rules, subset = rhs %in% 'ptratio=Less students with more teachers' & lift>2.5)
summary(ruleslowptratio)
inspect(head(sort(ruleslowptratio, by='confidence', decreasing = TRUE), n=10))

# Part E - Regression
fit_regression <- lm(ptratio~., data=Boston)
summary(fit_regression)
# zn, indus, nox, rad, lstat, medv
#nox highly correlated, next rad, medv, indus and zn
