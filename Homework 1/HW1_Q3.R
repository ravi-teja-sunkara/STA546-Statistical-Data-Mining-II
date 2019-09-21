rm(list=ls())
install.packages('rpart')
install.packages('rpart.plot')
install.packages('rattle')
library(ElemStatLearn)
library(rpart)
library(rpart.plot)
library(rattle)

#data
data("marketing")
View(marketing)
summary(marketing)
str(marketing)

# Class, training and testing
training_sample = marketing
training_sample$class <- 1


ref_sample <- training_sample
for(i in 1:ncol(ref_sample)){
  ref_sample[,i] = sample(ref_sample[,i], nrow(ref_sample), replace=T)
}

ref_sample$class <- 0

# combined
combined_data <- rbind(training_sample, ref_sample)
rm(ref_sample, training_sample)
str(combined_data)

# Model
model.control <- rpart.control(minbucket = 2, minsplit = 100, xval=10, cp=0.02, maxdepth=4)
fit_combined <- rpart(class~., data=combined_data, method='class', control=model.control)
summary(fit_combined)

# Visualization
fancyRpartPlot(fit_combined)
combined_data[,15] <- NULL
# plot(fit_combined, uniform = T, compress=T)
# text(fit_combined, cex=.8, use.n=T, all=T)
# prp(fit_combined , fallen.leaves = FALSE, type=4, extra=1, varlen=0, faclen=0, yesno.yshift=-1,cex =.5)

#prediction
predicted = predict(fit_combined, combined_data[, -c(15)])
predicted

# Even the decision tree doesn't have any nodes..

# # . From the above model, we can observe that features do not have any predictive power to
# do a classification as it has only single root indicating that.
# . To cross verify this, we can predict the model on the training set itself and observe that
# the probability for every row is one-half for both the classes.
# . Therefore, we can conclude that it doesn't have any predictive power.
