set.seed(3)
library(kernlab)
library(e1071)
library(ROCR)

data(spam)

# split into train 80%, test 20%
testidx <- which(1:length(spam[,1])%%5 == 0)
train <- spam[-testidx,]
test <- spam[testidx,]

# SVM training
model <- svm(type~.,kernel="linear", data=train, probability=TRUE)
prediction_svm <- predict(model, test, probability=TRUE)
table(test$type, prediction_svm)
svm_scores <- attr(prediction_svm,"probabilities")

# ROCR prediction & plotting
prediction_rocr <- prediction(svm_scores[1:920,1], test$type)
performance_rocr <- performance(prediction_rocr,"tpr","fpr")
plot(performance_rocr)

# AUC
performance_rocr.auc <- performance(prediction_rocr,"auc")
auc <- performance_rocr.auc@y.values
mtext(auc,side=3)

# 10-fold cross validation
folds <- function(x,n) split(x, sort(rep(1:n,len=length(x)))) # function to split data into 10
portion <- folds(spam,10)

# Create a list of 2: prediction and label. Each has 10 lists with all numbers
predictions <- list(1)
labels <- list(1)
results_cross <- list("predictions"=predictions,"labels"=labels)

# SVM in each portion of data
for (i in (1:length(portion))) {
  test <- portion[[i]]
  train <- data.frame()
  train_portion <- portion[-i]
  for (j in (1:length(train_portion))) {
    train <- rbind(train,train_portion[[i]])
  }
  
  model <- svm(type~.,kernel="linear", data=train, probability=TRUE)
  prediction_cross <- predict(model, test, probability=TRUE)
  cross_scores <- attr(prediction_cross,"probabilities")

  # store all results in a list
  results_cross[["predictions"]][[i]] <- cross_scores[1:(length(cross_scores)/2),1]
  results_cross[["labels"]][[i]] <- test$type
}

# plot ROCR with 10 portions
prediction_cross_rocr <- prediction(results_cross$predictions, results_cross$labels)
performance_cross_rocr <- performance(prediction_cross_rocr,"tpr","fpr")

plot(performance_cross_rocr,col="grey82",lty=3)
plot(performance_cross_rocr, lwd=3,avg="vertical",spread.estimate="boxplot",add=T)