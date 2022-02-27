set.seed(3)
library(kernlab)
library(e1071)
library(ROCR)

setwd("C:/Users/toshiba/Desktop/Sofía/Spring 2022/FCBB/Week 3")
dataset <- read.table("p03_Kato_P53_mutants_200.txt", header=TRUE)
dataset$CLASS <- factor(dataset$CLASS)

# 10-fold cross validation
folds <- function(x,n) split(x, sort(rep(1:n,len=length(x)))) # function to split data into 10
num <- 3
portion <- folds(dataset[,-1],num)

# Create a list of 2: prediction and label. Each has 10 lists with all numbers
predictions <- list(1)
labels <- list(1)
results_cross <- list("predictions"=predictions,"labels"=labels)

# SVM in each portion of data
for (i in (1:num)) {
  test <- portion[[i]]
  train <- data.frame()
  train_portion <- portion[-i]
  for (j in (1:length(train_portion))) {
    train <- rbind(train,train_portion[[j]])
  }
  
  model <- svm(CLASS~.,kernel="linear", data=train, probability=TRUE)
  prediction_cross <- predict(model, test, probability=TRUE)
  cross_scores <- attr(prediction_cross,"probabilities")
  
  # store all results in a listi=
  results_cross[["predictions"]][[i]] <- cross_scores[1:(length(cross_scores)/2),2]
  results_cross[["labels"]][[i]] <- test$CLASS
}

# plot ROCR with 10 portions
prediction_cross_rocr <- prediction(results_cross$predictions, results_cross$labels)
performance_cross_rocr <- performance(prediction_cross_rocr,"tpr","fpr")
# get opposite curve?

plot(performance_cross_rocr,col="grey82",lty=3)
plot(performance_cross_rocr, lwd=3,avg="vertical",spread.estimate="boxplot",add=T)

# AUC
performance_cross_rocr.auc <- performance(prediction_cross_rocr,"auc")
auc <- performance_cross_rocr.auc@y.values
mtext(auc,side=3)
