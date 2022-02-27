set.seed(3)
library(klaR)

data(iris)

sub_iris <- iris[1:100, c(1, 3, 5)]
names(sub_iris) <- c("sepal", "petal", "species")
head(sub_iris)
sub_iris$species <- factor(sub_iris$species)

testidx <- which(1:length(sub_iris[,1])%%5 == 0)
train <- sub_iris[-testidx,]
test <- sub_iris[testidx,] 

nbmodel <- NaiveBayes(species~., data=train, type='raw')
prediction <- predict(nbmodel, test[,-3])
table(prediction$class, test[,3])

score <- prediction$posterior[1:(length(score)/2),1]
actual_class <- test$species == 'setosa'
pred <- prediction(score, actual_class)
nbperf <- performance(pred, "tpr", "fpr")
nbauc <- performance(pred, "auc")
nbauc <- unlist(slot(nbauc, "y.values"))
plot(nbperf)
mtext("AUC is",nbauc,side=3)
mtext(nbauc,side=3)
#legend(0.5,3,c(c(paste("AUC is", nbauc)),"\n"),border="white",cex=1.0, box.col = "white")
