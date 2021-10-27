#### validation of SURFACER clustering
# https://machinelearningmastery.com/machine-learning-in-r-step-by-step/
setwd("F:/projects/Breast/")

load("results/000_spmra.rda")
load("surfacer_2020.rda")
library(caret)
library(kernlab)
library(randomForest)
library(arm)
#library(RSNNS)
#validation_index <- createDataPartition(dataset$Species, p=0.80, list=FALSE)
load("results/000_patient_mapping.rda")
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="plum2"]<-"Basal-enriched"
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="tomato"]<-"Lum1"
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="cyan"]<-"Lum2"
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="seagreen3"]<-"Lum3"
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="darkkhaki"]<-"Mixed"

dataset<-as.data.frame(t(spmra))
identical(rownames(my_patients),rownames(dataset)) ### TRUE
dataset$subtype<-factor(my_patients$`Surface Genes`,levels<-c("Basal-enriched","Lum1","Lum2","Lum3","Mixed"))


## See available algorithms in caret
modelnames <- paste(names(getModelInfo()), collapse=',  ')
print(modelnames)

#### now use the caret package for prediction
# create a list of 80% of the rows in the original dataset we can use for training
validation_index <- createDataPartition(dataset$subtype, p=0.80, list=FALSE)
# select 20% of the data for validation
validation <- dataset[-validation_index,]
# use the remaining 80% of data to training and testing the models
dataset <- dataset[validation_index,]
# dimensions of dataset
dim(dataset)
# [1] 638   2055
dim(validation)
# [1]  271 2055
# list types for each attribute
sapply(dataset, class)
# list the levels for the class
levels(dataset$subtype)
### Summarize class distribution
# summarize the class distribution
percentage <- prop.table(table(dataset$subtype)) * 100
cbind(freq=table(dataset$subtype), percentage=percentage)

#                 freq percentage
# Basal-enriched  112   15.36351
# Lum1            209   28.66941
# Lum2            128   17.55830
# Lum3            164   22.49657
# Mixed           116   15.91221

# summarize the class distribution
percentage <- prop.table(table(validation$subtype)) * 100
cbind(freq=table(validation$subtype), percentage=percentage)

# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"

#### Test several algorithms

# a) linear algorithms
set.seed(3)
fit.lda <- caret::train(subtype~., data=dataset, method="lda", metric=metric, trControl=control)

# Bayesian generalized Linear Models
set.seed(3)
fit.bglm <- caret::train(subtype~., data=dataset, method="bayesglm", metric=metric, trControl=control)

# b) nonlinear algorithms
# kNN
set.seed(3)
fit.knn <- caret::train(subtype~., data=dataset, method="knn", metric=metric, trControl=control)
# c) advanced algorithms
# SVM
set.seed(3)
fit.svm <- caret::train(subtype~., data=dataset, method="svmRadial", metric=metric, trControl=control)

# Random Forest
set.seed(3)
fit.rf <- caret::train(subtype~., data=dataset, method="rf", metric=metric, trControl=control)

#Multi Layer perceptron
set.seed(3)
fit.mlp <- caret::train(subtype~., data=dataset, method="mlp", metric=metric, trControl=control)

# Nearest shrunken centroids
library(pamr)
set.seed(3)
fit.pam <- caret::train(subtype~., data=dataset, method="pam", metric=metric, trControl=control)

# DeepBoost
library(deepboost)
set.seed(3)
fit.deepboost <- caret::train(subtype~., data=dataset, method="deepboost", metric=metric, trControl=control)

# GBM
library(gbm)
set.seed(3)
fit.gbm <- caret::train(subtype~., data=dataset, method="gbm", metric=metric, trControl=control)

# Greedy Prototype Selection	protoclass
library(protoclass)
set.seed(3)
fit.protoclass <- caret::train(subtype~., data=dataset, method="protoclass", metric=metric, trControl=control)

# Neural Networks with Feature Extraction	pcaNNet
set.seed(3)
fit.pcan <- caret::train(subtype~., data=dataset, method="pcaNNet", metric=metric, trControl=control)

# Stabilized Nearest Neighbor Classifier	snn
library(snn)
set.seed(3)
fit.snn <- caret::train(subtype~., data=dataset, method="snn", metric=metric, trControl=control)

# Variational Bayesian Multinomial Probit Regression	vbmpRadial
# library(vbmp)
# set.seed(3)
# fit.vbmpR <- caret::train(subtype~., data=dataset, method="vbmpRadial", metric=metric, trControl=control)

# select best model
# summarize accuracy of models
results <- resamples(list(lda=fit.lda, bgm=fit.bglm, knn=fit.knn, svm=fit.svm, rf=fit.rf,mlp=fit.mlp,pam=fit.pam,
                          gbm=fit.gbm, protoclass=fit.protoclass, pcan=fit.pcan, snn=fit.snn))
summary(results)
save(results,file = "results/001_ML_performance.rda")
# compare accuracy of models
png("plots/000_validation_models.png",w=2000,h=1500, res=300)
dotplot(results)
dev.off()
#summarize best model
# summarize Best Model
print(fit.svm)

# Prediction
# estimate skill of LDA on the validation dataset
predictions <- predict(fit.svm, validation)
cm<-confusionMatrix(predictions, validation$subtype)
library(ztable)
library(magrittr)
options(ztable.type="html")
z=ztable(cm$table) 
library(tidyverse)
z %>% makeHeatmap() %>% print(caption="Table 4. Heatmap Table")
getwd()
z %>% makeHeatmap() %>%
  print(caption="Table 6. Heatmap table with user-defined palette")
