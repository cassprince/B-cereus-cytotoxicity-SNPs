## RF Regression - SNPs and group
## 01.22.24 

##import data with cytotoxicity, phylogenetic group, and SNPs
data <- read.csv("~/Desktop/PhD/biomarker/SNP_hits_sheet.csv", row.names = 1)

# Count NA values in each row and add the result as a new column
data$count_NA <- apply(data, 1, function(row) sum(is.na(row)))

##1. change "NA" as "2" and run all samples (DO NOT run Ln 11 - 12 if decide this approach)
data[is.na(data)] <- "2"

##2. Remove isolates with at least one "NA" (DO NOT run Ln 8 if decide this approach)
data_no_NA <- subset(data, count_NA == 0)
data <- data_no_NA

##Remove NA count
data <- data[,-24]
for (i in 2:23){
  data[,i] <- as.factor(data[,i])
}
head(data)
str(data)

## Load required packages
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3viz)
library(paradox)

set.seed(3)

## Define data (i.e., outcome variable, input variable, data object)
task <- as_task_regr(Average.Cell.Viability~., data=data)
## Define learner (Random forest regression), and parameters to tune (i.e., mtry and ntrees)
learner <- lrn("regr.ranger", 
               mtry=to_tune(p_int(2,ncol(data)-1)),
               num.trees = to_tune(1e1, 1e4,logscale=FALSE),
               importance = "permutation")
## Define tuning parameters (using 10-fold cross validation to minimize mean squared error)
instance = ti(task= task,
              learner = learner,
              resampling = rsmp("cv", folds=10),
              measures = msr("regr.mse"),
              terminator = trm("none"))
tuner = tnr("grid_search", resolution = 5)
tuner

## Hyperparameter tuning
tuner$optimize(instance)

## Check tuning results 
as.data.table(instance$archive)[, .(mtry, num.trees, regr.mse, batch_nr, resample_result)]
autoplot(instance, type = "surface")

## Running model and export variable importance
split = partition(task, ratio=0.7)
measure = msr("regr.mse")
learner_fin = lrn("regr.ranger",
                  mtry = instance$result_x_domain$mtry,
                  num.trees = instance$result_x_domain$num.trees,
                  importance = "permutation")
learner_fin$train(task, split$train)
prediction = learner_fin$predict(task, split$test)
imp <- as.data.frame(learner_fin$importance())
colnames(imp) <- "importance"
library(ggplot2)
p <- ggplot(data=imp, aes(x=reorder(rownames(imp), importance), y=importance)) + geom_bar(stat="identity")
p + coord_flip() + theme_minimal() +
  xlab("Variable (phylogenetic group and SNPs") +
  ylab("Permutational Variable Importance")
