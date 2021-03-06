## @knitr loading
library(knitr)
library(gdata)
library(ranger)
library(caret)
library(foreach)
library(doMC)
library(e1071)
library(pROC)
library(grid)
library(gridExtra)
library(stats)



## done in class, loads all the input data
load("./DenBoerData_loaded.Rdata")

## register a paralell backend and define number of Threads
## Note the the numberThreads variable gets passed later to the ranger function
## this is neccessary, because ranger uses all available threads as a default and I
## did not intend to overload the cluster.
numberThreads <- 6
doMC::registerDoMC(numberThreads)
## Here I store all the different parameters we want to use for the different
## classification methods
numberOfTrees <- c(200, 500, 1000)
kernels <- c("radial", "linear", "polynomial", "sigmoid")
## for reusability provide the name of the class variable for this dataset
## I also save the minimum number allowed cases per class in an extra variable
classVariable <- "sample.labels"
minNumberofCasesPerClass <- 30
## this is a helper variable for the feature selection 
## I only allow the first 200 genes of the random forest variable importance
## to be considerred for my final selection
preFiltered <- 200
## Expressive Naming for the result table later on
parameters <- c("RF 200 Trees all Genes",
                "RF 500 Trees all Genes",
                "RF 1000 Trees all Genes",
                "RF 200 Trees RF Selection",
                "RF 500 Trees RF Selection",
                "RF 1000 Trees RF Selection",
                "RF 200 Trees Variance Selection",
                "RF 500 Trees Variance Selection",
                "RF 1000 Trees Variance Selection",
                "SVM Radial Kernel all Genes",
                "SVM Linear Kernel all Genes",
                "SVM Ploynomial Kernel all Genes",
                "SVM Sigmoid Kernel all Genes",
                "SVM Radial Kernel RF Selection",
                "SVM Linear Kernel RF Selection",
                "SVM Ploynomial Kernel RF Selection",
                "SVM Sigmoid Kernel RF Selection",
                "SVM Radial Kernel Variance Selection",
                "SVM Linear Kernel Variance Selection",
                "SVM Ploynomial Kernel Variance Selection",
                "SVM Sigmoid Kernel Variance Selection")


####  Define Local directories and files ####
## Find the home directory
myHome <- Sys.getenv("HOME")

## Define the main directory for the  data and results
mainDir <- file.path(myHome, "uni/4_sem/biostatistics/report")
dir.create(path = mainDir, showWarnings = FALSE)
message("Main dir: ", mainDir)

## Define a file where we wills store the memory image
memImageFile <- file.path(mainDir, "DenBoerData_loaded.Rdata")

## Define the dir where we will download the data
destDir <- file.path(mainDir, "data")
message("Local data dir: ", destDir)

## Define a file where we will store all results
resultDir <- file.path(mainDir, "resultDir")
message("Result direcotry", resultDir)

## check if pheno data and expression data have same order of patients
check <- colnames(exprTable) == rownames(phenoTable)
message("Do Pheno and Expressiondata have the same order of patients")
message(all(check = TRUE))

## Transposing the exprTable do that it has the right format for the package
exprTableTransposed <- t(exprTable)

## this prevents an error in the ranger function due to ilegal names
colnames(exprTableTransposed) <- make.names(colnames(exprTableTransposed))


## @knitr preFilter
#### Classs Filtering ####
## beofre filtering I dsiplay the class frequency in a table
classesbeforeFiltered <- as.data.frame(sort(table(phenoTable$Sample.title),
                                            decreasing = TRUE))
colnames(classesbeforeFiltered) <- c("Class", "Freq")
kable(classesbeforeFiltered, caption = "Table 1: Class Sizes before Filtering")


## @knitr filter
## first we filter the classes and include only classes with more than 30 patients
# Function for class filtering
FilterClass <- function(minimumNumberofCases){
  classTable <- table(phenoTable[, classVariable])
  releventClasses <- names(which( classTable > minimumNumberofCases))
  helper <- phenoTable[, classVariable] %in% releventClasses
  return(helper)
}

## gets a logical vector of which samples to include
relevantPatients <- FilterClass(minNumberofCasesPerClass)
phenoTable <- phenoTable[relevantPatients, ]
#remove unused factors
phenoTable[,classVariable] <- factor(phenoTable[,classVariable])
phenoTable$Sample.title <- factor(phenoTable$Sample.title)

# filter also exprTable
exprTableTransposed <- exprTableTransposed[relevantPatients, ]

## for displaying in markdown
classesAfterFiltered <- as.data.frame(sort(table(phenoTable$Sample.title),
                                           decreasing = TRUE))
colnames(classesAfterFiltered) <- c("Class", "Freq")
kable(classesAfterFiltered, caption = "Table 2: Class Sizes after Filtering")
#grid.arrange(tableGrob(as.data.frame(classesAfterFiltered)))


## @knitr prepareSelect
#### Prepare Feature Selection ####
## Function that tunes the mytry parameter of RF

tuneMtry <- function(relevantGenes) {
  ## first we tune the mtry parameter of the random forest model with caret
  ## 10 folds repeat 3 times
  control <- trainControl(method = 'repeatedcv',
                          number = 10,
                          repeats = 3,
                          search = 'random')
  set.seed(123)
  rf.random <- train(y = phenoTable[,classVariable],
                     x = exprTableTransposed[, relevantGenes],
                     method = "ranger",
                     metric = "Accuracy",
                     trControl = control)
  
  tunedMtry <- rf.random$finalModel$mtry
  return (tunedMtry)
}


tunedMtryAllFeatures <- tuneMtry(colnames(exprTableTransposed))
## this is done becuase the tunedMtry parameter is used globally and gets overwritten later by the default value
tunedMtry <- tunedMtryAllFeatures

## Function that selects genes based on the random forest variable importance
## It takes a vector as index variable that indicates which patients should be used
SelectFeaturesWithRandomForest <- function(trainIndex,
                                           numFeatures,
                                           verbose = TRUE) {
  if(verbose) message("Fitting a Random Forest for feature selection")
  fit <- ranger(y = phenoTable[trainIndex, classVariable],
                x = exprTableTransposed[trainIndex,],
                importance = "impurity",
                mtry = tunedMtry,
                num.threads = numberThreads,
  )
  
  if(verbose) message("Finished fitting, now extracting variable importance")
  varImportance <- fit$variable.importance
  selectedGenes <- sort(varImportance, na.last = TRUE, decreasing = TRUE)
  return(selectedGenes[1:numFeatures])
}



## get number of features
getNumberofFeatures <- function(FUN, preFiltered) {
  selectedGenes <- FUN(c(1:nrow(exprTableTransposed)),
                       ncol(exprTableTransposed))
  ## helper for smoothing the plot
  lo <- loess(selectedGenes[1:preFiltered] ~ c(1:preFiltered))
  smoothed = predict(lo)  # actal smoothing
  secondDer <- diff(diff(smoothed))
  maximalChangePoint <- max(secondDer)
  maximalChangeIndex <- match(maximalChangePoint, secondDer)
  numberFeatures <- maximalChangeIndex
  return(list(numberFeatures = numberFeatures,
              smoothedLine = smoothed,
              genes = selectedGenes,
              maximalChangeIndex = maximalChangeIndex
  ))
}



## get the number of genes that I will select later on
helperRFSelection <- getNumberofFeatures(SelectFeaturesWithRandomForest, preFiltered)
numberFeaturesRF <- helperRFSelection[["numberFeatures"]]



## Function for ploting the variable importance
plotVariableImportance <- function(smoothed,
                                   selectedGenes,
                                   maximalChangeIndex,
                                   preFiltered) {
  plotData <- data.frame(`Variable Importance` = selectedGenes[1:preFiltered],
                         Index = c(1:preFiltered),
                         Smoothed = smoothed,
                         check.names = FALSE)
  colors <- c("Sepal Width" = "blue", "Petal Length" = "red",
              "Petal Width" = "orange")
  ggplot(plotData,
         aes(x = Index)) +
    geom_point(aes(y = `Variable Importance`, color = "Empiric")) +
    geom_line(aes(y = smoothed, color = "smoothed"))  +
    geom_vline(aes(xintercept = maximalChangeIndex, color= "Turning Point"),
               size = 1.5)
}

## execute the plotting function
plotVariableImportance(helperRFSelection[["smoothedLine"]],
                       helperRFSelection[["genes"]],
                       helperRFSelection[["maximalChangeIndex"]],
                       preFiltered = preFiltered)


## @knitr plotVariance
## Function for selecting genes based on Variance
SelectGenesByVariance <- function(trainIndex,
                                  numFeatures,
                                  verbose = TRUE) {
  if(verbose) message("Gettin genewise Variance")
  genewiseVar <- apply(exprTableTransposed[trainIndex, ], 2, var)
  sortedVar <- sort(genewiseVar, decreasing = TRUE, na.last = TRUE)
  return(sortedVar[1:numFeatures])
}

helperVarSelection <- getNumberofFeatures(SelectGenesByVariance, 500) 
numberFeaturesVar <- helperVarSelection[["numberFeatures"]]

## execute the plotting function
plotVariableImportance(helperVarSelection[["smoothedLine"]],
                       helperVarSelection[["genes"]],
                       helperVarSelection[["maximalChangeIndex"]],
                       500)

## @knitr select
#### Feature Selection ####
## Selection with leave one out cross validation
SelectLOOCV <- function(FUN,
                        numFeatures,
                        verbose = FALSE) {
  res <- foreach(i = 1:nrow(exprTableTransposed)) %dopar% {
    curVariables <- FUN(-i,
                        numFeatures,
                        verbose)
    return(curVariables)
    
  }
  names(res) <- rownames(phenoTable)
  return(res)
}

## perform actual feature selection and save selctions
loocvSelectionsRF <- SelectLOOCV(SelectFeaturesWithRandomForest, numberFeaturesRF)

## perfomr variance based selection
loocvSelectionsVar <- SelectLOOCV(SelectGenesByVariance, numberFeaturesVar)

## @knitr execute
#### Classifier Functions ####
## Function that classifies patients with random forest
RandomForestClassifier <- function(numberTrees,
                                   testIndex,
                                   trainIndex,
                                   selectedCovariates,
                                   verbose = TRUE) {
  if(verbose) message("Fitting the Random Forest")
  if(verbose) message(paste0("Using ",numberTrees," trees"))
  rf.fit <- ranger(y = phenoTable[trainIndex, classVariable],
                   x = exprTableTransposed[trainIndex, selectedCovariates],
                   num.trees = numberTrees,
                   num.threads = numberThreads,
                   mtry = tunedMtry)
  
  testData <- t(as.data.frame(exprTableTransposed[testIndex,
                                                  selectedCovariates]))
  if(length(testIndex) !=1) {
    testData <- exprTableTransposed[testIndex, selectedCovariates]
  }
  if(verbose) message("Starting prediction based on the fitted model")
  predicted <- predict(rf.fit, testData)
  if(verbose) message("Finished predicting, now returning the predictions")
  predicted$predictions
  return(predicted$predictions)
}



## Function that classifies patients with SVM
SvmClassifier <- function(myKernel,
                          testIndex,
                          trainIndex,
                          selectedCovariates,
                          verbose = TRUE) {
  if(verbose) message("Starting fitting a SVM Mode")
  svm.fit <- svm(y = phenoTable[trainIndex, classVariable],
                 x = exprTableTransposed[trainIndex, selectedCovariates],
                 kernel = myKernel,
                 gamma = 0.1,
                 cost = 10,
                 type = "C-classification")
  # neccessary in case the test index has length one (loocv)
  testData <- t(as.data.frame(exprTableTransposed[testIndex,
                                                  selectedCovariates]))
  if(length(testIndex) !=1) {
    testData <- exprTableTransposed[testIndex, selectedCovariates]
  }
  predicted <- predict(svm.fit, newdata = testData)
  return(predicted)
}



#### Function for LOOCV ####
LOOCV <- function(FUN,
                  parameter,
                  selection,
                  verbose = FALSE) {
  res <- foreach(i = 1:nrow(exprTableTransposed)) %dopar% {
    curVariables <- names(selection[
      rownames(exprTableTransposed)[i]][[1]])
    trainIndex <- c(1:nrow(exprTableTransposed))[-i]
    if(verbose) message("Fitting Loocv")
    res <- FUN(parameter,
               i, -i,
               curVariables,
               verbose)
    res <- droplevels(res)
    names(res) <- NULL
    return(res)
  }
  if(verbose) message('returning')
  names(res) <- rownames(exprTableTransposed)
  return(unlist(res))
}

# for using all genes in classifiers helper variable
selectedGenes <- rep(1, ncol(exprTableTransposed))
names(selectedGenes) <- colnames(exprTableTransposed)
## helper for perfoming the classification with two different sets of selcted genes
allGenesSelected <- rep(list(selectedGenes), nrow(exprTableTransposed))
names(allGenesSelected) <- rownames(exprTableTransposed)
selections <- list(
                   allGenes = allGenesSelected,
                   rfSelection = loocvSelectionsRF,
                   varSelection = loocvSelectionsVar)


## this loop saves all prediction results in a list
## might be a bit over-engineered for our case, but as soon as one want to compare
## more parameters and methods and feature selections this becomes very useful
resultVector <- foreach(selection = names(selections), .combine = "c") %do% {
  selData <- selections[[selection]]
  tunedMtry <- ceiling(sqrt(length(selData[[1]])))
  
  message(tunedMtry)
  
  ## combines all the results with the random forest
  rf.comb <- foreach(numTree = numberOfTrees, .combine = "c") %do% {
    rf.loocv <- LOOCV(RandomForestClassifier, numTree, selData)
    helperLoocvFile <- paste0("rf_loocv_Selection_", selection,
                              "_numTrees_", numTree, ".csv")
    curLoocvFile <- file.path(resultDir, helperLoocvFile)
    write.csv(rf.loocv, curLoocvFile)
    res <- list(numTree = rf.loocv)
    names(res) <- paste0(numTree)
    return(res)
  }
  ## combines all results of the svm
  svm.comb <- foreach (kern = kernels, .combine = "c") %do% {
    svm.loocv <- LOOCV(SvmClassifier, kern, selData)
    helperLoocvFile <- paste0("SVM_loocv_Selection_",
                              selection, "_kernel_", kern, ".csv")
    curLoocvFile <- file.path(resultDir, helperLoocvFile)
    write.csv(svm.loocv, curLoocvFile)
    res <- list(kern = svm.loocv)
    names(res) <- paste0(kern)
    return(res)
  }
  svm.name <- paste0("svm ", selection)
  rf.name <- paste0("rf ", selection)
  res <- c(rfRes = rf.comb, svmRes = svm.comb)
  names(res) <- paste0(selection, names(res))
  return(res)
}


## Reordering the results to have the same methods together
index <- c(1,2,3,8,9, 10, 15, 16, 17,4,5,6,7,11,12,13,14, 18, 19, 20, 21)
resultVector <- resultVector[order(index)]

names(resultVector) <- parameters
response <- phenoTable[, classVariable]
response <- factor(response, levels = rev(levels(response)))
names(response) <- rownames(phenoTable)



## @knitr displayResult
#### Functions for displaying the results ####
## Functions for plotting a heatmap of a confusion table
plotResults <- function(res, response) {
  cm <- confusionMatrix(res, response)
  cm.table <- as.data.frame(cm$table)
  cm.stats <-data.frame(Stats = cm$overall)
  cm.stats$Stats <- round(cm.stats$Stats,2)
  cm.percentage <- as.data.frame(prop.table(cm$table))
  cm.table$Perc <- round(cm.percentage$Freq*100,2)
  
  cm.plot <- ggplot(data = cm.table, aes(x = Prediction , y =  Reference, fill = Freq)) +
    geom_raster(aes(fill = Freq)) +
    geom_text(aes(label = paste("",Freq,",",Perc,"%")),
              color = 'white', size = 3, fontface = "bold") +
    theme_light()
  cm.stats <- cm.stats[-c(5,7), ,drop = F]
  cm.statsTable <-  tableGrob(cm.stats)
  
  grid.arrange(cm.plot, cm.statsTable,nrow = 1, ncol = 2, widths=c(0.7, 0.3))
  #top = textGrob(paste0("Confusion Table Heatmap \n",title), gp = gpar(fontsize=20,font=1)))
  
}


plotResults(resultVector[[3]], response)
## @knitr rfSecond
plotResults(resultVector[[7]], response)

## @knitr displaySVM
plotResults(resultVector[[14]], response)
## @knitr displaySVM2
plotResults(resultVector[[18]], response)
## @knitr appendix
plotResults(resultVector[[5]], response)
plotResults(resultVector[[6]], response)

plotResults(resultVector[[8]], response)
plotResults(resultVector[[9]], response)
plotResults(resultVector[[10]], response)
plotResults(resultVector[[11]], response)
plotResults(resultVector[[12]], response)
plotResults(resultVector[[13]], response)


plotResults(resultVector[[1]], response)
plotResults(resultVector[[2]], response)
## @knitr overviewResultsRF
getResultOverview <- function (results) {
  evaluationResults <- foreach(i = c(1:length(results)), .combine = 'rbind') %do% {
    res <- results[[i]]
    cm <- confusionMatrix(res, response)
    check <- names(response) == names(res)
    message("Patients in Same Order:")
    message(all(check,TRUE))
    cm.table <- cm$table
    hits <- sum(diag(cm.table))
    errors <- sum(cm.table) - hits
    misclError <-  errors / sum(cm.table)
    accuracy <- 1 - misclError
    ci.upper <- accuracy + 1.96 *
      sqrt( (accuracy * (1 - accuracy)) / nrow(exprTableTransposed))
    ci.lower <- accuracy - 1.96 *
      sqrt( (accuracy * (1 - accuracy)) / nrow(exprTableTransposed))
    ci.interval <- paste0(round(ci.lower,4), "-", round(ci.upper,4))
    message("Misclassifiction Error for current predictions:")
    message(round(misclError,4))
    return(data.frame(MER = round(misclError,4),
                      Accuracy = round(accuracy,4),
                      `95% CI` = ci.interval,
                      check.names = FALSE))
  }
  rownames(evaluationResults) <- names(results)
  return(evaluationResults)
}


resultTable <- getResultOverview(resultVector)
rfResultTable <- resultTable[1:9, ]
kable(rfResultTable, caption = "Table 3: Overview of the Random Forest Classifier Performance with different Parameters")

## @knitr overviewResultsSVM
svmResultTable <- resultTable[10:21, ]
kable(svmResultTable, caption = "Table 4: Overview of the SVM Classifier Performance with different Parameters")


## @knitr inDepth
cm.svmBest <- t(confusionMatrix(response, resultVector[[12]])$byClass)
kable(round(cm.svmBest,4), caption = "Table 5: Classwise statistics of the SVM Classifier trained with the polynomial kernel and allgenes. And the Random Forest Classifier trained with 1000 trees and all Genes")
#grid.arrange(tableGrob(round(cm.svmBest,4)))
#cm.rfBest <- t(confusionMatrix(response, resultVector[[3]])$byClass)
#kable(round(cm.svmBest,4), caption = "Table 6: Classwise statistics of the Random Forest Classifier trained with 1000 trees and all Genes.")
#kable(round(cm.rfBest, 4), caption = "Table 6: Classwise performance of the Random Forest Classifier when trainied with 1000 trees and all genes")
#grid.arrange(tableGrob(round(cm.rfBest, 4)),
#             tableGrob(round(cm.svmBest,4)), textGrob("A"), textGrob("B"),
#            nrow = 2, ncol = 2)
## @knitr saving
#### Saving the image file
message("Saving Image file")
message(paste(memImageFile))
save.image(memImageFile)
