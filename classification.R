## @knitr loading
## done in class, loads all the input data
load("DenBoerData_loaded.Rdata")

## register a paralell backend and define number of Threads
## Note the the numberThreads variable gets passed later to the ranger function
## this is neccessary, because ranger uses all available threads as a default and I
## did not intend to overload the cluster.
numberThreads <- 4
doMC::registerDoMC(numberThreads)
## Here I store all the different parameters we want to use for the different
## classification methods
numberOfTrees <- c(200, 500, 1000)
kernels <- c("radial", "linear", "polynomial", "sigmoid")
## for reusability provide the name of the class variable for this dataset
## I also save the minimum number allowed cases per class in an extra variable
classVariable <- "sample.labels"
minNumberofCasesPerClass <- 30
preFiltered <- 200
parameters <- c("RF 200 Trees all Genes",
                "RF 500 Trees all Genes",
                "RF 1000 Trees all Genes",
                "RF 200 Trees Feature Selection",
                "RF 500 Trees Feature Selection",
                "RF 1000 Trees Feature Selection",
                "SVM Radial Kernel all Genes",
                "SVM Linear Kernel all Genes",
                "SVM Ploynomial Kernel all Genes",
                "SVM Sigmoid Kernel all Genes",
                "SVM Radial Kernel Feature Selection",
                "SVM Linear Kernel Feature Selection",
                "SVM Ploynomial Kernel Feature Selection",
                "SVM Sigmoid Kernel Feature Selection")


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


## @knitr filter
## first we filter the classes and include only classes with more than 30 patients
# Function for class filtering
FilterClass <- function(minimumNumberofCases){
  releventClasses <- names(which(table(
    phenoTable[, classVariable]) > minimumNumberofCases))
  helper <- phenoTable[, classVariable] %in% releventClasses
  return(helper)
}

relevantPatients <- FilterClass(minNumberofCasesPerClass)
phenoTable <- phenoTable[relevantPatients, ]
#remove unused factors
phenoTable[,classVariable] <- factor(phenoTable[,classVariable])
exprTableTransposed <- exprTableTransposed[relevantPatients, ]
classesAfterFilerd <- table(phenoTable[, classVariable])
grid.arrange(tableGrob(as.data.frame(classesAfterFilerd)))


## @knitr prepareSelect
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
getNumberofFeatures <- function() {
  selectedGenes <- SelectFeaturesWithRandomForest(c(1:nrow(exprTableTransposed)),
                                                    ncol(exprTableTransposed))

  lo <- loess(selectedGenes[1:preFiltered] ~ c(1:200))
  smoothed = predict(lo)
  secondDer <- diff(diff(smoothed))
  maximalChangePoint <- max(secondDer)
  maximalChangeIndex <- match(maximalChangePoint, secondDer)
  numberFeatures <- maximalChangeIndex
  return(list(numberFeatures = numberFeatures,
           smoothedLine = smoothed,
           genes = selectedGenes))
}


helper <- getNumberofFeatures()
numberFeatures <- helper[["numberFeatures"]]

plotVariableImportance <- function(smoothed, selectedGenes) {
  plotData <- data.frame(`Variable Importance` = selectedGenes[1:preFiltered],
                         Index = c(1:preFiltered),
                         Smoothed = out,
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


plotVariableImportance(helper[["smoothedLine"]], helper[["genes"]])



## @knitr select
#### Feature Selection ####
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



## Selection with leave one out cross validation
SelectRandomForestLOOCV <- function(numFeatures,
                                    verbose = FALSE) {
  res <- foreach(i = 1:nrow(exprTableTransposed)) %dopar% {
    curVariables <- SelectFeaturesWithRandomForest(-i,
                                                   numFeatures,
                                                   verbose)
    return(curVariables)

  }
  names(res) <- rownames(phenoTable)
  return(res)
}


#### Classifier Functions ####
## ----------------------------------------------------------------------------
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
## ----------------------------------------------------------------------------
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



#### Functions for displaying the results ####
## ----------------------------------------------------------------------------
getResultOverview <- function (results) {
  evaluationResults <- foreach(i = c(1:length(results)), .combine = 'rbind') %do% {
    res <- results[[i]]
    cm <- confusionMatrix(res, response)
    check <- names(response) == names(res)
    message("Patients in Same Order:")
    print(all(check,TRUE))
    cm.table <- cm$table
    print(cm.table)
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
                      CI = ci.interval))
  }
  rownames(evaluationResults) <- names(results)
  return(evaluationResults)
}



## Functions for plotting a heatmap of a confusion table
plotResults <- function(res, response, title) {
  cm <- confusionMatrix(res, response)
  cm.table <- as.data.frame(cm$table)
  cm.stats <-data.frame(Statistics = cm$overall)
  cm.stats$Statistics <- round(cm.stats$Statistics,2)
  cm.percentage <- as.data.frame(prop.table(cm$table))
  cm.table$Perc <- round(cm.percentage$Freq*100,2)

  cm.plot <- ggplot(data = cm.table, aes(x = Prediction , y =  Reference, fill = Freq)) +
             geom_raster(aes(fill = Freq)) +
             geom_text(aes(label = paste("",Freq,",",Perc,"%")),
             color = 'white', size = 3, fontface = "bold") +
             theme_light()

  cm.statsTable <-  tableGrob(cm.stats)

  grid.arrange(cm.plot, cm.statsTable,nrow = 1, ncol = 2,
  top = textGrob(paste0("Confusion Table Heatmap \n",title), gp = gpar(fontsize=20,font=1)))

}




#### EXECUTION ####
## ----------------------------------------------------------------------------

## @knitr execution
## finally we can perform the acutal selection
loocvSelections <- SelectRandomForestLOOCV(numberFeatures)
write.csv(loocvSelections, 'loocvSelections.csv')


allGenesSelected <- rep(list(selectedGenes), nrow(exprTableTransposed))
names(allGenesSelected) <- rownames(exprTableTransposed)
selections <- list(allGenes = allGenesSelected,
                   rfSelection = loocvSelections)
tunedMtry = tunedMtryAllFeatures
resultVector <- foreach(selection = names(selections), .combine = "c") %do% {
  if(selection == "rfSelection"){
    tunedMtry = tunedMtryFeatureSelection
  }
  message(tunedMtry)
  selData <- selections[[selection]]
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


## @knitr evaluation
## Reordering the results to have the same methods together
index <- c(1,2,3,8,9,10,4,5,6,7,11,12,13,14)
resultVector <- resultVector[order(index)]

names(resultVector) <- parameters
response <- phenoTable[, classVariable]
response <- factor(response, levels = rev(levels(response)))
names(response) <- rownames(phenoTable)
resultTable <- getResultOverview(resultVector)

plotResults(resultVector[[1]], response, names(resultVector)[1])
plotResults(resultVector[[2]], response, names(resultVector)[2])
plotResults(resultVector[[3]], response, names(resultVector)[3])
plotResults(resultVector[[4]], response, names(resultVector)[4])
plotResults(resultVector[[5]], response, names(resultVector)[5])
plotResults(resultVector[[6]], response, names(resultVector)[6])
plotResults(resultVector[[7]], response, names(resultVector)[7])
plotResults(resultVector[[8]], response, names(resultVector)[8])
plotResults(resultVector[[9]], response, names(resultVector)[9])
plotResults(resultVector[[10]], response, names(resultVector)[10])
plotResults(resultVector[[11]], response, names(resultVector)[11])
plotResults(resultVector[[12]], response, names(resultVector)[12])
plotResults(resultVector[[13]], response, names(resultVector)[13])
plotResults(resultVector[[14]], response, names(resultVector)[14])


#### Saving the image file
## ----------------------------------------------------------------------------
message("Saving Image file")
message(paste(memImageFile))
save.image(memImageFile)
