library(knitr)
library(ranger)
library(caret)
library(foreach)
library(doMC)
library(e1071)


## register a paralell backend and define number of Threads
numberThreads <- 4
registerDoMC(numberThreads)

####  Define Local directories and files ####
## Find the home directory
myHome <- Sys.getenv("HOME")

## Define the main directory for the  data and results
mainDir <- file.path(myHome, "uni/4_sem/biostatistics/report/")
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


#### Define variables for the location and name of the input files ####
## URL fo the folder containing the data fle
## Maybe refactor and pass via command line arguments later
dataURL <- "https://github.com/jvanheld/stat1/raw/master/data/DenBoer_2009/"

files <- c(
  expr = "GSE13425_Norm_Whole.tsv.gz",
  pheno = "phenoData_GSE13425.tsv.gz",
  groups = "GSE13425_group_descriptions.tsv.gz"
)



## Function that downloads input files in case they don't exist locally yet
DownloadMyFiles <- function(files,
                            dataURL,
                            destDir) {

  ## Create local directory
  dir.create(destDir, recursive = TRUE, showWarnings = FALSE)

  for (f in files) {

    ## Destination file
    destFile <- file.path(destDir, f)

    ## Check if file exists
    if (file.exists(destFile)) {
      message("skipping download because file exists: \n", destFile)
    } else {
      sourceURL <- file.path(dataURL, f)
      download.file(url = sourceURL, destfile = destFile)
    }
  }
}




## Call the function to download the files
DownloadMyFiles(files = files, dataURL = dataURL, destDir = destDir)

kable(data.frame(list.files(destDir)),
      caption = "Content of the destination directory after file download. ")



## ---------------------------------------------------------------------
#### Load data tables ####

## Load expression table
exprTable <- read.table(file.path(destDir, files["expr"]),
                        sep = "\t",
                        header = TRUE,
                        quote = "",
                        row.names = 1)
dim(exprTable)
names(exprTable)
head(exprTable)

## Load metadescriptions (pheno table)
phenoTable <- read.table(file.path(destDir, files["pheno"]),
                         sep = "\t",
                         header = TRUE,
                         quote = "",
                         row.names = 1
                         )
dim(phenoTable)
names(phenoTable)
head(phenoTable)
## for reusability provide the name of the class variable for this dataset
classVariable <- "sample.labels"

## Load group descriptions
groupDescriptions <- read.table(file.path(destDir, files["groups"]),
                                sep = "\t",
                                header = TRUE,
                                quote = "",
                                row.names = 1)
dim(groupDescriptions)
print(groupDescriptions)
kable(groupDescriptions)



#### Prepare data for classification task ####
## check if pheno data and expression data have same order of patients
check <- names(exprTable) == rownames(phenoTable)
message("Do Pheno and Expressiondata have the same order of patients")
print(all(check = T))
exprTableTransposed <- t(exprTable)

## this prevents an error in the ranger function due to illliegal names
colnames(exprTableTransposed) <- make.names(colnames(exprTableTransposed))



#### Feature Selection ####
## Function that selects genes based on the random forest variable importance

SelectFeaturesWithRandomForest <- function(featureNumber) {
    fit <- ranger(y = phenoTable[, classVariable],
                  x = exprTableTransposed[,1:50],
                  importance = "impurity")
    varImportance <- fit$variable.importance
    print(varImportance)
    selectedGenes <- sort(varImportance, na.last = TRUE, decreasing = TRUE)
    plot(sort(varImportance, na.last = TRUE, decreasing = TRUE))
    selectedGenes <- selectedGenes[1:featureNumber]
    return(selectedGenes)
}



## Function that classifies patients with random forest
RandomForestClassifier <- function(numberTrees,
                                   tunedMtry,
                                   trainIndex,
                                   testIndex,
                                   selectedCovariates) {

    fit <- ranger(y = phenoTable[trainIndex, classVariable],
                     x = exprTableTransposed[trainIndex, selectedCovariates],
                     num.trees = numberTrees,
                     num.threads = numberThreads,
                     mtry = tunedMtry)

    predicted <- predict(fit,
                         exprTableTransposed[testIndex, selectedCovariates])
    return(predicted$predictions)
}



#### Function for LOOCV ####
## --------------------------------------------------------------------

RandomForestLOOCV <- function(numberTrees,
                              my.mtry,
                              selectedCovariates) {
    inputData <- exprTableTransposed[,selectedCovariates]
    res <- foreach(i = 1:nrow(inputData)) %dopar% {
        curTrain <- inputData[-i,]
        curTest <- t(as.data.frame(inputData[i,]))
        curTrainGroup <- phenoTable[-i,classVariable]

        rf.out <- ranger(y = curTrainGroup, x = curTrain)
        predicted <- predict(rf.out, curTest)
        res <- predicted$predictions
        res <- droplevels(res)
        return(res)
    }
    names(res) <- rownames(phenoTable)
    return(unlist(res))
}



#### SVM sandbox ####

svm1 <- svm(y = phenoTable[, classVariable], x = exprTableTransposed, 
            method="C-classification", kernal="radial", 
            gamma=0.1, cost=10)
plot(svm1)
pred <- predict(svm1, exprTableTransposed)
plot(svm1, data = data.frame(exprTableTransposed))
xtab <- table(phenoTable[, classVariable], pred)
xtab
table(phenoTable[,classVariable])
#### EXECUTION ####t
## --------------------------------------------------------------------

pr <- RandomForestClassifier(200, 2, c(1:190), c(1:190),
                             colnames(exprTableTransposed)[1:8])
pr
train.control <- trainControl(method = "LOOCV")
rf.fit <- train(method = "ranger", x = testing, y = classVariableVector,
                trainControl = ctrain.control)

prLoocv <- RandomForestLOOCV(50, 2, colnames(exprTableTransposed)[1:5])
prLoocv
selected <- SelectFeaturesWithRandomForest(10);
save(memImageFile);
