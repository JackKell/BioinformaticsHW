library(class)
library (e1071)
library(rpart)

# Read the data from a table and converts a column to the row names
readTableWithRowNames = function(filename, nameCol = 1) {
  # Reads the data
  data = read.delim(filename);
  # Set rownames
  dimnames(data)[[nameCol]] = data[,nameCol];
  # Remove rowname column from the data
  data[nameCol] = NULL;
  return(data);
}

pprint = function(...) {
  print(paste(...))
}

noneKNN = function(xData, yData) {
  kFolds = 10
  repeats = 5

  folds = cut(seq(1, nrow(xData)), breaks=kFolds, labels=FALSE);

  accuracies = c()

  for (r in 1:repeats) {
    # Feature Selection of none

    # Shuffle Data
    xData = xData[sample(nrow(xData)),];
    yData = yData[rownames(xData), , drop=FALSE];

    predictions = c()

    # 11 - KNN
    for (k in 1:kFolds) {
      testIndexes = which(folds == k, arr.ind=TRUE)
      xTestData = xData[testIndexes, ]
      xTrainData = xData[-testIndexes, ]
      yTrainData = yData[-testIndexes,,drop=FALSE]
      yTestData = yData[testIndexes,,drop=FALSE]

      prediction = as.numeric(knn(xTrainData, xTestData, unlist(yTrainData), k = 11)) - 1
      predictions = c(predictions, prediction)

    }
    accuracy = mean(predictions == unlist(yData))
    accuracies = c(accuracies, accuracy)
  }
  averageAccuracy = mean(accuracies)
  return(averageAccuracy)
}

correlation11KNN = function(xDataRaw, yDataRaw) {
  kFolds = 10
  repeats = 5

  folds = cut(seq(1, nrow(xDataRaw)), breaks=kFolds, labels=FALSE);
  yDataRaw = yDataRaw[rownames(xDataRaw), , drop=FALSE]

  accuracies = c()

  for (r in 1:repeats) {
    # Feature Selection based on correaltion
    correlationVector = abs(apply(xDataRaw, 2, function(x) {cor(unlist(x), unlist(yDataRaw))}))
    correlatedColumns = names(correlationVector[correlationVector >= quantile(correlationVector, prob=1-25/100)])

    # Shuffle Data
    xData = xDataRaw[sample(nrow(xDataRaw)), correlatedColumns, drop=FALSE];
    # print(xData)
    yData = yDataRaw[rownames(xData), , drop=FALSE];

    predictions = c()

    # 11 - KNN
    for (k in 1:kFolds) {
      testIndexes = which(folds == k, arr.ind=TRUE)
      xTestData = xData[testIndexes, , drop=FALSE]
      xTrainData = xData[-testIndexes, , drop=FALSE]
      yTrainData = yData[-testIndexes, , drop=FALSE]
      yTestData = yData[testIndexes, , drop=FALSE]

      results = as.numeric(knn(xTrainData, xTestData, as.numeric(unlist(yTrainData)), k = 11)) - 1

      predictions = c(predictions, results)
    }
    accuracy = mean(predictions == unlist(yData))
    accuracies = c(accuracies, accuracy)
  }
  averageAccuracy = mean(accuracies)
  return(averageAccuracy)
}

noneSvm = function(xDataRaw, yDataRaw) {
  kFolds = 10
  repeats = 5

  folds = cut(seq(1, nrow(xDataRaw)), breaks=kFolds, labels=FALSE);
  yDataRaw = yDataRaw[rownames(xDataRaw), , drop=FALSE]

  accuracies = c()

  for (r in 1:repeats) {
    # Shuffle Data
    xData = xDataRaw[sample(nrow(xDataRaw)), , drop=FALSE];
    # print(xData)
    yData = yDataRaw[rownames(xData), , drop=FALSE];

    predictions = c()

    for (k in 1:kFolds) {
      testIndexes = which(folds == k, arr.ind=TRUE)
      xTestData = xData[testIndexes, , drop=FALSE]
      xTrainData = xData[-testIndexes, , drop=FALSE]
      yTrainData = yData[-testIndexes, , drop=FALSE]
      yTestData = yData[testIndexes, , drop=FALSE]

      model = svm(xTrainData, yTrainData)
      prediction = round(as.numeric(predict(model, xTestData)))
      predictions = c(predictions, prediction)
    }

    accuracy = mean(predictions == unlist(yData))
    accuracies = c(accuracies, accuracy)
  }
  averageAccuracy = mean(accuracies)
  return(averageAccuracy)
}

correlatedSvm = function(xDataRaw, yDataRaw) {
  kFolds = 10
  repeats = 5

  folds = cut(seq(1, nrow(xDataRaw)), breaks=kFolds, labels=FALSE);
  yDataRaw = yDataRaw[rownames(xDataRaw), , drop=FALSE]

  accuracies = c()

  for (r in 1:repeats) {
    # Feature Selection based on correaltion
    correlationVector = abs(apply(xDataRaw, 2, function(x) {cor(unlist(x), unlist(yDataRaw))}))
    correlatedColumns = names(correlationVector[correlationVector >= quantile(correlationVector, prob=1-25/100)])

    # Shuffle Data
    xData = xDataRaw[sample(nrow(xDataRaw)), correlatedColumns, drop=FALSE];
    # print(xData)
    yData = yDataRaw[rownames(xData), , drop=FALSE];

    predictions = c()

    for (k in 1:kFolds) {
      testIndexes = which(folds == k, arr.ind=TRUE)
      xTestData = xData[testIndexes, , drop=FALSE]
      xTrainData = xData[-testIndexes, , drop=FALSE]
      yTrainData = yData[-testIndexes, , drop=FALSE]
      yTestData = yData[testIndexes, , drop=FALSE]

      model = svm(xTrainData, yTrainData)
      prediction = round(as.numeric(predict(model, xTestData)))
      predictions = c(predictions, prediction)
    }

    accuracy = mean(predictions == unlist(yData))
    accuracies = c(accuracies, accuracy)
  }
  averageAccuracy = mean(accuracies)
  return(averageAccuracy)
}

noneRPart = function(xDataRaw, yDataRaw) {
  kFolds = 10
  repeats = 5

  folds = cut(seq(1, nrow(xDataRaw)), breaks=kFolds, labels=FALSE);
  yDataRaw = yDataRaw[rownames(xDataRaw), , drop=FALSE]

  accuracies = c()

  for (r in 1:repeats) {
    # Shuffle Data
    xData = xDataRaw[sample(nrow(xDataRaw)), , drop=FALSE];
    # print(xData)
    yData = yDataRaw[rownames(xData), , drop=FALSE];

    predictions = c()

    for (k in 1:kFolds) {
      testIndexes = which(folds == k, arr.ind=TRUE)
      xTestData = xData[testIndexes, , drop=FALSE]
      xTrainData = xData[-testIndexes, , drop=FALSE]
      yTrainData = yData[-testIndexes, , drop=FALSE]
      yTestData = yData[testIndexes, , drop=FALSE]

      model = rpart(as.numeric(unlist(yTrainData))~., data=xTrainData, method="class")
      prediction = round(as.numeric(predict(model, xTestData)))
      predictions = c(predictions, prediction)
    }

    accuracy = mean(predictions == unlist(yData))
    accuracies = c(accuracies, accuracy)
  }
  averageAccuracy = mean(accuracies)
  print(averageAccuracy)
  return(averageAccuracy)
}

correlatedRPart = function(xDataRaw, yDataRaw) {
  kFolds = 10
  repeats = 5

  folds = cut(seq(1, nrow(xDataRaw)), breaks=kFolds, labels=FALSE);
  yDataRaw = yDataRaw[rownames(xDataRaw), , drop=FALSE]

  accuracies = c()

  for (r in 1:repeats) {
    # Feature Selection based on correaltion
    correlationVector = abs(apply(xDataRaw, 2, function(x) {cor(unlist(x), unlist(yDataRaw))}))
    correlatedColumns = names(correlationVector[correlationVector >= quantile(correlationVector, prob=1-25/100)])

    # Shuffle Data
    xData = xDataRaw[sample(nrow(xDataRaw)), correlatedColumns, drop=FALSE];
    # print(xData)
    yData = yDataRaw[rownames(xData), , drop=FALSE];

    predictions = c()

    for (k in 1:kFolds) {
      testIndexes = which(folds == k, arr.ind=TRUE)
      xTestData = xData[testIndexes, , drop=FALSE]
      xTrainData = xData[-testIndexes, , drop=FALSE]
      yTrainData = yData[-testIndexes, , drop=FALSE]
      yTestData = yData[testIndexes, , drop=FALSE]

      model = rpart(as.numeric(unlist(yTrainData))~., data=xTrainData, method="class")
      prediction = round(as.numeric(predict(model, xTestData)))
      predictions = c(predictions, prediction)
    }

    accuracy = mean(predictions == unlist(yData))
    accuracies = c(accuracies, accuracy)
  }
  averageAccuracy = mean(accuracies)
  print(averageAccuracy)
  return(averageAccuracy)
}

main = function() {
  # adrData = readTableWithRowNames("ADRs_HLGT_n438x232.txt")[1:20,]
  # geneData = readTableWithRowNames("gene_expression_n438x978.txt")[1:20, 1:4]

  adrData = readTableWithRowNames("ADRs_HLGT_n438x232.txt")
  geneData = readTableWithRowNames("gene_expression_n438x978.txt")

  # print(adrData[1:5, 1:3])
  # print(geneData[1:5, 1:3])

  methodNames = c("noneKNN", "correlation11KNN", "noneSvm", "correlatedSvm")
  methods = c(noneKNN, correlation11KNN, noneSvm, correlatedSvm)

  for (i in 1:3) {
    print(i)
    noneRPart(geneData, adrData[i])
    # for (method in methods) {
    #   print(method(geneData, adrData[i]))
    # }
  }
  return("Finished Without Problems")
}

main()
