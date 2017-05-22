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

methodExample = function(xData, yData) {
  kFolds = 10
  repeats = 5

  folds = cut(seq(1, nrow(xData)), breaks=kFolds, labels=FALSE);
  colnames = colnames(adrData)

  for (r in 1:repeats) {
    # Shuffle Data
    xData = xData[sample(nrow(geneData)),];
    yData = yData[rownames(geneData),][i];

    predictions = c();

    for (k in 1:kFolds) {
      testIndexes = which(folds == k, arr.ind=TRUE);
      xTestData = xData[testIndexes, ];
      xTrainData = xData[-testIndexes, ];
      adrTest = adrData[rownames(testData),][,i];
      adrTrain = adrData[rownames(trainData),][,i];

      nullModel = glm(adrTrain ~ 1, data = trainData, family = "binomial");
      fullModel = glm(adrTrain ~ ., data = trainData, family = "binomial");

      forwardModel = step(nullModel, scope = list(upper = fullModel), data = trainData, direction = "forward", trace = FALSE);
      backwardModel = step(fullModel, data = trainData, direction = "backward", trace = FALSE);

      fmPrediction = as.integer(round(predict(forwardModel, testData, type = "response")));
      bmPrediction = as.integer(round(predict(backwardModel, testData, type = "response")));

      currentFmPredictions = c(currentFmPredictions, fmPrediction);
      currentBmPredictions = c(currentBmPredictions, bmPrediction);
    }
    fmCm = table(adr[,1], currentFmPredictions);
    bmCm = table(adr[,1], currentBmPredictions);
  }
  # return the cross validation value for your predictions
}

methodOne = function(message) {
  # print(paste("Hey,", message))
  return(0)
}

methodTwo = function(message) {
  # print(paste("Hi,", message))
  return(1)
}

main = function() {
  set.seed(12345)
  adrData = readTableWithRowNames("ADRs_HLGT_n438x232.txt")
  geneData = readTableWithRowNames("gene_expression_n438x978.txt")

  print(adrData[1:5, 1:3])
  print(geneData[1:5, 1:3])

  methods = c(methodOne, methodTwo)
  for(i in 1:232) {
    for(method in methods) {
      method("Bob")
    }
  }
  return("Finished Without Problems")
}

main()
