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

# Converts a correaltion matrix to a dataframe while removing redundant data
correlationToDataFrame = function(correlationData) {
  # Reference http://stackoverflow.com/questions/7074246/show-correlations-as-an-ordered-list-not-as-a-large-matrix
  # Flag redundant as NA
  correlationData[lower.tri(correlationData,diag=TRUE)] = NA;
  # Convert correlation table to a dataframe
  correlationData = as.data.frame(as.table(correlationData));
  # Remove all NA values
  correlationData = na.omit(correlationData);
  # Reset the number of the rows
  rownames(correlationData) = NULL;
  return(correlationData);
}

# Returns most similar drug as dataFrame
getMostSimilarDrug = function(correlationMatrix, drugName) {
  # Filters correlationMatrix for given drugName
  drugCorrelation = correlationMatrix[correlationMatrix[1] == drugName | correlationMatrix[2] == drugName, ];
  # Takes absolute value of correlations
  drugCorrelation[,3] = abs(drugCorrelation[,3]);
  # Sort by correlation
  drugCorrelation = drugCorrelation[order(drugCorrelation[,3], decreasing=TRUE),][1,];
  rownames(drugCorrelation) = NULL;
  return(drugCorrelation);
}

questionOne = function(geneExpressionFile) {
  drugGeneExpressionData = readTableWithRowNames(geneExpressionFile);
  drugGeneExpressionData = t(drugGeneExpressionData);
  drugCorrelation = correlationToDataFrame(cor(drugGeneExpressionData));

  png("drugCorrelationHistogram.png");
  hist(drugCorrelation[,3], xlab="Drug Correlation", main="Drug Correlation Histogram");

  medianDrugCorrelation = median(drugCorrelation[,3]);
  print(paste("Median Drug Correlation:", medianDrugCorrelation));

  # Get the top most positively correlated drugs
  topPositivelyCorrelatedDrugs = drugCorrelation[order(drugCorrelation$Freq, decreasing=TRUE),][1:10,];
  rownames(topPositivelyCorrelatedDrugs) = NULL;

  topPositivelyCorelatedDrug1 = topPositivelyCorrelatedDrugs[1,1];
  topPositivelyCorelatedDrug2 = topPositivelyCorrelatedDrugs[1,2];
  topPositivelyCorrelatedDrugCor = topPositivelyCorrelatedDrugs[1,3];
  print(paste("Top Positively Corelated Drug Pair: (", topPositivelyCorelatedDrug1, ",", topPositivelyCorelatedDrug2, ")", "with a correaltion of", topPositivelyCorrelatedDrugCor));

  print("Top 10 Most Positively Correlated Drugs");
  print(topPositivelyCorrelatedDrugs);

  png("topPositivelyCorrelatedDrugs.png");
  plot(
    drugGeneExpressionData[, topPositivelyCorelatedDrug1],
    drugGeneExpressionData[, topPositivelyCorelatedDrug2],
    xlab=topPositivelyCorelatedDrug1,
    ylab=topPositivelyCorelatedDrug2,
    main=paste(topPositivelyCorelatedDrug1, "over", topPositivelyCorelatedDrug2)
  );

  # Get the top most negatively correlated drugs
  topNegativelyCorrelatedDrugs = drugCorrelation[order(drugCorrelation$Freq, decreasing=FALSE),][1:10,];
  rownames(topNegativelyCorrelatedDrugs) = NULL;

  topNegativelyCorelatedDrug1 = topNegativelyCorrelatedDrugs[1,1];
  topNegativelyCorelatedDrug2 = topNegativelyCorrelatedDrugs[1,2];
  topNegativelyCorrelatedDrugCor = topNegativelyCorrelatedDrugs[1,3];

  print(paste("Top Negatively Corelated Drug Pair: (", topNegativelyCorelatedDrug1, ",", topNegativelyCorelatedDrug2, ")", "with a correlation of", topNegativelyCorrelatedDrugCor));

  png("topNegativelyCorrelatedDrugs.png");
  plot(
    drugGeneExpressionData[, topNegativelyCorelatedDrug1],
    drugGeneExpressionData[, topNegativelyCorelatedDrug2],
    xlab=topNegativelyCorelatedDrug1,
    ylab=topNegativelyCorelatedDrug2,
    main=paste(topNegativelyCorelatedDrug1, " over ", topNegativelyCorelatedDrug2)
  )

  print("Top 10 Most Negatively Correlated Drugs");
  print(topNegativelyCorrelatedDrugs);

  clofarabineCorr = getMostSimilarDrug(drugCorrelation, "CLOFARABINE");
  daunorubicinCorr = getMostSimilarDrug(drugCorrelation, "DAUNORUBICIN");
  fludarabineCorr = getMostSimilarDrug(drugCorrelation, "FLUDARABINE");

  print(clofarabineCorr);
  print(daunorubicinCorr);
  print(fludarabineCorr);
}

questionTwo = function(geneExpressionFile, adrHlgtFile) {
  geneData = readTableWithRowNames(geneExpressionFile)[,1:50];
  adrData = readTableWithRowNames(adrHlgtFile);

  fmSideEffectSE = NA;
  fmSE = .Machine$integer.max;
  fmModelSE = NA;
  bmSideEffectSE = NA;
  bmSE = .Machine$integer.max;
  bmModelSE = NA;

  fmSideEffectAIC = NA;
  fmSAIC = .Machine$integer.max;
  fmModelAIC = NA;
  bmSideEffectAIC = NA;
  bmSAIC = .Machine$integer.max;
  bmModelAIC = NA;
  warnings()
  colnames = colnames(adrData)
  for(i in 1:232) {
    sideEffectName = colnames[i];
    nullModel = glm(adrData[,i] ~ 1, data = geneData, family = "binomial");
    fullModel = glm(adrData[,i] ~ ., data = geneData, family = "binomial");
    forwardModel = step(nullModel, scope = list(upper = fullModel), data = geneData, direction = "forward", trace = FALSE);
    fmPrediction = predict(forwardModel, type = "response");
    fmCm = table(adrData[,i], round(fmPrediction));
    fmError = 0;
    if (ncol(fmCm) == 1) {
      fmError = fmCm[2,1];
    } else {
      fmError = fmCm[2,1] + fmCm[1,2];
    }
    fmAIC = AIC(forwardModel);

    backwardModel = step(fullModel, data = geneData, direction = "backward", trace = FALSE);
    bmPrediction = predict(backwardModel, type = "response");
    bmCm = table(adrData[,i], round(bmPrediction));
    bmError = 0;
    if (ncol(bmCm) == 1) {
      bmError = bmCm[2,1];
    } else {
      bmError = bmCm[2,1] + bmCm[1,2];
    }
    bmAIC = AIC(backwardModel);

    if(fmSAIC > fmAIC) {
      fmSideEffectAIC = sideEffectName;
      fmSAIC = fmAIC;
      fmModelAIC = forwardModel;
    }

    if(fmSE > fmError) {
      fmSideEffectSE = sideEffectName;
      fmSE = fmError;
      fmModelSE = forwardModel;
    }

    if(bmSAIC > bmAIC) {
      bmSideEffectAIC = sideEffectName;
      bmSAIC = bmAIC;
      bmModelAIC = backwardModel;
    }

    if(bmSE > bmError) {
      bmSideEffectSE = sideEffectName;
      bmSE = bmError;
      bmModelSE = backwardModel;
    }
  }

  print("Forward Model with lowest error:")
  print(fmSideEffectSE);
  print(paste("Error:", fmSE));
  print(coef(fmModelSE));
  print("backward Model with lowest error:")
  print(bmSideEffectSE);
  print(paste("Error:", bmSE));
  print(coef(bmModelSE));

  print("Forward Model with lowest AIC:")
  print(fmSideEffectAIC);
  print(paste("AIC:", fmSAIC));
  print(coef(fmModelAIC));
  print("Backward Model with lowest AIC:")
  print(bmSideEffectAIC);
  print(paste("AIC:", bmSAIC));
  print(coef(bmModelAIC));
}

questionThree = function(geneExpressionFile, adrHlgtFile) {
  geneData = readTableWithRowNames(geneExpressionFile)[,1:20];
  adrData = readTableWithRowNames(adrHlgtFile);
  kFolds = 10
  repeats = 3
  folds = cut(seq(1, nrow(geneData)), breaks=kFolds, labels=FALSE);


  SideEffectWithAverageMinimumFmError = NA;
  SideEffectWithAverageMinimumFmErrorValue = .Machine$integer.max;

  SideEffectWithAverageMinimumBmError = NA;
  SideEffectWithAverageMinimumBmErrorValue = .Machine$integer.max;

  SideEffectWithAverageMinimumFmAIC = NA;
  SideEffectWithAverageMinimumFmAICValue = .Machine$integer.max;

  SideEffectWithAverageMinimumBmAIC = NA;
  SideEffectWithAverageMinimumBmAICValue = .Machine$integer.max;

  colnames = colnames(adrData)
  for (i in 1:ncol(adrData)) {
    sideEffectName = colnames[i];
    fmAverageError = 0;
    bmAverageError = 0;
    fmAverageAIC = 0;
    bmAverageAIC = 0;

    for (r in 1:repeats) {
      print(paste(i, ":", r))
      # Shuffle Data
      geneData = geneData[sample(nrow(geneData)),];
      # print(geneData)
      temp = adrData[rownames(geneData),][i];
      # print(adrData)

      currentFmPredictions = c();
      currentBmPredictions = c();

      for (k in 1:kFolds) {
        testIndexes = which(folds == k, arr.ind=TRUE);
        testData = geneData[testIndexes, ];
        # print(testData)
        trainData = geneData[-testIndexes, ];
        # print(trainData)
        adrTest = adrData[rownames(testData),][,i];
        # print(adrTest)
        adrTrain = adrData[rownames(trainData),][,i];
        # print(adrTrain)

        nullModel = glm(adrTrain ~ 1, data = trainData, family = "binomial");
        fullModel = glm(adrTrain ~ ., data = trainData, family = "binomial");

        forwardModel = step(nullModel, scope = list(upper = fullModel), data = trainData, direction = "forward", trace = FALSE);
        backwardModel = step(fullModel, data = trainData, direction = "backward", trace = FALSE);

        fmPrediction = as.integer(round(predict(forwardModel, testData, type = "response")));
        # print(fmPrediction)
        bmPrediction = as.integer(round(predict(backwardModel, testData, type = "response")));
        # print(bmPrediction)

        currentFmPredictions = c(currentFmPredictions, fmPrediction);
        currentBmPredictions = c(currentBmPredictions, bmPrediction);
        # print(currentFmPredictions)
        # print(currentBmPredictions)

        # fmAverageAIC = fmAverageAIC + AIC(forwardModel);
        # bmAverageAIC = bmAverageAIC + AIC(backwardModel);
      }
      fmCm = table(temp[,1], currentFmPredictions);
      bmCm = table(temp[,1], currentBmPredictions);

      # print(fmCm)
      # print(bmCm)

      if (nrow(fmCm) == 2) {
        fmAverageError = fmAverageError + fmCm[2,1];
      }

      if (ncol(fmCm) == 2) {
        fmAverageError = fmAverageError + fmCm[1,2];
      }

      if (nrow(bmCm) == 2) {
        bmAverageError = bmAverageError + bmCm[2,1];
      }

      if (ncol(bmCm) == 2) {
        bmAverageError = bmAverageError + bmCm[1,2];
      }

      # print(fmAverageError)
      # print(bmAverageError)
    }

    fmAverageError = fmAverageError / repeats;
    bmAverageError = bmAverageError / repeats;
    # print(fmAverageError)
    # print(bmAverageError)
    fmAverageAIC = fmAverageAIC / repeats;
    bmAverageAIC = bmAverageAIC / repeats;

    if(SideEffectWithAverageMinimumFmAICValue > fmAverageAIC) {
      SideEffectWithAverageMinimumFmAIC = sideEffectName;
      SideEffectWithAverageMinimumFmAICValue = fmAverageAIC;
    }

    if(SideEffectWithAverageMinimumFmErrorValue > fmAverageError) {
      SideEffectWithAverageMinimumFmError = sideEffectName;
      SideEffectWithAverageMinimumFmErrorValue = fmAverageError;
    }

    if(SideEffectWithAverageMinimumBmAICValue > bmAverageAIC) {
      SideEffectWithAverageMinimumBmAIC = sideEffectName;
      SideEffectWithAverageMinimumBmAICValue = bmAverageAIC;
    }

    if(SideEffectWithAverageMinimumBmErrorValue > bmAverageError) {
      SideEffectWithAverageMinimumBmError = sideEffectName;
      SideEffectWithAverageMinimumBmErrorValue = bmAverageError;
    }
  }
  print("Forward Model with lowest error:")
  print(SideEffectWithAverageMinimumFmError);
  print(paste("Average Error:", SideEffectWithAverageMinimumFmErrorValue));

  print("Backward Model with lowest error:")
  print(SideEffectWithAverageMinimumBmError);
  print(paste("Average Error:", SideEffectWithAverageMinimumBmErrorValue));

  # print("Forward Model with lowest AIC:")
  # print(SideEffectWithAverageMinimumFmAIC);
  # print(paste("AIC:", SideEffectWithAverageMinimumFmAICValue));
  #
  # print("Backward Model with lowest AIC:")
  # print(SideEffectWithAverageMinimumBmAIC);
  # print(paste("Average AIC:", SideEffectWithAverageMinimumBmAICValue));
}

confunctionMatrix = function(x, y) {
  both = union(x, y);
  inX = both %in% x;
  inY = both %in% y;
  return(table(inX, inY))
}

# Main funciton to answer all of the question posed in the R_hw1-1.pdf file
main = function() {
  # data directories
  geneExpressionFile = "gene_expression_n438x978.txt";
  adrHlgtFile = "ADRs_HLGT_n438x232.txt";

  # questionOne(geneExpressionFile);
  # questionTwo(geneExpressionFile, adrHlgtFile);
  questionThree(geneExpressionFile, adrHlgtFile);

}


main();
