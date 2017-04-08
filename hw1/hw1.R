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

questionOne = function(geneExpressionFile){
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

# Main funciton to answer all of the question posed in the R_hw1-1.pdf file
main = function() {
  # data directories
  geneExpressionFile = "gene_expression_n438x978.txt";
  adrHlgtFile = "ADRs_HLGT_n438c232.txt";

  questionOne(geneExpressionFile);

}


main();
