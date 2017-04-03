main = function() {
  geneExpressionFile = "gene_expression_n438x978.txt";
  adrHlgtFile = "ADRs_HLGT_n438c232.txt";

  # Format geneData
  geneData = read.delim(geneExpressionFile);
  dimnames(geneData)[[1]] = geneData[,1];
  geneData = geneData[,-1];

  # Reference http://stackoverflow.com/questions/7074246/show-correlations-as-an-ordered-list-not-as-a-large-matrix
  # Get the correlation of the genes
  geneCorrelation = cor(geneData);
  # Mark the values in the bottom half and the diagnal of the correlation as NA because those values are either redundant or meaningless.
  geneCorrelation[lower.tri(geneCorrelation,diag=TRUE)]=NA;
  # Convert the matrix to a dataframe
  geneCorrelation=as.data.frame(as.table(geneCorrelation));
  # Remove any rows that have a frequency of NA
  geneCorrelation=na.omit(geneCorrelation);
  # Get the top most positively correlated genes
  topPositivelyCorrelatedGenes=geneCorrelation[order(geneCorrelation$Freq, decreasing=TRUE),][1:10,]
  rownames(topPositivelyCorrelatedGenes) = NULL;
  topNegativelyCorrelatedGenes=geneCorrelation[order(geneCorrelation$Freq, decreasing=FALSE),][1:10,]
  rownames(topNegativelyCorrelatedGenes) = NULL;
  print(topPositivelyCorrelatedGenes);
  print(topNegativelyCorrelatedGenes);
}

main();
