library(cluster)

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

questionOne = function() {
  # Read in the data
  data = t(readTableWithRowNames("gene_expression_n438x978.txt"))
  print(dim(data))
  print(data[1:4, 1:3])

  roughOutput = NULL
  # Test all of the possible values of k incrementing by a rough range to find where the best cluster size is in general
  stepSize = 20
  roughClusterSizes = seq(2, nrow(data), stepSize)
  # roughClusterSizes = seq(2, 40, stepSize)
  # The number of trails done to get the averagered average silhouette value for the current groups size (k)
  trails = 3.0
  for (k in roughClusterSizes) {
    averageSilhouetteOnAverage = 0
    for (i in 1:trails) {
      print(paste(k, i))
      currentKmeansResults = kmeans(data, k)
      currentSilhouetteValues = silhouette(currentKmeansResults$cluster, dist(data, method="euclidean"))
      averageSilhouetteValue = mean(currentSilhouetteValues[, "sil_width"])
      averageSilhouetteOnAverage = averageSilhouetteOnAverage + averageSilhouetteValue
    }
    averageSilhouetteOnAverage = averageSilhouetteOnAverage / trails
    roughOutput = c(roughOutput, averageSilhouetteOnAverage)
  }

  png ("kmeans_sil.png")
  plot (roughClusterSizes, roughOutput, type="l", main="The Average Silhouetee Per Cluster for a\nGiven Number of Clusters on Average", xlab="Number of Clusters", ylab="Average Silhouette Value Per Cluster On Average")
  names(roughOutput) = roughClusterSizes
  print(roughOutput)

  # From the results in the output we find that the max is some where around
  roughBestClusterSize = names(roughOutput)[which(roughOutput==max(roughOutput))]

  # print(roughBestClusterSize)
  fineOutput = NULL
  fineClusterSizes = seq(2, 20)
  for (k in fineClusterSizes) {
    averageSilhouetteOnAverage = 0
    for (i in 1:trails) {
      print(paste(k, i))
      currentKmeansResults = kmeans(data, k)
      currentSilhouetteValues = silhouette(currentKmeansResults$cluster, dist(data, method="euclidean"))
      averageSilhouetteValue = mean(currentSilhouetteValues[, "sil_width"])
      averageSilhouetteOnAverage = averageSilhouetteOnAverage + averageSilhouetteValue
    }
    averageSilhouetteOnAverage = averageSilhouetteOnAverage / trails
    fineOutput = c(fineOutput, averageSilhouetteOnAverage)
  }

  png ("kmeans_sil_fine.png")
  plot (fineClusterSizes, fineOutput, type="l", main="The Average Silhouetee Per Cluster for a\nGiven Number of Clusters on Average", xlab="Number of Clusters", ylab="Average Silhouette Value Per Cluster On Average")

  names(fineOutput) = fineClusterSizes
  print(fineOutput)

  kmeansResults = kmeans(data, 2)
  clusters = table(kmeansResults$cluster)
  png("histtest.png")
  plot(clusters)
  print(clusters)

  sil = silhouette(kmeansResults$cluster, dist(data, method="euclidean"));
  png("clusterVis.png");
  plot(sil, main="Silhouette Plot");

}

questionTwo = function() {
  # Read in the data
  data = t(readTableWithRowNames("gene_expression_n438x978.txt"))
  clusterSize = 2
  kmeansResults = kmeans(data, 2)

  # Calculate the size of the geneset universe
  con = file("gene_lists.txt",  open="r")
  count = 1
  universe = unlist(row.names(data))
  while (TRUE) {
    line = readLines(con, n=1)
    if (length(line) == 0) {
      break
    }
    geneSet = unlist(strsplit(line, "\t"))
    universe = union(geneSet, universe)
  }
  universeSize = length(universe)
  print(universeSize)
  close(con)

  # Do bio calculations
  con = file("gene_lists.txt",  open="r")
  output = NULL
  for (i in 1:5000) {
    print(i)
    line = readLines(con, n=1)
    currentGenes = unlist(strsplit(line, "\t"))
    if (length(line) == 0) {
      break
    }
    for (clusterIndex in 1:2) {
      currentCluster = row.names(data[kmeansResults$cluster==clusterIndex,])
      tp = length(intersect(currentGenes, currentCluster))
      fp = length(currentCluster) - tp
      fn = length(currentGenes) - tp
      tn = universeSize - tp - fp - fn
      confusionMatrix = matrix(c(tp, fn, fp, tn), nr=2, dimnames=list(c("DE", "notDE"), c("InGeneSet", "NotInGeneSet")))
      pvalue = fisher.test(confusionMatrix, alternative="greater")$p.value
      newOutputRow = c(clusterIndex, i, tp, fp, fn, tn, pvalue)
      output = rbind(output, newOutputRow)
    }
    count = count + 1
  }
  close(con)
  output = data.frame(output)
  colnames(output) = c("cluster", "genesetrow", "tp", "fp", "fn", "tn", "pvalue")
  # print(output)
  print(output[which(output["pvalue"]==max(output["pvalue"]))[1],])
  print(output[which(output["pvalue"]==min(output["pvalue"]))[1],])

}

main = function() {
  set.seed(12345)
  # questionOne()
  questionTwo()
}

main()
