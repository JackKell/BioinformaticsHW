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
  # Test all of the possible value of k from 10 to 435 incrementing by 10 to find a rough range of where the best cluster size is
  stepSize = 20
  roughClusterSizes = seq(stepSize, nrow(data), stepSize)
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
  roughBestClusterSize = names(output)[which(output==max(output))]
  print(roughBestClusterSize)
  fineOutput = NULL
  fineClusterSizes = seq(roughBestClusterSize - stepSize, roughBestClusterSize + stepSize)
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



  kmeansResults = kmeans(data, 170)
  clusters = table(kmeansResults$cluster)
  print(clusters)
  print(sum(clusters))


}

main = function() {
  set.seed(12345)
  questionOne()
}

main()
