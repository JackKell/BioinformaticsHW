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

main = function() {
  data = readTableWithRowNames("gene_expression_n438x978.txt")
  print(data[1:4, 1:3])

  print(dim(data))
  kmeansResults = kmeans(data, 170)
  clusters = table(kmeansResults$cluster)
  print(clusters)

  output = NULL

  groupsSeq = seq(150, 220, 5)
  for (k in groupsSeq) {
    print(k)
    curr.kmeans <- kmeans(data, k)
    curr.sil <- silhouette(curr.kmeans$cluster, dist(data, method="euclidean"))
    output <- c(output, mean(curr.sil[, "sil_width"]))
  }

  png ("kmeans_sil.png")
  plot (groupsSeq, output, type="l", main="kmeans", xlab="number of clusters", ylab="average silhouette")
}

main()
