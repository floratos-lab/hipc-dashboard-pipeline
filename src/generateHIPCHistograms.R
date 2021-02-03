# BiocManager::install("ggplot2")
# BiocManager::install("gridExtra")

# Plot histograms of response agents per signature

library(ggplot2)
library(gridExtra)

inputFileBase <- "hipc_2020-08-11"
versionNum    <- "30" # download processsing iteration
version       <- paste0("v", versionNum) 
outFile       <- "hipc-uniqueCountBySignature.png"

geneSheetNameOut <- "gene_expression"

cellSheetNameOut <- "cell_type_frequency"
 
geneFilename <- paste(inputFileBase, geneSheetNameOut, version, sep = "_")
geneFilename <- paste(geneFilename, "responseAgentCountsByRow.csv", sep = "-")

cellFilename <- paste(inputFileBase, cellSheetNameOut, version, sep = "_")
cellFilename <- paste(cellFilename, "responseAgentCountsByRow.csv", sep = "-")

genesCntDF <- read.csv(geneFilename)
cellsCntDF <- read.csv(cellFilename)

png(outFile)

p1 <- ggplot(genesCntDF, aes(count)) + geom_histogram(bins = 10) + scale_x_log10() +
  labs(x = "Number of genes") + labs(y = "Count") + 
  labs(title = "Genes per signature", tag = "A")  +
  annotation_logticks(sides = "b") # +
  # theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))

p2 <- ggplot(cellsCntDF, aes(count)) + geom_histogram(bins = 9) + scale_x_continuous(breaks = seq(1, 9, 1)) +
  labs(x = "Number of cell types") + labs(y = "Count") + 
  labs(title = "Cell types per signature", tag = "B") # +
  # theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))

grid.arrange(p1, p2, nrow = 1)

dev.off()
