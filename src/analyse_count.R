args <- commandArgs(trailingOnly=TRUE)
path <- args[[1]]
count <- read.csv(path, stringsAsFactors = F, comment.char = "")

# total number of clusters
cat("\nTotal number of consensus spectra:\n")
cat(length(unique(count$msclusterID)))
cat("\n")

# total number of not blank clusters
if ('BLANKS_TOTAL' %in% names(count)) {
  cat("\nTotal number of consensus spectra NOT blank:\n")
  cat(length(unique(count$msclusterID[count$BLANKS_TOTAL == 0])))
  cat("\n")
}

# qtd of clusters by size
cat("\nNumber of consensus spectra by size (number of clustered spectra):\n")
prettyNum(table(count$numSpectra))

# total number of spectra
cat("\nTotal number of processed spectra from the raw MS/MS data:\n")
cat(sum(count$numSpectra))

# number of unique masses
cat("\nNumber of unique masses rounded in digit 2:\n")
cat(length(unique(round(count$mzConsensus, 2))))

# number of clusters by masses
cat("\nMean number of consensus spectra by m/z: \n")
cat(round(mean(table(round(count$mzConsensus, 2))), 3))
