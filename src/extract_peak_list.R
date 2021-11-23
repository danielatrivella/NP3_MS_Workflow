## ----load-libs, message = FALSE--------------------------------------------
cat("Loading Rcpp functions...\n")
library(purrr)
script_path <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  } else {
    # 'source'd via R console
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
}
Rcpp::sourceCpp(file.path(script_path(), "read_mgf_peak_list_R.cpp"))

bin_size <- 0.05
scale_factor <- 0.5
trim_mz = -1 # always return full list of peaks

inverse_scale <- function(ints, scale_factor)
{
  if (scale_factor == 0) # log scale was applied
  {
    # inverse of 1.0 + log(1.0 + multVal * ints[i]))
    ints <- expm1(ints-1)
  } else {
    ints <- ints^(1/scale_factor)
  }
  return(ints)
}

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Six arguments must be supplied to create the pairwise similarity table:\n", 
       " 1 - Job name;\n",
       " 2 - Path to the MGF file with the spectra to be compared pairwise;\n",
       " 3 - Path to the result count table;\n", 
       " 4 - The bin size to consider two fragmented peaks m/z's the same;\n",
       " 5 - The spectra fragmented peaks scaling method: sqrt, ln or NULL;\n",
       call.=FALSE)
} else {
  data_name <- args[[1]]
  path_mgf_dir <- file.path(args[[2]])
  counts_path <- file.path(args[[3]]) 
  if (!file.exists(counts_path))
  {
    stop("The count file '", counts_path,
         "' do not exists. Provide a valid path to the counts file.")
  }
  bin_size <- as.numeric(args[[4]])
  if (is.na(bin_size) || bin_size < 0) {
    stop("Wrong bin size value '", args[[4]],
         "'. It must be a numeric value greater or equal to zero.")
  }
  scale_factor <- as.numeric(args[[5]])
  if (!is.na(scale_factor) || !is.null(scale_factor))
  {
    scale_factor <- as.numeric(scale_factor[[1]])
    if (scale_factor < 0)
    {
      warning("Invalid scale factor provided, it must be a numeric value greater or equal to 0: \n", 
              "  - 0 : ln natural logarithm of the intensities\n",
              "  - 1 : no scale\n",
              "  - x : x > 0 pow of the intensities to x. (e.g. x = 0.5 square root)\n",
              "Scale factor equals 0.5 (sqrt) will be selected by default.", call. = FALSE)
      scale_factor <- 0.5
    }
  } else {
    scale_factor <- 1 # no scale
  }
}

if (!dir.exists(path_mgf_dir) || !any(grepl(".mgf$", list.files(path_mgf_dir))))
{
  stop("The MGF dir '", path_mgf_dir,
       "' do not exists or do not contains MGF files. ", 
       "Provide a valid path to where the MGF files are located.")
}
# get the path to each mgf file inside the mgf dir
path_mgf <- file.path(path_mgf_dir, list.files(path_mgf_dir)[grepl("_[0-9]+.mgf$", list.files(path_mgf_dir))])
# order paths by mass
path_mgf <- path_mgf[order(as.numeric(gsub(paste0(data_name, "_[0-9]+_[0-9]+_|_[0-9]+.mgf"), "", basename(path_mgf))))]
n_mgf <- length(path_mgf)
if (n_mgf == 0) {
  stop("The MGF dir '", path_mgf_dir,
       "' do not contains MGF files from ", data_name,
       ". Provide a valid path to where the MGF files are located.")
}

ms_count <- read.csv(counts_path, stringsAsFactors = FALSE, comment.char = "", strip.white = TRUE)
ms_count$peaksInt <- ms_count$peaksList <- NA

cat("\n  * Extracting, normalizing and scaling fragment peaks list of", data_name, 
        "clustered spectra to concatenate in the counts files *\n")
# add progress
progress_comp <- seq(from = 1, to = n_mgf, by = 1)
cat("        |", rep("", length(progress_comp)-1), "|\n        |")
ti <- Sys.time()
for (i in seq_along(path_mgf)) {
  cat("=")
  # read the ms2 data and get the scan indexes
  ms2_sample <- readMgfPeaksList(path_mgf[[i]], bin_size = bin_size, 
                                 trim_mz = trim_mz, scale_factor = scale_factor,
                                 join_isotopic_peaks = 0)
  
  ms2_sample$MZS <- lapply(ms2_sample$MZS, function(x) paste0(round(x, 4), collapse = ";"))
  ms2_sample$INTS <- lapply(ms2_sample$INTS, function(x) paste0(round(x, 4), collapse = ";"))
  
  match_scans_order <- match(ms2_sample$SCANS, ms_count$msclusterID)
  ms_count$peaksList[match_scans_order] <- unlist(ms2_sample$MZS, recursive = FALSE)
  ms_count$peaksInt[match_scans_order] <- unlist(ms2_sample$INTS, recursive = FALSE)
}

cat("|\n")
tf <- Sys.time()
cat("    * Done extracting peak list in", round(tf-ti, 2), units(tf-ti), "*\n")

if (any(is.na(ms_count$peaksList)) || any(is.na(ms_count$peaksInt)))
  stop("Bad extraction of peak lists. Some scans are missing, check if the paths are OK and retry.")

# check if the inverse scaled intensities sum 1000
check_ints_sum <- sapply(ms_count$peaksInt, function(x) {
  ints <- as.numeric(strsplit(x, ";")[[1]])
  # remove scale from intensities
  ints <- inverse_scale(ints, scale_factor)
  round(sum(ints)) == 1000
})
if (!all(check_ints_sum)) {
  stop("Bad scaling of the peak lists' intensities. Some spectra do not have the inverse scales intensities summing 1000 (the normalized value).")
}

write.csv(ms_count, file = counts_path, row.names = FALSE)

# concate result to the other count table
if (length(args) == 6)
{
  counts_path2 <- file.path(args[[6]]) 
  if (!file.exists(counts_path2))
  {
    stop("The count file '", counts_path2,
         "' do not exists. Provide a valid path to the counts file.")
  }
  
  ms_count2 <- read.csv(counts_path2, stringsAsFactors = FALSE, comment.char = "", strip.white = TRUE)
  
  ms_count2$peaksList <-  ms_count$peaksList
  ms_count2$peaksInt <-  ms_count$peaksInt
  write.csv(ms_count2, file = counts_path2, row.names = FALSE)
}