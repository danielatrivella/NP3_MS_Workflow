# default parameters and options
bin_size <- 0.05
scale_factor <- 0.5
trim_mz <- TRUE
parallel_cores <- 1
join_isotopic_peaks <- 1 # always on
# ion_mode <- "+"
# e_mass <- 1.00783

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("Seven arguments must be supplied to create the pairwise similarity table:\n",
       " 1 - Job name;\n",
       " 2 - Path to the MGF file with the spectra to be compared pairwise or directory where the mgfs from the clustering job exists;\n",
       " 3 - Path to the output folder to save the resulting similarity table;\n", 
       " 4 - The bin size to consider two fragmented peaks m/z's the same;\n",
       " 5 - The spectra fragmented peaks scaling method: 0 - ln, 1 - no scale and x > 0 power scale;\n",
       " 6 - A logical indicating if the spectra should be trimmed by the precursor mass;\n",
       " 7 - The number of cores to use for parallel processing. At least 2 are needed for parallellization. If 1 disable parallel processing.\n",
       call.=FALSE)
} else {
  data_name <- args[[1]]
  path_mgf_dir <- file.path(args[[2]])
  output_path <- file.path(args[[3]]) 
  bin_size <- as.numeric(args[[4]])
  scale_factor <- as.numeric(args[[5]])
  trim_mz <- as.logical(args[[6]])
  # ion_mode <- as.numeric(args[[7]])
  
  if (length(args) == 7)
  {
    parallel_cores <- as.numeric(args[[7]])
    if (parallel_cores < 1) {
      warning("Invalid parallel_cores value (", parallel_cores, 
              "). The number of cores to be used for parallel processing must be at least 2 ",
              "or 1 for disabling parallelization. The number of cores was set to 1 by default.")
      parallel_cores <- 1    
    } else if (parallel_cores > 1 && !require(parallel)) {
      warning("Parallel library is not available, disabling parallelization. ",
              "The number of cores was set to 1 by default.")
      parallel_cores <- 1 
    }
  }
}

## ----load-libs, message = FALSE--------------------------------------------
cat("Loading package dplyr and Rcpp functions...\n")
suppressPackageStartupMessages(library(dplyr))
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
Rcpp::sourceCpp(file.path(script_path(), 'read_mgf_peak_list_R.cpp'))
Rcpp::sourceCpp(file.path(script_path(), 'dot_product_list.cpp'))

compareSpectraNormDotProductRow <- function(i)
{
  c(rep(0, i-1), 1, round(normDotProductShiftList(peaks_A = ms2_sample$MZS[[i]], 
                                             ints_A = ms2_sample$INTS[[i]], 
                                             mz_A = ms2_sample$PREC_MZ[[i]],
                                             peaks_B = ms2_sample$MZS[(i+1):n_scans],
                                             ints_B = ms2_sample$INTS[(i+1):n_scans], 
                                             mzs_B = ms2_sample$PREC_MZ[(i+1):n_scans],
                                             bin_size = bin_size), 3))
}

compareSpectraNormDotProductSample <- function(i, peaks_sample_B, ints_sample_B, mzs_sample_B)
{
  round(normDotProductShiftList(peaks_A = ms2_sample$MZS[[i]], 
                           ints_A = ms2_sample$INTS[[i]], 
                           mz_A = ms2_sample$PREC_MZ[[i]],
                           peaks_B = peaks_sample_B,
                           ints_B = ints_sample_B, 
                           mzs_B = mzs_sample_B,
                           bin_size = bin_size), 3)
}

# if a file was passed, read the mgf from the clean step
# else if a directory was passed, parse the mgfs from the clustering step
if (file.exists(path_mgf_dir) && !dir.exists(path_mgf_dir) && grepl(".mgf$",path_mgf_dir)) {
  path_mgf <- path_mgf_dir
} else {
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
}
n_mgf <- length(path_mgf)
if (n_mgf == 0) {
  stop("The MGF dir '", path_mgf_dir,
       "' do not contains MGF files from ", data_name,
       ". Provide a valid path to where the MGF files are located.")
}

if (!dir.exists(output_path))
{
  stop("The output folder '", output_path, 
       "' do not exists. Provide a valid path to where the pairwise similarity table should be saved.")
}
if (!is.logical(trim_mz) || is.na(trim_mz))
{
  warning("Invalid trim_mz parameter, it must be a logical indicating if ",
          "the spectra fragmentation should be trimmed by the precursor mass. ",
          "Trim mz set to TRUE by default.", call. = FALSE)
  trim_mz <- TRUE
}
if (trim_mz) {
  trim_mz <- 1
} else {
  trim_mz <- -1
}
if (!is.numeric(bin_size) || bin_size < 0)
  stop("Invalid bin_size parameter value (", bin_size ,
       "). The bin size parameter must be a non negative numeric value. ",
       "Execution aborted.")

if (!is.na(scale_factor) || !is.null(scale_factor))
{
  scale_factor <- as.numeric(scale_factor[[1]])
  if (scale_factor < 0)
  {
    warning("Invalid scale factor provided, it must be a numeric value greater or equal to 0: \n", 
            "  - 0 : ln natural logarithm of the intensities\n",
            "  - 1 : no scale\n",
            "  - x : x > 0 pow of the intensities to x. (e.g. x = 0.5 square root)\n",
            "Factor 0.5 (sqrt) will be selected by default.", call. = FALSE)
    scale_factor <- 0.5
  }
} else {
  scale_factor <- 1 # no scale
}

# if (!(ion_mode %in% c(-1,1)))
#   stop("The ion_mode arg must be a numeric value indicating the precursor ", 
#        "ion mode. One of the following valid numeric values corresponding ", 
#        "to a ion adduct type: '1' = positive, or '-1' = negative",  call. = FALSE)
# # Set ion mode appropriately
# ion_mode <- ifelse(ion_mode > 0, e_mass, -e_mass)

# create pairwise similarity file, first cell with the parameters
# write.table(t(c(paste0("parameters sim pairwise - scale_factor:", scale_factor, 
#                        ";bin_size:", bin_size, ";trim_mz:", trim_mz), 
#                 as.character(sort(scanIndex(ms2_sample))))), 
#             file = file.path(output_path, paste0("similarity_table_", data_name, ".csv")), sep = ",", 
#             row.names = FALSE, col.names = FALSE)
total_spectra <- 0
cat("\n  * Comparing", data_name, "result spectra pairwise *\n")
# add progress
progress_comp <- seq(from = 1, to = n_mgf, by = 1)
cat(paste0("        |", paste0(rep(" ", length(progress_comp)), collapse = ""), "|\n        |", collapse = ""))
ti <- Sys.time()
scans_num <- c() # will store the name of the scans
for (i in seq_along(path_mgf)) {
  cat("=")
  # read the ms2 data and get the scan indexes
  ms2_sample <- readMgfPeaksList(path_mgf[[i]], bin_size = bin_size, 
                                 trim_mz = trim_mz, scale_factor = scale_factor,
                                 join_isotopic_peaks = join_isotopic_peaks)
  
  # get scan order
  scan_index <- order(ms2_sample$SCANS) # get the indexes sorted by scan
  n_scans <- length(ms2_sample$SCANS)
  
  # order mgf data
  ms2_sample <- lapply(ms2_sample, function(x) x[scan_index])
  scan_index <- seq_along(ms2_sample$SCANS)
  
  if (i == 1) {
    scans_num <- c(scans_num, ms2_sample$SCANS)
    total_spectra <- total_spectra + n_scans
  }
  
  # compare all spectra of mgf[[i]] pairwise
  if (parallel_cores > 1 && require(parallel) && n_scans > 1)
  {
    # parallelized pairwise comparisions
    cl <- makeCluster(parallel_cores)
    script_path_ <- script_path()
    clusterExport(cl, c("ms2_sample", "n_scans","script_path_",
                        "bin_size"), envir=environment())
    
    clusterEvalQ(cl, {
      invisible(Rcpp::sourceCpp(file.path(script_path_, 
                                          "dot_product_list.cpp")))})
    
    write.table(cbind(parSapply(cl, 1:(n_scans-1), 
                          compareSpectraNormDotProductRow), 
                      c(rep(0.00,n_scans-1), 1.00)), 
                file = file.path(output_path, 
                                 paste0("similarity_table_", data_name, "_tmp.csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE)
    #stopCluster(cl)
  } else if (n_scans > 1) {
    # sequential pairwise comparisions
    write.table(cbind(sapply(1:(n_scans-1), compareSpectraNormDotProductRow), 
                      c(rep(0.00,n_scans-1), 1.00)), 
                file = file.path(output_path, 
                                 paste0("similarity_table_", data_name, "_tmp.csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE)
  } else {
    # only one scan
    write.table(c(1.00), 
                file = file.path(output_path, 
                                 paste0("similarity_table_", data_name, "_tmp.csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE)
  }
  
  if (i < n_mgf)
  {
    # compare all spectra of mgf[[i]] with all spectra of the others mgf
    for (j in (i+1):n_mgf) {
      # read the ms2 data and get the scan indexes
      ms2_sample_j <- readMgfPeaksList(path_mgf[[j]], bin_size = bin_size, 
                                       trim_mz = trim_mz, scale_factor = scale_factor,
                                       join_isotopic_peaks = join_isotopic_peaks)
      
      # get scan order
      scan_index_j <- order(ms2_sample_j$SCANS) # get the indexes sorted by scan
      n_scans_j <- length(ms2_sample_j$SCANS)
      
      # order mgf data
      ms2_sample_j <- lapply(ms2_sample_j, function(x) x[scan_index_j])
      
      if (i == 1) {
        total_spectra <- total_spectra + n_scans_j
        scans_num <- c(scans_num, ms2_sample_j$SCANS)
      }
      
      if (parallel_cores > 1 && require(parallel) && n_scans > 1)
      {
        # use existing cluster
        if (n_scans_j == 1) #transpose result if only one member
          write.table(t(parSapply(cl, scan_index, 
                                      compareSpectraNormDotProductSample, 
                                      ms2_sample_j$MZS, ms2_sample_j$INTS, ms2_sample_j$PREC_MZ)), 
                      file = file.path(output_path, 
                                       paste0("similarity_table_", data_name, "_tmp.csv")), 
                      sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
        else
          write.table(parSapply(cl, scan_index, 
                                compareSpectraNormDotProductSample, 
                                ms2_sample_j$MZS, ms2_sample_j$INTS, ms2_sample_j$PREC_MZ), 
                      file = file.path(output_path, 
                                       paste0("similarity_table_", data_name, "_tmp.csv")), 
                      sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      } else {
        # sequential pairwise comparisions
        if (n_scans_j == 1) 
          write.table(t(sapply(scan_index, compareSpectraNormDotProductSample, 
                               ms2_sample_j$MZS, ms2_sample_j$INTS, ms2_sample_j$PREC_MZ)),
                      file = file.path(output_path, 
                                       paste0("similarity_table_", data_name, "_tmp.csv")), 
                      sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
        else
          write.table(sapply(scan_index, compareSpectraNormDotProductSample, 
                             ms2_sample_j$MZS, ms2_sample_j$INTS, ms2_sample_j$PREC_MZ),
                      file = file.path(output_path, 
                                       paste0("similarity_table_", data_name, "_tmp.csv")), 
                      sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
    rm(ms2_sample_j)
  }
  
  if (parallel_cores > 1 && require(parallel) && n_scans > 1)
    stopCluster(cl)
  
  sim_table_tmp <- read.csv(file.path(output_path, 
                                      paste0("similarity_table_", data_name, "_tmp.csv")), 
                            header = FALSE, stringsAsFactors = FALSE)
  
  if (i > 1)
  {
    sim_table_tmp <- rbind(matrix(0, 
                                  nrow = total_spectra-nrow(sim_table_tmp), 
                                  ncol = ncol(sim_table_tmp)), 
                           sim_table_tmp)
    
    write.table(t(sim_table_tmp), file.path(output_path, 
                                            paste0("similarity_table_", data_name, ".csv")),
                row.names = ms2_sample$SCANS, col.names = FALSE, append = TRUE, sep = ",")
  } else {
    # add heading i = 1 and row.names
    write.table(t(c(paste0("parameters sim pairwise - scale_factor:", scale_factor, 
                           ";bin_size:", bin_size, ";trim_mz:", trim_mz), scans_num)),
                  file = file.path(output_path, paste0("similarity_table_", data_name, ".csv")), 
                  sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(t(sim_table_tmp), file.path(output_path, 
                                            paste0("similarity_table_", data_name, ".csv")),
                row.names = ms2_sample$SCANS, col.names = FALSE, append = TRUE, sep = ",")
  }
  
  rm(sim_table_tmp, ms2_sample)
}

cat("|\n")
tf <- Sys.time()
cat("    * Done making", ((total_spectra*total_spectra - total_spectra)/2),
        "pairwise comparisions in", round(tf-ti, 2), units(tf-ti), "*\n")
# remove tmp sim table
unlink(file.path(output_path,
                 paste0("similarity_table_", data_name, "_tmp.csv")))
#  cat("\n  * Making the pairwise similarity table symmetric *\n")
#  ti <- Sys.time()
#  sim_pairwise <- read.csv(file.path(output_path, paste0("similarity_table_", data_name, ".csv")),
#                           header = FALSE, stringsAsFactors = FALSE)
#  sim_pairwise[lower.tri(sim_pairwise)] <- t(sim_pairwise)[lower.tri(sim_pairwise)]
#  write.table(sim_pairwise, file.path(output_path,
#                                      paste0("similarity_table_", data_name, ".csv")),
#            row.names = FALSE, col.names = FALSE, sep = ",")
#  tf <- Sys.time()
#  cat("    * Done in ", round(tf-ti, 2), " ", units(tf-ti), " *\n")