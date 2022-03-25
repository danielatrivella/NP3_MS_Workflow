# default parameters and options
bin_size <- 0.05
scale_factor <- 0.5
trim_mz <- TRUE
parallel_cores <- 1
join_isotopic_peaks <- 1 # always on
max_shift <- 200 # maximum mass difference allowed to search for shifted m/z fragmented ion in the cosine computation
# ion_mode <- "+"
# e_mass <- 1.00783

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  stop("Seven arguments must be supplied to create the pairwise similarity table:\n",
       " 1 - Job name;\n",
       " 2 - Path to the MGF file with the spectra to be compared pairwise or directory where the mgfs from the clustering job exists;\n",
       " 3 - Path to the output folder to save the resulting similarity table;\n", 
       " 4 - The bin size to consider two fragmented peaks m/z's the same;\n",
       " 5 - The spectra fragmented peaks scaling method: 0 - ln, 1 - no scale and x > 0 power scale;\n",
       " 6 - A logical indicating if the spectra should be trimmed by the precursor mass;\n",
       " 7 - The maximum difference allowed between precursor m/zs to search for shifted m/z fragments in the cosine computation (max_shift);\n",
       " 8 - (optional) The number of cores to use for parallel processing. At least 2 are needed for parallellization. If 1 disable parallel processing (default).\n",
       call.=FALSE)
} else {
  data_name <- args[[1]]
  path_mgf_dir <- file.path(args[[2]])
  output_path <- file.path(args[[3]]) 
  bin_size <- as.numeric(args[[4]])
  scale_factor <- as.numeric(args[[5]])
  trim_mz <- as.logical(args[[6]])
  max_shift <- as.numeric(args[[7]])
  # ion_mode <- as.numeric(args[[7]])
  
  if (length(args) == 8)
  {
    parallel_cores <- as.numeric(args[[8]])
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
  # retrieve a list with cosine and the matched peaks separated
  # return list with cosine row and the matched peaks row
  np3_cos_matches_list <- normDotProductShiftList(peaks_A = ms2_sample$MZS[[i]], 
                                                  ints_A = ms2_sample$INTS[[i]], 
                                                  mz_A = ms2_sample$PREC_MZ[[i]],
                                                  peaks_B = ms2_sample$MZS[(i+1):n_scans],
                                                  ints_B = ms2_sample$INTS[(i+1):n_scans], 
                                                  mzs_B = ms2_sample$PREC_MZ[(i+1):n_scans],
                                                  bin_size = bin_size,
                                                  max_shift = max_shift)
  list(c(rep(NA, i-1), 1, round(np3_cos_matches_list[[1]], 3)),
       c(rep(NA, i-1), length(ms2_sample$MZS[[i]]), round(np3_cos_matches_list[[2]])))
}

compareSpectraNormDotProductSample <- function(i, peaks_sample_B, ints_sample_B, mzs_sample_B)
{
  np3_cos_matches_list <- normDotProductShiftList(peaks_A = ms2_sample$MZS[[i]], 
                                                  ints_A = ms2_sample$INTS[[i]], 
                                                  mz_A = ms2_sample$PREC_MZ[[i]],
                                                  peaks_B = peaks_sample_B,
                                                  ints_B = ints_sample_B, 
                                                  mzs_B = mzs_sample_B,
                                                  bin_size = bin_size,
                                                  max_shift = max_shift)
  list(round(np3_cos_matches_list[[1]], 3),
       round(np3_cos_matches_list[[2]]))
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
  # order paths by mz and mz count (number of files with that mz range)
  path_mgf <- path_mgf[order(sapply(strsplit(gsub(paste0(data_name, "_[0-9]+_[0-9]+_|.mgf"), 
                                                  "", basename(path_mgf)), "_"),
                                    function(x) { as.numeric(x[1])+as.numeric(x[2])*0.1 }))]
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
if (!is.numeric(max_shift) || max_shift < 0)
  stop("Invalid max_shift parameter value (", max_shift ,
       "). The maximum allowed shift between precursor m/z parameter must be a non negative numeric value. ",
       "Execution aborted.")

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


# start pairwise comparison between the mgf samples
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
                        "bin_size", "max_shift"), envir=environment())
    
    clusterEvalQ(cl, {
      invisible(Rcpp::sourceCpp(file.path(script_path_, 
                                          "dot_product_list.cpp")))})
    # parallel pairwise comparisions
    # separate matched peaks from cosine and save in different tables
    comp_row_sim_matches <- parSapply(cl, 1:(n_scans-1), 
                                      compareSpectraNormDotProductRow)
    # save cosine values
    write.table(bind_cols(comp_row_sim_matches[1,], list(c(rep(NA,n_scans-1), 1.000))), 
                file = file.path(output_path, 
                                 paste0("similarity_table_", data_name, "_tmp.csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE, na = "")
    # save matches, with the num of spectra in the diagonal
    write.table(bind_cols(comp_row_sim_matches[2,], 
                          list(c(rep(NA,n_scans-1), length(ms2_sample$MZS[[n_scans]])))), 
                file = file.path(output_path, 
                                 paste0("similarity_table_matches_", data_name, "_tmp.csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE, na = "")
    rm(comp_row_sim_matches)
    #stopCluster(cl)
  } else if (n_scans > 1) {
    # sequential pairwise comparisions
    comp_row_sim_matches <- sapply(1:(n_scans-1), compareSpectraNormDotProductRow)
    # save cosine values
    write.table(bind_cols(comp_row_sim_matches[1,], list(c(rep(NA,n_scans-1), 1.000))), 
                file = file.path(output_path, 
                                 paste0("similarity_table_", data_name, "_tmp.csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE, na = "")
    # save matches, with the num of spectra in the diagonal
    write.table(bind_cols(comp_row_sim_matches[2,], 
                          list(c(rep(NA,n_scans-1), length(ms2_sample$MZS[[n_scans]])))), 
                file = file.path(output_path, 
                                 paste0("similarity_table_matches_", data_name, "_tmp.csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE, na = "")
    rm(comp_row_sim_matches)
  } else {
    # only one scan, save identical cosine and matches equals the number of peaks
    write.table(c(1.00), 
                file = file.path(output_path, 
                                 paste0("similarity_table_", data_name, "_tmp.csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE, na = "")
    write.table(c(length(ms2_sample$MZS[[n_scans]])), 
                file = file.path(output_path, 
                                 paste0("similarity_table_matches_", data_name, "_tmp.csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE, na = "")
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
      
      # make the pairwise comparisions
      if (parallel_cores > 1 && require(parallel) && n_scans > 1)
      {
        # use existing cluster
        comp_sample_sim_matches <- parSapply(cl, scan_index, 
                                             compareSpectraNormDotProductSample, 
                                             ms2_sample_j$MZS, ms2_sample_j$INTS,
                                             ms2_sample_j$PREC_MZ)
      } else {
        # sequential pairwise comparisions
        comp_sample_sim_matches <- sapply(scan_index, compareSpectraNormDotProductSample, 
                                          ms2_sample_j$MZS, ms2_sample_j$INTS,
                                          ms2_sample_j$PREC_MZ)
      }
      # separate matched peaks from cosine and save in different tables
      if (n_scans_j == 1) #transpose result if only one member
      {
        write.table(t(comp_sample_sim_matches[1,]), 
                    file = file.path(output_path, 
                                     paste0("similarity_table_", data_name, "_tmp.csv")), 
                    sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE, na = "")
        write.table(t(comp_sample_sim_matches[2,]), 
                    file = file.path(output_path, 
                                     paste0("similarity_table_matches_", data_name, "_tmp.csv")), 
                    sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE, na = "")
     } else {
        write.table(comp_sample_sim_matches[1,], 
                    file = file.path(output_path, 
                                     paste0("similarity_table_", data_name, "_tmp.csv")), 
                    sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE, na = "")
       write.table(comp_sample_sim_matches[2,], 
                   file = file.path(output_path, 
                                    paste0("similarity_table_matches_", data_name, "_tmp.csv")), 
                   sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE, na = "")
     }
      rm(comp_sample_sim_matches)
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
    sim_table_tmp <- rbind(matrix(NA, 
                                  nrow = total_spectra-nrow(sim_table_tmp), 
                                  ncol = ncol(sim_table_tmp)), 
                           sim_table_tmp)
    
    write.table(t(sim_table_tmp), file.path(output_path, 
                                            paste0("similarity_table_", data_name, ".csv")),
                row.names = ms2_sample$SCANS, col.names = FALSE, append = TRUE, sep = ",", na = "")
  } else {
    # add heading i = 1 and row.names
    write.table(t(c(paste0("parameters sim pairwise - scale_factor:", scale_factor, 
                           ";bin_size:", bin_size, ";trim_mz:", trim_mz,";max_shift:",max_shift), scans_num)),
                  file = file.path(output_path, paste0("similarity_table_", data_name, ".csv")), 
                  sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(t(sim_table_tmp), file.path(output_path, 
                                            paste0("similarity_table_", data_name, ".csv")),
                row.names = ms2_sample$SCANS, col.names = FALSE, append = TRUE, sep = ",", na = "")
  }
  rm(sim_table_tmp)
  
  # write final number of matches
  sim_table_matches_tmp <- read.csv(file.path(output_path, 
                                      paste0("similarity_table_matches_", data_name, "_tmp.csv")), 
                            header = FALSE, stringsAsFactors = FALSE)
  
  if (i > 1)
  {
    sim_table_matches_tmp <- rbind(matrix(NA, 
                                  nrow = total_spectra-nrow(sim_table_matches_tmp), 
                                  ncol = ncol(sim_table_matches_tmp)), 
                                  sim_table_matches_tmp)
    
    write.table(t(sim_table_matches_tmp), file.path(output_path, 
                                            paste0("similarity_table_matches_", data_name, ".csv")),
                row.names = ms2_sample$SCANS, col.names = FALSE, append = TRUE, sep = ",", na = "")
  } else {
    # add heading i = 1 and row.names
    write.table(t(c(paste0("parameters sim pairwise - scale_factor:", scale_factor, 
                           ";bin_size:", bin_size, ";trim_mz:", trim_mz,";max_shift:",max_shift), scans_num)),
                file = file.path(output_path, paste0("similarity_table_matches_", data_name, ".csv")), 
                sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(t(sim_table_matches_tmp), file.path(output_path, 
                                            paste0("similarity_table_matches_", data_name, ".csv")),
                row.names = ms2_sample$SCANS, col.names = FALSE, append = TRUE, sep = ",", na = "")
  }
  rm(ms2_sample, sim_table_matches_tmp)
}

cat("|\n")
tf <- Sys.time()
cat("    * Done making", ((total_spectra*total_spectra - total_spectra)/2),
        "pairwise comparisions in", round(tf-ti, 2), units(tf-ti), "*\n")
# remove tmp sim table and matches table
unlink(file.path(output_path,
                 paste0("similarity_table_", data_name, "_tmp.csv")))
unlink(file.path(output_path,
                 paste0("similarity_table_matches_", data_name, "_tmp.csv")))

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