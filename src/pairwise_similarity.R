# default parameters and options
bin_size <- 0.05
scale_factor <- 0.5
trim_mz <- TRUE
parallel_cores <- 1
join_isotopic_peaks <- 1 # always on
max_shift <- 200 # maximum mass difference allowed to search for shifted m/z fragmented ion in the cosine computation

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  stop("Seven arguments must be supplied to create the pairwise similarity table:\n",
       " 1 - output_name: Job name;\n",
       " 2 - path_mgf_dir: Path to the MGF file with the spectra to be compared pairwise or directory where the mgfs from the clustering job exists;\n",
       " 3 - output_path: Path to the output folder to save the resulting similarity table;\n", 
       " 4 - bin_size: The bin size to consider two fragmented peaks m/z's the same;\n",
       " 5 - scale_factor: The spectra fragmented peaks scaling method: 0 - ln, 1 - no scale and x > 0 power scale;\n",
       " 6 - trim_mz: A logical indicating if the spectra should be trimmed by the precursor mass;\n",
       " 7 - max_shift: The maximum difference allowed between precursor m/zs to search for shifted m/z fragments in the cosine computation (max_shift);\n",
       " 8 - parallel_cores: (optional) The number of cores to use for parallel processing. At least 2 are needed for parallellization. If 1 disable parallel processing (default).\n",
       call.=FALSE)
} else {
  data_name <- args[[1]]
  path_mgf_dir <- file.path(args[[2]])
  output_path <- file.path(args[[3]]) 
  bin_size <- as.numeric(args[[4]])
  scale_factor <- as.numeric(args[[5]])
  trim_mz <- as.logical(args[[6]])
  max_shift <- as.numeric(args[[7]])
  
  if (length(args) == 8)
  {
    parallel_cores <- as.numeric(args[[8]])
    if (parallel_cores < 1) {
      cat("Invalid parallel_cores value (", parallel_cores, 
              "). The number of cores to be used for parallel processing must be at least 2 ",
              "or 1 for disabling parallelization. The number of cores was set to 1 by default.\n")
      warning("Invalid parallel_cores parameter value.")
      parallel_cores <- 1    
    } else if (parallel_cores > 1 && !require(parallel)) {
      cat("Parallel library is not available, disabling parallelization. ",
              "The number of cores was set to 1 by default.\n")
      warning("Invalid parallel_cores parameter value.")
      parallel_cores <- 1 
    } else if (parallel_cores > 1 && !is.na(detectCores(logical = FALSE)) && parallel_cores >= detectCores(logical = FALSE)) {
      cat("The number of paralell_cores =",parallel_cores,
              " is greater or equal than the number of physical cores =",detectCores(logical = FALSE),
              " in your machine. This is not recommended,",
              "the number of parallel clusters will be limited by the number of physical cores minus two =",
              detectCores(logical = FALSE)-2, ".\n")
      warning("Invalid parallel_cores parameter value.")
      parallel_cores <- detectCores(logical = FALSE)-2
    }
  }
}


## ----load-libs, message = FALSE--------------------------------------------
cat("Loading package data.table and Rcpp functions...\n")
suppressPackageStartupMessages(library(data.table))
# data.table multithreading (fwrite) follows the requested number of cores
setDTthreads(as.integer(max(1, parallel_cores)))
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

# column names of the similarity table, in the order they are written
sim_columns <- c("msclusterID_source", "msclusterID_target", 
                 "similarity_value", "num_matched_peaks")

# read a mgf file and return its ms2 data sorted by scan
# (the data is only reordered when it is not already sorted)
read_mgf_sorted <- function(path)
{
  ms2 <- readMgfPeaksList(path, bin_size = bin_size, 
                          trim_mz = trim_mz, scale_factor = scale_factor,
                          join_isotopic_peaks = join_isotopic_peaks)
  # get scan order
  if (!identical(is.unsorted(ms2$SCANS), FALSE)) {
    scan_index <- order(ms2$SCANS) # get the indexes sorted by scan
    ms2 <- lapply(ms2, function(x) x[scan_index])
  }
  return(ms2)
}

# bind the results of the comparisons into a single data.table with rbindlist
# the results are lists with the columns (or NULL if there is no valid similarity), 
# the NULL results are ignored. Returns a data.table with zero rows if there is no similarity.
bind_sim_results <- function(results, stage)
{
  if (any(vapply(results, inherits, logical(1), what = "try-error"))) {
    stop("ERROR: a failure happened in the pairwise comparison of ", stage, ".")
  }
  return(rbindlist(results))
}

# returns a list with the columns of the sim values >= 0.1 from source msclusterID in
# ms2_sample[[i]] and the targets from i+1 to n_scans (or NULL if there is no valid sim)
# the list have 4 columns: "msclusterID_source","msclusterID_target","similarity_value","num_matched_peaks"
compareSpectraNormDotProductRow <- function(i)
{
  # gets the list of valid cosine >= 0.1, their number of matched peaks, and their idx in [(i+1):n_scans]
  # the cosine is already rounded in 3 decimal places
  np3_cos_matches_list <- normDotProductShiftList(peaks_A = ms2_sample$MZS[[i]], 
                                                  ints_A = ms2_sample$INTS[[i]], 
                                                  mz_A = ms2_sample$PREC_MZ[[i]],
                                                  peaks_B = ms2_sample$MZS[(i+1):n_scans],
                                                  ints_B = ms2_sample$INTS[(i+1):n_scans], 
                                                  mzs_B = ms2_sample$PREC_MZ[(i+1):n_scans],
                                                  bin_size = bin_size,
                                                  max_shift = max_shift)
  # if at least one valid sim, also set source msclusterID - equals SCANS[[i]] -, 
  # convert the target idx (+ i) to their msclusterID and return the columns as a list
  # (no matrix is created, the lists are bound with rbindlist)
  # otherwise return NULL
  number_valid_sim <- length(np3_cos_matches_list[[1]])
  if (number_valid_sim > 0) {
    return(list(msclusterID_source = rep.int(ms2_sample$SCANS[[i]], number_valid_sim),
                msclusterID_target = ms2_sample$SCANS[np3_cos_matches_list[[3]] + i],
                similarity_value = np3_cos_matches_list[[1]],
                num_matched_peaks = np3_cos_matches_list[[2]]))
  } else {
    return(NULL)
  }
}

# call the normDotProductShifList for spec[[i]] of ms2_sample and all spec in sample_B list
# returns a list with the columns of the valid sim (or NULL if there is no valid sim)
compareSpectraNormDotProductSample <- function(i, peaks_sample_B, ints_sample_B, mzs_sample_B, scans_sample_B)
{
  # gets the list of valid cosine >= 0.1, their number of matched peaks, and their idx in sample_B
  # the cosine is already rounded in 3 deciaml places
  np3_cos_matches_list <- normDotProductShiftList(peaks_A = ms2_sample$MZS[[i]], 
                                                  ints_A = ms2_sample$INTS[[i]], 
                                                  mz_A = ms2_sample$PREC_MZ[[i]],
                                                  peaks_B = peaks_sample_B,
                                                  ints_B = ints_sample_B, 
                                                  mzs_B = mzs_sample_B,
                                                  bin_size = bin_size,
                                                  max_shift = max_shift)
  # if at least one valid sim, also set source msclusterID - equals SCANS -, 
  # convert the target idx to their msclusterID and return the columns as a list
  # otherwise return NULL
  number_valid_sim <- length(np3_cos_matches_list[[1]])
  if (number_valid_sim > 0) {
    return(list(msclusterID_source = rep.int(ms2_sample$SCANS[[i]], number_valid_sim),
                msclusterID_target = scans_sample_B[np3_cos_matches_list[[3]]],
                similarity_value = np3_cos_matches_list[[1]],
                num_matched_peaks = np3_cos_matches_list[[2]]))
  } else {
    return(NULL)
  }
}

cat("Number of parallel cores to use in the pairwise comparison:",parallel_cores,"\n")


# if a file was passed, read the mgf from the clean step
# else if a directory was passed, parse the mgfs from the clustering step
if (file.exists(path_mgf_dir) && !dir.exists(path_mgf_dir) && grepl(".mgf$",path_mgf_dir)) {
  path_mgf <- normalizePath(path_mgf_dir)
} else {
  if (!dir.exists(path_mgf_dir) || !any(grepl(".mgf$", list.files(path_mgf_dir))))
  {
    stop("The MGF dir '", path_mgf_dir,
         "' do not exists or do not contains MGF files. ", 
         "Provide a valid path to where the MGF files are located.")
  }
  # get the path to each mgf file inside the mgf dir
  path_mgf <- normalizePath(file.path(path_mgf_dir, list.files(path_mgf_dir)[grepl("_[0-9]+.mgf$", list.files(path_mgf_dir))]))
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
  cat("Invalid trim_mz parameter, it must be a logical indicating if ",
          "the spectra fragmentation should be trimmed by the precursor mass. ",
          "Trim mz set to TRUE by default.\n")
  warning("Invalid trim_mz parameter", call. = FALSE)
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
    cat("Invalid scale factor provided, it must be a numeric value greater or equal to 0: \n", 
            "  - 0 : ln natural logarithm of the intensities\n",
            "  - 1 : no scale\n",
            "  - x : x > 0 pow of the intensities to x. (e.g. x = 0.5 square root)\n",
            "Factor 0.5 (sqrt) will be selected by default.\n")
    warning("Invalid scale factor parameter.", call. = FALSE)
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

# the similarity table is created here with the header, all the similarities are appended to it
sim_file <- file.path(output_path, paste0("similarity_table_", data_name, ".csv"))
writeLines(paste(sim_columns, collapse = ","), con = sim_file)

for (i in seq_along(path_mgf)) {
  cat("=")
  # read the ms2 data (sorted by scan)
  ms2_sample <- read_mgf_sorted(path_mgf[[i]])
  n_scans <- length(ms2_sample$SCANS)
  scan_index <- seq_len(n_scans)
  
  if (i == 1) {
    scans_num <- c(scans_num, ms2_sample$SCANS)
    total_spectra <- total_spectra + n_scans
  }
  
  # compare all spectra of mgf[[i]] pairwise
  # (if there is only one scan, there is nothing to compare)
  if (n_scans > 1)
  {
    if (parallel_cores > 1)
    {
      # parallelized pairwise comparisions
      # run paralellized verison o lapply, use forking
      comp_row_sim_matches <- bind_sim_results(mclapply(seq_len(n_scans-1), 
                                                        compareSpectraNormDotProductRow, 
                                                        mc.cores = parallel_cores),
                                               paste0("the spectra of '", path_mgf[[i]], "'"))
    } else {
      # sequential pairwise comparisions
      comp_row_sim_matches <- bind_sim_results(lapply(seq_len(n_scans-1), 
                                                      compareSpectraNormDotProductRow),
                                               paste0("the spectra of '", path_mgf[[i]], "'"))
    }
    
    if (nrow(comp_row_sim_matches) > 0) {
      # check if all target ID is valid
      if (anyNA(comp_row_sim_matches$msclusterID_target))
      {
        stop("ERROR: a wrong index mapping was detected in the similarity matrix for mgf '",
             path_mgf[[i]],"'. ")
      }
      # save cosine values, appending to the similarity table
      fwrite(comp_row_sim_matches, file = sim_file, 
             append = TRUE, col.names = FALSE, showProgress = FALSE)
    }
    rm(comp_row_sim_matches)
  }
  
  if (i < n_mgf)
  {
    # compare all spectra of mgf[[i]] with all spectra of the others mgf
    for (j in (i+1):n_mgf) {
      # read the ms2 data (sorted by scan)
      ms2_sample_j <- read_mgf_sorted(path_mgf[[j]])
      n_scans_j <- length(ms2_sample_j$SCANS)
      
      if (i == 1) {
        total_spectra <- total_spectra + n_scans_j
        scans_num <- c(scans_num, ms2_sample_j$SCANS)
      }
      
      # make the pairwise comparisons
      # with forking the data of the sample j is shared with the parallel cores (no copy)
      if (parallel_cores > 1 && n_scans > 1)
      {
        # run paralellized verison o lapply, uses forking
        comp_sample_sim_matches <- bind_sim_results(mclapply(scan_index, 
                                                             compareSpectraNormDotProductSample, 
                                                             peaks_sample_B = ms2_sample_j$MZS,
                                                             ints_sample_B = ms2_sample_j$INTS,
                                                             mzs_sample_B = ms2_sample_j$PREC_MZ,
                                                             scans_sample_B = ms2_sample_j$SCANS,
                                                             mc.cores = parallel_cores),
                                                    paste0("the spectra of '", path_mgf[[i]], 
                                                           "' and '", path_mgf[[j]], "'"))
      } else {
        # sequential pairwise comparison
        comp_sample_sim_matches <- bind_sim_results(lapply(scan_index, 
                                                           compareSpectraNormDotProductSample, 
                                                           peaks_sample_B = ms2_sample_j$MZS,
                                                           ints_sample_B = ms2_sample_j$INTS,
                                                           mzs_sample_B = ms2_sample_j$PREC_MZ,
                                                           scans_sample_B = ms2_sample_j$SCANS),
                                                    paste0("the spectra of '", path_mgf[[i]], 
                                                           "' and '", path_mgf[[j]], "'"))
      }
      
      if (nrow(comp_sample_sim_matches) == 0) {
        next
      }
      # save new sim and matches appending to the created sim table
      fwrite(comp_sample_sim_matches, file = sim_file, 
             append = TRUE, col.names = FALSE, showProgress = FALSE)
      rm(comp_sample_sim_matches)
    }
    rm(ms2_sample_j)
  }
  
  rm(ms2_sample)
}

cat("|\n")
tf <- Sys.time()
cat("    * Done making", ((total_spectra*total_spectra - total_spectra)/2),
        "pairwise comparisions in", round(tf-ti, 2), units(tf-ti), "*\n")
