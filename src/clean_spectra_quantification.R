##
# Step 5 - clean optimal clustering of the fragmented clusters
##
options(digits=10) # increase precision
options(readr.show_progress = FALSE)

# default args values
mz_tol <- 0.025
sim_tol <- 0.55 # sim to join
rt_tol <- 2
bin_size <- 0.05
scale_factor <- 0.5
trim_mz <- TRUE 
max_shift <- 200 
parallel_cores <- 1 # no paralleliztion
mz_rt_digits <- 4  # number of digits to round the mzConsensus and the retention times
ion_mode <- "+" # + or -
# BFLAG_cutoff_factor is A positive numeric value to scale the interquartile range (IQR) of the 
# blank spectra basePeakInt distribution and allow spectra with a basePeakInt
# value below this distribution median plus IQR*BFLAG_cutoff_factor to be joined with 
# a blank spectrum without relying on the similarity value.
# Or FALSE to disable it.
# The IQR is the range between the 1st quartile (25th quantile) and the 3rd 
# quartile (75th quantile). The spectra with a basePeakInt value <= median + 
# IQR*bflag_cutoff of the blank spectra basePeakInt distribution and BFLAG TRUE 
# will be joined to a blank spectrum in the clean step. This cutoff will affect 
# the spectra with BFLAG TRUE that would not get joined to a blank spectra when 
# relying only on the similarity cutoff. This is a turn around to the fact that 
# blank spectra have low quality spectra and thus can not fully rely on the similarity values.
BFLAG_cutoff_factor <- 1.5 
NOISE_cutoff <- 0 # the absolute base peak intesity minimum value that a spectra must have to be kept for cleaning

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  cat("Three arguments are mandatory to execute Step 5 - clean the clustering counts:\n", 
       " 1 - output_path: Path to the output data folder, inside the outs directory of the clustering result folder. ", 
       "It should contain the mgf folder, the counts_table folder with the peak area count CSV and the spectra count CSV. The data name will be extracted from here.;\n", 
       " 2 - path_batch_metadata: Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks; if this is a join_jobs result the metadata must be set as '-1';\n",
       " 3 - processed_data_path: Path to the pre processed data folder were the MGFs were created. Used to compute the peak areas. If this is a joining job this parm will be ignored;\n",
       " 4 - mz_tol: The precursor m/z tolerance in Da;\n",
       " 5 - sim_tol: The similarity tolerance to JOIN mass dissipation spectra;\n",
       " 6 - rt_tol: The retention time tolerance to enlarge the peaks before comparisions;\n",
       " 7 - bin_size: The bin size in Da to consider two fragmented peaks m/z's the same;\n",
       " 8 - scale_factor: The scale factor to be used to aggregated the peaks of joined clusters;\n",
       " 9 - trim_mz: A logical indicating if the spectra should be trimmed by the precursor mass for pairwise comparison;\n",
       " 10 - max_shift: The maximum difference allowed between precursor m/zs to search for shifted m/z fragments in the cosine computation of the pairwise comparion;\n", 
       " 11 - ion_mode: The ionization mode, one of -1 or 1 (default to 1);\n",
       " 12 - BFLAG_cutoff_factor: The bflag_cutoff a positive numeric or FALSE to disable it. If positive numeric, allow to join the spectra with BFLAG=T and a low basePeakInt (<= median + bflag_cutoff*IQR) with a blank spectra without relying on the similarity value;\n",
       " 13 - NOISE_cutoff: The NOISE_cutoff an absolute positive numeric or 0 to disable it. If positive numeric, allow to remove the spectra with a low basePeakInt  (basePeakInt < NOISE_cutoff) before the clean step;\n",
       " 14 - parallel_cores: (optional) The number of cores to use for parallel processing. At least 2 are needed for parallellization. If 1 disable parallel processing (default).\n\n")
  stop("Wrong number of arguments. Minimum of 3 (default values to the rest) or all arguments must be informed.",call.=FALSE)
} else {
  output_path <- file.path(args[[1]])
  if (!dir.exists(output_path))
  {
    stop("The job output folder '", output_path, 
         "' do not exists. Provide a valid path to where the job final result is located.")
  }
  output_name <- basename(output_path)
  
  path_area_count <- file.path(output_path, "count_tables", paste0(output_name,"_peak_area.csv"))
  if (!file.exists(path_area_count))
  {
    stop("The count file '", path_area_count,
         "' do not exists. Provide a valid output path to where csv file with the count ",
         "of peak area is located.")
  }
  
  path_spectra_count <- file.path(output_path, "count_tables", paste0(output_name,"_spectra.csv"))
  if (!file.exists(path_spectra_count))
  {
    stop("The count file '", path_spectra_count,
         "' do not exists. Provide a valid output path to where csv file with the count ",
         "of spectra is located.")
  }
  
  path_batch_metadata <- args[[2]]
  if (path_batch_metadata == "-1") {
    # this is a joinig job process, the metadata should be extracted from the 
    # default path to the original samples metadata - located in the output path
    joining_jobs <- TRUE
    path_batch_metadata <- file.path(output_path, "../..", 
                                     "original_samples_METADATA.csv")
    if (!file.exists(path_batch_metadata))
    {
      stop("The metadata with the original samples from the joined jobs '", path_batch_metadata,
           "' do not exists. Provide a valid output path to where the csv file with  ",
           "the default samples metadata is located.")
    }
  } else {
    joining_jobs <- FALSE
    path_batch_metadata <- file.path(args[[2]])
    if (!file.exists(path_batch_metadata))
    {
      stop("The CSV batch metadata file '", path_batch_metadata, 
           "' do not exists. Provide a valid path to where the metadata is located.")
    }
  }
  
  if (!joining_jobs) {
    # ignore the pre processed data dir when joining jobs, 
    # this info will be extracted from the jobs metadata
    processed_data_path <- file.path(args[[3]])
    if (!dir.exists(processed_data_path))
    {
      stop("The processed data folder '", processed_data_path, 
           "' do not exists. Provide a valid path to where the pre processed MGFs are located.")
    }
  }
  
  if (length(args) > 13) {
    mz_tol <- as.numeric(args[[4]])
    if (is.na(mz_tol))
      stop("The m/z tolerance must be a numeric value. Wrong value informed: ",mz_tol)
    sim_tol <- as.numeric(args[[5]]) # sim to join
    if (is.na(sim_tol))
      stop("The similarity tolerance must be a numeric value. Wrong value informed: ",sim_tol)
    rt_tol <- as.numeric(args[[6]])
    if (is.na(rt_tol))
      stop("The retention time tolerance must be a numeric value. Wrong value informed: ",rt_tol)
    bin_size <- as.numeric(args[[7]])
    if (is.na(bin_size))
      stop("The bin size must be a numeric value. Wrong value informed: ",bin_size)
    
    scale_factor <- as.numeric(args[[8]])
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
    
    trim_mz <- as.logical(args[[9]])
    if (!is.logical(trim_mz) || is.na(trim_mz))
    {
      cat("Invalid trim_mz parameter, it must be a logical indicating if ",
          "the spectra peaks should be trimmed by the precursor mass for comparison. ",
          "Trim mz set to TRUE by default.\n")
      warning("Invalid trim_mz parameter", call. = FALSE)
      trim_mz <- TRUE
    }
    if (trim_mz) {
      trim_mz <- 1
    } else {
      trim_mz <- -1
    }
    max_shift <- as.numeric(args[[10]])
    if (!is.numeric(max_shift) || max_shift < 0) {
      stop("Invalid max_shift parameter value (", max_shift ,
           "). The maximum allowed shift between precursor m/z parameter must be a non negative numeric value. ",
           "Execution aborted.")
    }
    
    ion_mode <- as.numeric(args[[11]])
    if (!(ion_mode %in% c(-1,1)))
      stop("The ion_mode arg must be a numeric value indicating the precursor ", 
           "ion mode. One of the following valid numeric values corresponding ", 
           "to a ion adduct type: '1' = positive, or '-1' = negative",  call. = FALSE)
    ion_mode <- ifelse(ion_mode > 0, "+", "-")
    
    BFLAG_cutoff_factor <- args[[12]]  # check if it is FALSE
    if (BFLAG_cutoff_factor != "FALSE") { 
      # check if it is numeric
      BFLAG_cutoff_factor <- as.numeric(args[[12]]) 
      if (is.na(BFLAG_cutoff_factor) || BFLAG_cutoff_factor < 0)
        stop("The BFLAG_cutoff arg must be a positive numeric value or FALSE to disable it. ",
             "If positive numeric value, indicates that the features with low basePeakInt and BFLAG TRUE should be joined with ",
             "a blank mz without relying on the similarity value.",  call. = FALSE)
    } else {
      BFLAG_cutoff_factor <- as.logical(args[[12]]) 
    }
    
    NOISE_cutoff <- as.numeric(args[[13]])   # convert to absolute value
    # check if it is numeric
    if (is.na(NOISE_cutoff) || NOISE_cutoff < 0)
      stop("The NOISE_cutoff arg must be a positive numeric value or 0 to disable it. ",
           "If positive numeric value, the features with a low basePeakInt < NOISE_cutoff will be removed before the clean step.",  
           call. = FALSE)
    
    
    parallel_cores <- as.numeric(args[[14]])
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
  } else {
    cat("- Using the default value in the rest of the arguments (4 to 14). Inform all the arguments if this is not wanted.")
  }
}
#print(args)

## ----load-libs, message = FALSE--------------------------------------------
cat("Loading packages Matrix, stringi, readr, dplyr, CPP functions...\n")
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
library(Matrix)
library(stringi)

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

#Rcpp::sourceCpp(file.path(script_path(), 'triangular_matrix_R.cpp'))
Rcpp::sourceCpp(file.path(script_path(), 'cluster_disjoint_set_union_find_R.cpp'))
Rcpp::sourceCpp(file.path(script_path(), 'read_mgf_peak_list_R.cpp'))
Rcpp::sourceCpp(file.path(script_path(), 'norm_dot_product.cpp'))
source(file.path(script_path(), "pairwise_similarity_comparison_sets.R"))
source(file.path(script_path(), "count_peak_area.R"))
source(file.path(script_path(), "read_metadata_table.R"))
source(file.path(script_path(), "writeMgfData_NP3.R"))
source(file.path(script_path(), "final_report",  "compute_fragmented_cluster.R"))
# commented for debug
# source(file.path("src/pairwise_similarity_comparison_sets.R"))
# Rcpp::sourceCpp(file.path('src/cluster_disjoint_set_union_find_R.cpp'))
# Rcpp::sourceCpp(file.path('src/read_mgf_peak_list_R.cpp'))
# Rcpp::sourceCpp(file.path('src/norm_dot_product.cpp'))
# source(file.path("src/count_peak_area.R"))
# source(file.path("src/read_metadata_table.R"))
# source(file.path("src/writeMgfData_NP3.R"))
# source(file.path("src/final_report",  "compute_fragmented_cluster.R"))

inverse_scale <- function(ints, scale_factor)
{
  if (scale_factor == 0) # log scale was applied
  {
    # inverse of 1.0 + log(1.0 + multVal * ints[i]))
    ints <- expm1(ints-1)
  } else {
    ints <- ints^(1/scale_factor)
  }
  return(round(ints,5))
}

scaleInts <- function(ints, scale_factor)
{
  if (scale_factor == 0) # log scale was applied
  {
    # following 1.0 + log(1.0 + multVal * ints[i]))
    ints <- 1.0 + log(1.0+ints)
  } else {
    ints <- ints ^ scale_factor
  }
  return(round(ints,5))
}

# function to merge the count table by column
# rules to merge each column
merge_counts <- function(col_name, x)
{
  switch(col_name,
         msclusterID = min(as.numeric(x$msclusterID)), # TODO get which max basePeakInt or keep the smaller id?
         origin_cluster = as.numeric(x$msclusterID[[1]]),
         numJoins =  sum(x$numJoins) + nrow(x) - 1,
         mzConsensus=,rtMean=
           round(weighted.mean(x[[col_name]], x$basePeakInt), mz_rt_digits),
         rtMin=ifelse(any(x$rtMin == 0 & x$rtMax == 1000000), 0, round(weighted.mean(x[[col_name]], x$basePeakInt), mz_rt_digits)),
         rtMax=ifelse(any(x$rtMin == 0 & x$rtMax == 1000000), 1000000, round(weighted.mean(x[[col_name]], x$basePeakInt), mz_rt_digits)),
         numSpectra =,BLANKS_TOTAL =,BEDS_TOTAL=,CONTROLS_TOTAL=,sumInts = sum(x[[col_name]]),
         basePeakInt=max(as.numeric(x[[col_name]])),
         BEDFLAG=,BFLAG =,CFLAG = any(as.logical(x[[col_name]])),
         HFLAG=,DESREPLICATION=,scans=,joinedOriginJobsID=,joinedJobsID=,msclusterID_integrative=
           ifelse(any(!is.na(x[[col_name]])), # if there is a not NA value paste it
                  paste(x[[col_name]][!is.na(x[[col_name]])], collapse = ";"), 
                  NA),
         peaksList=,peaksInt= # only keep the most intense spectra as the representative
           x[[col_name]][which.max(as.numeric(x$basePeakInt))],
         peakIds= # if there is a not NA value paste it, remove duplicates
           ifelse(any(!is.na(x[[col_name]])), 
                  paste(unique(unlist(strsplit(x[[col_name]][!is.na(x[[col_name]])], ";"))), collapse = ";"),
                  NA),
         joinedIDs = paste(c(x$msclusterID[is.na(x[[col_name]])],
                             x[[col_name]][!is.na(x[[col_name]])]), collapse = ";"),
         ifelse(is.numeric(x[[col_name]]) && !startsWith(col_name, 'gnps_'), # do not sum gnps results, concatenate them
                sum(x[[col_name]], na.rm = TRUE), # sum all the counts cols (_area and _spectra) AND the fragmented_clusters
                ifelse(any(!is.na(x[[col_name]])), # if there is a not NA value paste it, delim = ;
                       paste(x[[col_name]][!is.na(x[[col_name]])], collapse = ";"), 
                       NA))) # cocatenate string fields (e.g. from gnps)
}

# TODO not used anymore, remove or use in another place
# function to remove from a similarity table file in input_sim_table_file 
# the rows containing any of the IDs present in the remove_msclsuterIDs list
# and store the remaining table to output_sim_table_file
# uses awk with regex for better performance and almost no RAM overload
sim_table_filtering_out <- function(remove_msclusterIDs = c(1,2,3), 
                                    input_sim_table_file, output_sim_table_file) {
  # concate the IDs to be removed using pipe in a regex format
  regex_pattern_removal <- paste0("'^(", paste(remove_msclusterIDs, 
                                              collapse = "|"), ")$'")
  # Build the awk arguments
  # -F, tells awk it is a CSV file
  # $1 !~ pat ensures column 1 does not match the list
  # $2 !~ pat ensures column 2 does not match the list
  awk_args <- c(
    "-F,",
    "-v", paste0("pat_rm=", regex_pattern_removal),
    "'NR == 1 || (($1 !~ pat_rm) && ($2 !~ pat_rm))'" ,
    input_sim_table_file)
  # Stream directly via system2
  status <- system2(
    command = "awk",
    args = awk_args,
    stdout = output_sim_table_file)
  # Check execution status
  if (status != 0) {
    stop("ERROR: Awk filtering failed: ", paste(status, collapse = "\n"),
         ". Could not remove some msclusterIDs from the similarity table '",
         input_sim_table_file, "'. Check for warnings.")
  }
}

# format a awk command to read the rows of a given similarity table CSV in sim_file
# that contains data between the msclusterID in the list_msclusterIDs_filter
awk_cmd_filter_sim_table_msclusterIDs <- function(list_msclusterIDs_filter,
                                                  sim_file) {
  # concate the IDs to be filtered using pipe in a regex format
  regex_pattern_keep <- paste0("'^(", paste(list_msclusterIDs_filter, 
                                               collapse = "|"), ")$'")
  # Build the awk arguments
  # -F, tells awk it is a CSV file
  # $1 ~ pat ensures column 1 matches the list
  # $2 ~ pat ensures column 2 matches the list
  awk_cmd <- paste("awk -F,","-v", paste0("pat_keep=", regex_pattern_keep),
    "'NR == 1 || (($1 ~ pat_keep) && ($2 ~ pat_keep))'",sim_file)
  return(awk_cmd)
}

# function to aggregate the sim table using awk to replace joined ids with the representative ID
# only return the rows that are not a selfloop and are not duplicated
# format a awk command to read the rows of a given similarity table CSV in sim_file
# that may contain data from one of the msclusterIDs in the list_msclusterIDs_sub
# and replace the matched ids with the given representative ID in the respective col
# return the rows with msclusterID_from != msclusterID_to (no selfloop) and not duplicated
# @ list_msclusterIDs_sub is a string with a list of msclusterIDs separated by pipe "|"
# @ representative_msclusterID is the msclusterID to be used to replace the ones present in list_msclusterIDs_sub
# @ sim_file is the path to the similarity file to make the replacements
sim_table_replace_msclusterIDs <- function(list_msclusterIDs_sub, 
                                           representative_msclusterID,
                                           sim_file) {
  # concate the IDs to be filtered using pipe in a regex format
  regex_pattern_sub <- paste0("'^(", list_msclusterIDs_sub, ")$'")
  # Build the awk arguments
  # -F, tells awk it is a CSV file
  # $1 ~ pat ensures column 1 matches the list
  # $2 ~ pat ensures column 2 matches the list
  # $1 != $2 && (!seen[$1,$2]++) ensures column 1 is different from column 2 and is not duplicated
  awk_cmd <- paste("-F,","-v", paste0("pat_sub=", regex_pattern_sub),
                   paste0("'BEGIN {FS=OFS=\",\"} {if ($1 ~ pat_sub) $1 = ",representative_msclusterID,
                          "; if ($2 ~ pat_sub) $2 = ",representative_msclusterID,
                          "; if ($1 != $2 && (!seen[$1,$2]++)) print}'"),
                   sim_file)
  
  # output to a tmp file
  sim_tmp_file <- sub(pattern=".csv", replacement="_tmp.csv", x=sim_file)
  # Stream directly via system2
  status <- system2(
    command = "awk",
    args = awk_cmd,
    stdout = sim_tmp_file)
  # Check execution status
  if (status != 0) {
    stop("ERROR: Awk replacement failed: ", paste(status, collapse = "\n"),
         ". Could not replace some msclusterIDs from the similarity table '",
         sim_file, "'. Check for warnings.")
  }
  # if status ok, mv tmp to sim_file
  file.rename(from=sim_tmp_file, sim_file)
  return(status)
}

# function to extract the pairs of msclusterIDs that share a sample occurrence
# from the df with the samples quantification columns of each mz in list_msclusterIDs
# return a df with the mzs ID that share at least one sample with two columns: 
# from and to with the detected pairs; use Matrix library for tcrossprod of a sparseMatrix
get_mz_pairs_share_samples <- function(df_samples_quantification, list_msclusterIDs)
{
  # Convert to a sparse matrix and calculate the row-by-row intersections
  sparse_mat <- as(df_samples_quantification, "sparseMatrix")
  shared_counts <- tcrossprod(sparse_mat)
  # Keep only the lower triangle to get unique pairs (removes self-matches and mirrors)
  shared_counts <- tril(shared_counts, k = -1)
  # Extract the row indices and coordinates where the shared count is >= 1
  summary_data <- summary(shared_counts)
  
  # Build the final clean data frame with your actual row names, 
  # from j to i to keep smaller IDs in the "from" list
  result_pairs <- data.frame(
    mz_ID_from = list_msclusterIDs[summary_data$j],
    mz_ID_to = list_msclusterIDs[summary_data$i]
  )
  return(result_pairs)
}

tstart <- Sys.time()
# read count files
ms_spectra_count <- suppressMessages(read_csv(path_spectra_count, 
                                              guess_max = 5000,
                                              col_types = cols(.default="?", 
                                                               msclusterID="i")))

# disable the BFLAG cutoff if there is no blank sample
if ('BLANKS_TOTAL' %in% names(ms_spectra_count)) {
  blanks_flag <- TRUE
} else {
  blanks_flag <- FALSE
}
# if no blanks, disable BFLAG cutoff
if (!blanks_flag) {
  BFLAG_cutoff_factor <- FALSE
  BFLAG_cutoff <- -1
  cat("\nBFLAG cutoff : disabled\n")
} else if (BFLAG_cutoff_factor != FALSE) {
  # compute the summary of the basePeakInt distribution for blank mzs if any, or using the complete distribution
  summary_basePeakInt <- summary(ms_spectra_count$basePeakInt[ms_spectra_count$BLANKS_TOTAL > 0])
  # BFLAG_cutoff is a numeric value, compute the cutoff as the basePeakInt median + bflag_cutoff_factor*(q75-q25) of the blank mzs
  BFLAG_cutoff <- summary_basePeakInt[['Median']] + BFLAG_cutoff_factor*(summary_basePeakInt[['3rd Qu.']]-summary_basePeakInt[['1st Qu.']]) 
  cat("\nBFLAG cutoff : MS2 basePeakInt median +",BFLAG_cutoff_factor,"* (q75-q25) =", BFLAG_cutoff,"\n")
  remove(summary_basePeakInt)
} else {
  if (!exists('BFLAG_cutoff')) {
    BFLAG_cutoff <- -1
    cat("\nBFLAG cutoff : disabled\n")
  }
}

if (NOISE_cutoff == 0) {
  cat("\nNOISE cutoff : disabled\n")
} else {
  cat("\nNOISE cutoff : ", NOISE_cutoff," absolute MS2 basePeakInt\n")
}

# apply the noise cutoff filter before the clean step or the similarities computation 
# based on the basePeakInt < noise_cutoff
if (NOISE_cutoff > 0 && any(ms_spectra_count$basePeakInt < NOISE_cutoff)) {
  cat("\n** Applying the noise cutoff and removing all spectra with a basePeakInt value <",
      NOISE_cutoff,"**\n")
  t1 <- Sys.time()
  spectra_to_keep <- (ms_spectra_count$basePeakInt >= NOISE_cutoff)
  msclusterID_to_remove <- ms_spectra_count$msclusterID[!spectra_to_keep]
  ms_spectra_count <- ms_spectra_count[spectra_to_keep,]
  # recompute number of scans
  nscans <- nrow(ms_spectra_count)
  tf <- Sys.time()
  cat("\n  * Done removing", sum(!spectra_to_keep), "noise spectra in", 
      round(tf-t1, 2), units(tf-t1), "*\n\n")
  remove(spectra_to_keep, msclusterID_to_remove)
  if (nscans == 0) {
    stop("Removed all spectra with the defined noise cutoff. Reset this parameter value if this is not wanted.")
  }
}

cat("\n** Retrieving the groups of putative isomers m/z with 3*mz_tol and ",
    "max retention time range - all the possible fragmented clusters to be compared pairwise **\n")

# first retrieve the groups of putative isomers with extended tolerances
mz_isomers_groups <- bind_cols(ms_spectra_count[,"msclusterID"], compute_fragmented_clusters(ms_spectra_count, 
                                                 max(ms_spectra_count$rtMax) - min(ms_spectra_count$rtMin), 
                                                 3*mz_tol,
                                                 putative_origin=TRUE, 
                                                 table_type="isomers from Step 3 - Clustering ",
                                                 add_space_verbose="  "))
number_mz_no_isomers <- sum(is.na(mz_isomers_groups$origin_cluster))
cat("\n    - Total number of m/z without isomers and percentage overall m/zs:", 
    number_mz_no_isomers, "(",round(100*number_mz_no_isomers/nrow(ms_spectra_count),2), "% )\n")

# make pairwise comparison of all possible isomers and store only the sim above the cutoff
cat("\n** Creating the similarity table for the cleaning step - filter similarities >=",
    sim_tol,"**\n")
# Define file paths
mgf_clust_file  <- file.path(output_path,"mgf/",paste0(output_name,"_all.mgf"))
if (!file.exists(mgf_clust_file))
{
  stop("The MGF file '", mgf_clust_file,
       "' from the clustering step does not exists. ",
       "Rerun the job or provide a valid path to the output directory.")
} 
# new aggMax sim file of clean
sim_file <- file.path(output_path, 
                      "molecular_networking/similarity_tables", 
                      paste0("similarity_table_", output_name, "_aggMax.csv"))
sim_clust_file <- file.path(output_path, 
                      "molecular_networking/similarity_tables", 
                      paste0("similarity_table_", output_name, ".csv"))
# call pairwise similarity for the clusters_origin
n_sim_computed <- pairwise_similarity_comparison_sets(df_mscluster_sets = mz_isomers_groups, 
                                    path_mgf=mgf_clust_file, 
                                    sim_output_path=sim_clust_file, 
                                    set_name="putative isomers groups from Step 3 Clustering",
                                    sim_cutoff=sim_tol,
                                    bin_size=bin_size, 
                                    scale_factor=scale_factor,
                                    trim_mz=(trim_mz==1),
                                    max_shift=max_shift,
                                    join_isotopic_peaks=1, # always on
                                    parallel_cores=parallel_cores)
# copy sim clust table to the sim aggMax table, which will be aggregated during clean
if (!file.copy(from=sim_clust_file, to=sim_file, overwrite=TRUE))
  stop("Could not copy the similarity file '",sim_clust_file,"' to the aggMax file.")
remove(sim_clust_file)
# read the metadata samples table
batch_metadata <- readMetadataTable(path_batch_metadata)

# remove any column from tremolo or gnps
ms_spectra_count <- ms_spectra_count[, !startsWith(names(ms_spectra_count), "tremolo_")]
ms_spectra_count <- ms_spectra_count[, !startsWith(names(ms_spectra_count), "gnps_")]
# reset columns used in the joins
ms_spectra_count$joinedIDs <- NA
ms_spectra_count$numJoins <- 0
# multiple the peakInts by the sumInts to correctly scale the joined spectra peaks
ms_spectra_count$peaksInt <- sapply(seq_len(nrow(ms_spectra_count)), function(i){
  peaksInt <- as.numeric(strsplit(ms_spectra_count$peaksInt[[i]], ";")[[1]])
  # remove scale from intensities
  peaksInt <- inverse_scale(peaksInt, scale_factor)
  # weight by the sumInts
  peaksInt <- peaksInt*ms_spectra_count$sumInts[[i]]
  paste0(peaksInt, collapse=";")
})

not_count_columns <- which(!(names(ms_spectra_count) %in% paste0(batch_metadata$SAMPLE_CODE, "_spectra"))) 
count_columns <- which(names(ms_spectra_count) %in% paste0(batch_metadata$SAMPLE_CODE, "_spectra"))

# order the spectra table by the msclusterID column
ms_spectra_count <- ms_spectra_count[order(ms_spectra_count$msclusterID),]
nscans <- nrow(ms_spectra_count)

cat("\n** Cleanning counts, removing redudancies and aggregating data **\n")

# compute the number of fragmented clusters in the clustering counts
ms_spectra_count$fragmented_clusters <- NULL
ms_spectra_count$origin_clusters <- NULL
ms_spectra_count <- bind_cols(ms_spectra_count, compute_fragmented_clusters(ms_spectra_count, 
                                                                            rt_tol, mz_tol, 
                                                                            putative_origin=TRUE, 
                                                                            table_type="Step 3 clustering"))

cat("\n  * Starting Clean Optimal Clustering with Greedy Heuristics in 10 Rounds *\n")
if (!dir.exists(file.path(output_path, "count_tables", "clean")))
  dir.create(file.path(output_path, "count_tables", "clean"), showWarnings = FALSE)
t0 <- Sys.time()
# start on step 1, analys blank and not blank together
step_join <- 1
num_joins_total <- 0
# stores the number of joins by step and the total number of joins by cluster
num_joins_last_step <- total_num_join_clusters <- ms_spectra_count$numJoins

# repeat the clean in rounds, after the first round look just for the joined clusters 
repeat
{
  t1 <- Sys.time()
  cat("\n  * Round", step_join, "*\n")
  num_joins <- 0
  
  # start each round retrieving the msclusters with fragmented clusters, blank or not blank
  if (step_join == 1)
  {
    # in the first round, retrieve all mz with fragmented clusters
    msclusterID_to_search <- ms_spectra_count$msclusterID[ms_spectra_count$fragmented_clusters > 0]
  } else {
    # after the first round, retrieve the fragmented clusters or
    # mz with fragmented clusters that have numJoins > 0 (which changed)
    # and use their unique origin cluster
    msclusterID_to_search <- unique(ms_spectra_count$origin_cluster[
      ms_spectra_count$fragmented_clusters != 0 & ms_spectra_count$numJoins > 0])
  }
  
  nsearch <- length(msclusterID_to_search)
  if (nsearch == 0) {
    cat("    * Done. No more fragmented clusters to join. *\n")
    break()
  }
  
  progress_joins <- msclusterID_to_search[unique(trunc(c(seq(from = 1, to = nsearch, by = nsearch/25), nsearch)))]
  n_progress <- length(progress_joins)
  cat(paste0("    |", paste0(rep(" ", n_progress), collapse = ""), "|\n    |", collapse = ""))

  # iterate directly on the m/z with fragmented clusters and 
  # after the first round, only iterate over the ones with numJoins>0
  # msclusterID_to_search is the list of msclusterID from the m/z with fragmented clusters or their origin clusters
  # and not their idx, because the idx changes during cleaning
  for (scan_num in msclusterID_to_search) 
  {
    i <- which(ms_spectra_count$msclusterID == scan_num)
    #cat("i: ", i, "\n")
    if (n_progress > 0 && scan_num == progress_joins[[1]]) {
      progress_joins <- progress_joins[-1]
      n_progress <- n_progress - 1
      cat("=")
    }
    if (is.na(i)) {
      # this msclusterID was already joined, go to the next
      next()
    }
    
    # get next cluster scan number and cluster info
    cluster <- ms_spectra_count[i,]
    origin_cluster <- ms_spectra_count[i,"origin_cluster"][[1]]
    # get the cluster peak from the fragmented clusters origin,
    # this includes the current cluster row
    cluster_peak <- ms_spectra_count[ms_spectra_count$origin_cluster %in% origin_cluster,]
    
    # create lists of m/z pairs to be joined
    join_from_msclusterIDs <- c()
    join_to_msclusterIDs <- c()
    
    # if BFLAG_cutoff is enabled, join blanks first if any
    # for the mz with BFLAG TRUE in the current peak and basePeakInt < BFLAG_cutoff, 
    # join them to the blank m/z (BLANKS_TOTAL > 0) with the greater basePeakInt 
    if (BFLAG_cutoff >= 0 && any(cluster_peak$BLANKS_TOTAL > 0)) {
      bflag_mzs <- (cluster_peak$BFLAG)
      blank_mzs <- (cluster_peak$BLANKS_TOTAL > 0)
      if (any(cluster_peak[bflag_mzs, "basePeakInt"] < BFLAG_cutoff)) {
        bflag_mz_to_join <- cluster_peak$msclusterID[bflag_mzs][cluster_peak[bflag_mzs, "basePeakInt"] < BFLAG_cutoff]
        blank_mz_to_receive <- cluster_peak$msclusterID[blank_mzs][which.max(cluster_peak$basePeakInt[blank_mzs])]
        join_from_msclusterIDs <- bflag_mz_to_join[bflag_mz_to_join != blank_mz_to_receive]
        join_to_msclusterIDs <- rep_len(blank_mz_to_receive, length(join_from_msclusterIDs))
        # check if all ids in the from list are smaller than the ones in the to list,
        # if not, reverse the pairs where this is not true
        if (any(join_from_msclusterIDs > join_to_msclusterIDs)) {
          reverse_ids <- (join_from_msclusterIDs > join_to_msclusterIDs)
          join_from_msclusterIDs_original <- join_from_msclusterIDs
          join_from_msclusterIDs[reverse_ids] <- join_to_msclusterIDs[reverse_ids]
          join_to_msclusterIDs[reverse_ids] <- join_from_msclusterIDs_original[reverse_ids]
          rm(join_from_msclusterIDs_original)
        }
      }
    }
    
    # if all mz from the current cluster peak were already selected for join 
    # based on the bflag cutoff, skip searching for new pairs of mz to join;
    # otherwise, if there is at least one not covered mz try to find a pair for it
    # based on the similarities or the shared peakIds
    if (length(unique(c(join_from_msclusterIDs, join_to_msclusterIDs))) == nrow(cluster_peak)) 
    {
      new_clusters_members <- list(cluster_peak$msclusterID)
      new_clusters_size <- nrow(cluster_peak)
      
    } 
    else
    {
      # get valid sim values between the mz from the the cluster_peak 
      sim_cluster_peak <- as_tibble(data.table::fread(
        cmd=awk_cmd_filter_sim_table_msclusterIDs(list_msclusterIDs_filter = cluster_peak$msclusterID,
                                                  sim_file=sim_file)))
      # sim_cluster_peak contains all mz pairs that should be joined
      if (nrow(sim_cluster_peak) > 0)
      {
        join_from_msclusterIDs <- c(join_from_msclusterIDs,
                                    sim_cluster_peak$msclusterID_source)
        join_to_msclusterIDs <- c(join_to_msclusterIDs,
                                  sim_cluster_peak$msclusterID_target)
      }
      
      # now try to find the pairs that share a MS1 peak
      # first get the pairs of mz that appear in at least one common sample
      # and for those check the peakIDs from the matching samples
      matrix_cluster_peak_samples_count <- as.matrix(cluster_peak[,count_columns])
      pairs_mz_matching_samples <- as_tibble(get_mz_pairs_share_samples(matrix_cluster_peak_samples_count, 
                                                           cluster_peak$msclusterID))
      if (nrow(pairs_mz_matching_samples) > 0) 
      {
        # first check if the mz pairs are already in the join list, 
        # and skip the ones that will already be joined
        matching_pair_in_join_list <- ((pairs_mz_matching_samples$mz_ID_from %in% join_from_msclusterIDs) & 
                                    (pairs_mz_matching_samples$mz_ID_to %in% join_to_msclusterIDs))
        if (any(matching_pair_in_join_list)) {
          matching_pair_idx_repeated <- which(matching_pair_in_join_list)
          not_new_matching_pair_samples <- sapply(matching_pair_idx_repeated, function(j) {
            (pairs_mz_matching_samples$mz_ID_to[j] %in% join_to_msclusterIDs[join_from_msclusterIDs==pairs_mz_matching_samples$mz_ID_from[j]])
          })
          pairs_mz_matching_samples <- pairs_mz_matching_samples[-matching_pair_idx_repeated[not_new_matching_pair_samples],]
        }
        if (nrow(pairs_mz_matching_samples) > 0) 
        {
        
          # pairs_mz_matching_samples contains all mz pairs that appear in at least one common sample
          # from pairs_mz_matching_samples, check if those mz pairs have matching peakIds, 
          # in the matching samples
          pairs_mz_matching_samples$samples_names <- mapply(function(a, b) {
              # Find column indices where both Row A and Row B are equal to 1
              sub("_spectra$", "", 
                  names(which(matrix_cluster_peak_samples_count[a, ] >= 1 & 
                                matrix_cluster_peak_samples_count[b, ] >= 1)))
            }, match(pairs_mz_matching_samples$mz_ID_from, cluster_peak$msclusterID), 
            match(pairs_mz_matching_samples$mz_ID_to, cluster_peak$msclusterID), 
            SIMPLIFY = FALSE)
          # create customized grepl to check if the pairs of mz sharing a sample also share a peakID, 
          # removing fake ids - which is not consistent to be used for a join
          pairs_mz_matching_samples$valid_pairs <- apply(pairs_mz_matching_samples, 1, function(x) {
            if (length(x[[3]]) == 0)
              return(FALSE)
            peakIds_from <- cluster_peak[cluster_peak$msclusterID == x[[1]], "peakIds"][[1]]
            peakIds_to <- cluster_peak[cluster_peak$msclusterID == x[[2]], "peakIds"][[1]]
            pattern_matching_samples <- paste0(paste0("(", "[0-9]+_", x[[3]][[1]],")"),collapse="|")
            peakIds_matches_from <- stri_extract_all_regex(peakIds_from, pattern_matching_samples)
            peakIds_matches_to <- stri_extract_all_regex(peakIds_to, pattern_matching_samples)
            # if nothing to extract, remove resulting NA to prevent a match
            peakIds_matches_from <- peakIds_matches_from[!is.na(peakIds_matches_from)]
            peakIds_matches_to <- peakIds_matches_to[!is.na(peakIds_matches_to)]
            return(!all(is.na(match(peakIds_matches_from, peakIds_matches_to))))
          })
          
          if (any(pairs_mz_matching_samples$valid_pairs)) {
            # concat the valid pairs of fragmented cluster with matching samples
            # to the sim_cluster_peak
            join_from_msclusterIDs <- c(join_from_msclusterIDs,
                                        pairs_mz_matching_samples$mz_ID_from[pairs_mz_matching_samples$valid_pairs])
            join_to_msclusterIDs <- c(join_to_msclusterIDs, 
                                      pairs_mz_matching_samples$mz_ID_to[pairs_mz_matching_samples$valid_pairs])
          } 
        }
      }
      
      # finally compute the clusters from all mz pairs, and for the ones
      # with more than 2 nodes, merge the lines 
      # if there is at least one pair of mz to join, get unique clusters from them
      if (length(join_from_msclusterIDs) > 0) {
        # convert the ids of the pairs to idx from 0 to the size of the cluster_peak-1
        new_clusters_members <- get_clusters_from_pairs(total_nodes=nrow(cluster_peak),
                                from=match(join_from_msclusterIDs, cluster_peak$msclusterID)-1,
                                to=match(join_to_msclusterIDs, cluster_peak$msclusterID)-1)
        new_clusters_members <- lapply(new_clusters_members, function(x) {
          cluster_peak$msclusterID[x+1]
        })
        new_clusters_size <- sapply(new_clusters_members, length)
      } else {
        # no joins, the clusters are the original cluster_peak, size equals one
        new_clusters_members <- as.list(cluster_peak$msclusterID)
        new_clusters_size <- rep_len(1, length(new_clusters_members))
      }
    }
    # check if there is a new cluster with size > 1 to join; otherwise go to the next m/z
    if (all(new_clusters_size == 1)) # no mz to join
    {
      next()
    } 
    
    # filter the new clusters with size > 1 - at least one join
    new_clusters_members <- new_clusters_members[new_clusters_size > 1]
    new_clusters_size <- new_clusters_size[new_clusters_size > 1]
    for (j in seq_along(new_clusters_members))
    {
      # update number of joins, remove the representative mz that will be kept
      num_joins <- num_joins + new_clusters_size[j] - 1
      new_cluster_j_members <- new_clusters_members[[j]]
      # get the new cluster members pos in the count table
      new_cluster_j_count_idx <- match(new_cluster_j_members, ms_spectra_count$msclusterID)
      # merge counts based on spectra counts
      new_cluster_j <- lapply(names(cluster_peak), merge_counts, cluster_peak[cluster_peak$msclusterID %in% new_cluster_j_members,])
      # get representative cluster ID - which for now is the lower msclusterID idx
      new_cluster_j_representative_ID <- new_cluster_j[[1]]
      new_cluster_j_representative_idx <- match(new_cluster_j_representative_ID, ms_spectra_count$msclusterID)
      # replace representative cluster with the new cluster counts
      ms_spectra_count[new_cluster_j_representative_idx,] <- new_cluster_j
      # remove merged mz rows from the count tables, keeping the representative mz row
      new_cluster_j_count_idx <- new_cluster_j_count_idx[new_cluster_j_count_idx != new_cluster_j_representative_idx]
      ms_spectra_count <- ms_spectra_count[-new_cluster_j_count_idx,] 
      # removed merged mz rows from the count of joins 
      num_joins_last_step <- num_joins_last_step[-new_cluster_j_count_idx]
      total_num_join_clusters <- total_num_join_clusters[-new_cluster_j_count_idx]
    }
    
    # check if the final number of rows in the ms_spectra_count match with
    # the original number of scans minus the number of joins
    if ((nscans - num_joins) != nrow(ms_spectra_count)) {
      stop("ERROR: wrong number of joins compared to the current size of the count table being cleaned.")
    }
  }
  cat("|\n")
  
  # reset number of joins by m/z with the joined clusters of last step 
  ms_spectra_count$numJoins <- ms_spectra_count$numJoins - num_joins_last_step
  
  # reset number of joined clusters by m/z in this step
  num_joins_last_step <- ms_spectra_count$numJoins
  # update number of total joins by cluster and in total after this step
  total_num_join_clusters <- total_num_join_clusters + num_joins_last_step
  num_joins_total <- num_joins_total + num_joins
  
  t2 <- Sys.time()
  cat("\n    * Joined", num_joins, "similar clusters in", 
      round(t2-t1, 2), units(t2-t1), "*\n")
  
  # aggregate similarity values of joined idxs
  if (num_joins > 0)
  {
    # recompute fragmented clusters for the mz with numJoins > 0
    ms_spectra_count[,c("fragmented_clusters","origin_cluster")] <- 
      compute_fragmented_clusters(ms_spectra_count, 
                                rt_tol, mz_tol, 
                                putative_origin=TRUE, 
                                table_type=paste0("m/z with joins from Round ", step_join),
                                with_numJoins=TRUE,
                                add_space_verbose="    ")
    
    # retrieve the new origin_clusters with joins
    new_origin_clusters <- unique(ms_spectra_count[ms_spectra_count$fragmented_clusters!=0 & 
                                                   ms_spectra_count$numJoins > 0, "origin_cluster"][[1]])
    if (length(new_origin_clusters) == 0) {
     # finish cleaning
      # no more joins to make, there is no new fragmented cluster from a 
      # previous joined result (which may have changed)
      break()
    }
    cat("      * Aggregating similarity table with the joined Ids *\n")
    t1 <- Sys.time()
    # for the mz involved in the new fragmented clusters with joins,
    # updated their similarities - replace old msclusterID present in the 
    # joinedIDs with the new msclusterID
    
    # filter clusters to update
    clusters_to_agg_sim <- ms_spectra_count[(
      (ms_spectra_count$origin_cluster %in% new_origin_clusters) & 
       ms_spectra_count$numJoins > 0), c("msclusterID", "numJoins", "joinedIDs")]
    # call aggregate similarities - replace the joinedIDs by the representative msclusterID
    # duplicates and selfloops are removed - overwrite sim_file 
    clusters_sub_status <- sapply(seq_len(nrow(clusters_to_agg_sim)), 
     function(i) {
       # first remove the reference msclusterID from the joinedIDs list
       joinedIDs_sub <- sub(x=clusters_to_agg_sim$joinedIDs[[i]], 
           pattern=paste0(clusters_to_agg_sim$msclusterID[[i]],";"), 
           replacement="")
       # then change ; to |
       joinedIDs_sub <- gsub(x=joinedIDs_sub, pattern=";", replacement="|")
       sim_table_replace_msclusterIDs(list_msclusterIDs_sub = joinedIDs_sub,
                                      representative_msclusterID = clusters_to_agg_sim$msclusterID[[i]],
                                      sim_file)
     })
    nscans <- nrow(ms_spectra_count)
    t2 <- Sys.time()
    cat("        * Done aggregating similarity table in ", 
        round(t2-t1, 2), units(t2-t1), "*\n")
  }
  
  # stop joining if no join was made in the last step and this is not the blanks step (step_join == 0)
  if (step_join >= 1 && (num_joins == 0 || step_join == 10))
    break()
  
  step_join <- step_join + 1
}

# remove last fragmented cluster and origin cluster
ms_spectra_count[,c("fragmented_clusters","origin_cluster")] <- NULL

# get total number of joins by cluster
ms_spectra_count$numJoins <- total_num_join_clusters

# scale the spectra peaksInt by dividing by the sumInts
ms_spectra_count$peaksInt <- sapply(seq_len(nrow(ms_spectra_count)), function(i){
  peaksInt <- as.numeric(strsplit(ms_spectra_count$peaksInt[[i]], ";")[[1]])
  peaksInt <- peaksInt/ms_spectra_count$sumInts[[i]]
  # rescale the peaksInt of the not joined clusters
  if (is.na(ms_spectra_count$joinedIDs[[i]])) {
    peaksInt <- scaleInts(peaksInt, scale_factor)
  }
  paste0(peaksInt, collapse=";")
})

tf <- Sys.time()
cat("\n  * Done reducing", num_joins_total, "similar clusters in", 
    round(tf-t0, 2), units(tf-t0), "*\n\n")
rm(cluster, cluster_peak, total_num_join_clusters, num_joins, num_joins_total, 
   num_joins_last_step)
# warnings()

cat("\n** Computing Peak Area by Sample ** \n\n")
t1 <- Sys.time()
ms_area_count <- ms_spectra_count
names(ms_area_count)[count_columns] <- sub(pattern = "_spectra",
                                           replacement = "_area", fixed = TRUE,
                                           x = names(ms_area_count)[count_columns])
if (!joining_jobs) 
{
  # this is the default flow of np3, compute the peak areas using the provided
  # pre processed data path
  peak_areas_base_peak_int <- compute_peak_area(processed_data_path,
                                                ms_area_count$msclusterID,
                                                lapply(ms_area_count$scans, function(x) strsplit(x, ";")[[1]]),
                                                lapply(ms_area_count$peakIds, function(x) strsplit(x, ";")[[1]]),
                                                batch_metadata)
} else {
  # this is the join_jobs flow, compute the peak areas using all the original samples
  # pre processed data
  peak_areas_base_peak_int <- compute_peak_areas_joined_jobs(output_path, 
                                                             ms_area_count$msclusterID,
                                                             ms_area_count$scans,
                                                             ms_area_count$peakIds, 
                                                             ms_area_count$joinedOriginJobsID)
}
# assign the peak areas following the count columns order 
# and also add the base peak intensity (last column)
ms_area_count[,count_columns] <- peak_areas_base_peak_int[,
                                                          match(names(ms_area_count)[count_columns],
                                                                names(peak_areas_base_peak_int)[-ncol(peak_areas_base_peak_int)])]
ms_spectra_count$basePeakInt <- ms_area_count$basePeakInt <- peak_areas_base_peak_int$basePeakInt
rm(peak_areas_base_peak_int)

# compute max area and the mean precursor intensity of the final clusters
ms_area_count$maxArea <- ms_spectra_count$maxArea <- apply(ms_area_count[,count_columns], 1, max)
ms_area_count$meanInt <- ms_spectra_count$meanInt <- ms_area_count$sumInts / ms_area_count$numSpectra
tf <- Sys.time()
cat("\n  * Done in", round(tf-t1, 2), units(tf-t1), "*\n\n")

cat("\n** Checking joined ids, recomputing samples types indicators and aggregating peak list of joined ids **\n")
t1 <- Sys.time()
# check if joined ID are ok: min ID is in msclusterID column and the others IDs are not
# and aggregate peaksList and peaksInt by ordering and joining adjacent peaks
joined_idx <- which(!is.na(ms_area_count$joinedIDs))
ms_spectra_count[joined_idx, c("peaksList", "peaksInt")] <- 
  ms_area_count[joined_idx, c("peaksList", "peaksInt")] <- 
  Reduce(rbind, lapply(joined_idx, function(i) {
    # check joined scans are ok  
    scan_nums <-  as.numeric(strsplit(ms_area_count$joinedIDs[[i]], split = ";")[[1]])
    min_scan <- min(scan_nums)
    scan_nums <- scan_nums[scan_nums != min_scan]
    
    if (!(min_scan %in% ms_area_count$msclusterID))
      stop("Error in the cleanning step. Missing scan number in the resulting counts table.")
    if (any(scan_nums %in% ms_area_count$msclusterID))
      stop("Error in the cleanning step. Cleaned scan number is still present in the counts table.")
    
    # join adj peaks
    peaks <- as.numeric(strsplit(ms_area_count$peaksList[[i]], split = ";")[[1]])
    ints <- as.numeric(strsplit(ms_area_count$peaksInt[[i]], split = ";")[[1]])
    # order peaks
    peaksOrder <- order(peaks)
    peaks <- peaks[peaksOrder]
    ints <- ints[peaksOrder]
    
    sapply(joinAdjPeaksScalee(peaks, ints, bin_size, -1, scale_factor, 0), 
           function(x) {
             x <- round(x,5)
             paste(x,collapse = ";")
           })
  }), init = NULL)
# check if the inverse scaled intensities sum 1000
check_ints_sum <- sapply(ms_area_count$peaksInt[joined_idx], function(x) {
  ints <- as.numeric(strsplit(x, ";")[[1]])
  # remove scale from intensities
  ints <- inverse_scale(ints, scale_factor)
  round(sum(ints)) == 1000
})
if (!all(check_ints_sum)) {
  stop("Bad scaling of the peak lists' intensities. Some spectra do not have the inverse scales intensities summing 1000 (the normalized value).")
}

# sort by msclsuterID
ms_area_count <- arrange(ms_area_count, msclusterID)
ms_spectra_count <- arrange(ms_spectra_count, msclusterID)

# compute Blanks total and controls total and beds total
blanks_code <- batch_metadata[batch_metadata$SAMPLE_TYPE == "blank", "SAMPLE_CODE"]
if (length(blanks_code) > 0)
{
  if (length(blanks_code) == 1) {
    ms_area_count$BLANKS_TOTAL <- ms_area_count[,  paste0(blanks_code,"_area")][[1]]
    ms_spectra_count$BLANKS_TOTAL <- ms_spectra_count[,  paste0(blanks_code,"_spectra")][[1]]
  } else {
    ms_area_count$BLANKS_TOTAL <- rowSums(ms_area_count[, paste0(blanks_code,"_area")])
    ms_spectra_count$BLANKS_TOTAL <- rowSums(ms_spectra_count[, paste0(blanks_code,"_spectra")])
  }
}
controls_code <- batch_metadata[batch_metadata$SAMPLE_TYPE == "control", "SAMPLE_CODE"]
if (length(controls_code) > 0)
{
  if (length(controls_code) == 1)
  {
    ms_area_count$CONTROLS_TOTAL <- ms_area_count[, paste0(controls_code,"_area")][[1]]
    ms_spectra_count$CONTROLS_TOTAL <- ms_spectra_count[, paste0(controls_code,"_spectra")][[1]]
  } else {
    ms_area_count$CONTROLS_TOTAL <- rowSums(ms_area_count[, paste0(controls_code,"_area")])
    ms_spectra_count$CONTROLS_TOTAL <- rowSums(ms_spectra_count[, paste0(controls_code,"_spectra")])
  }
  rm(controls_code)
}
bed_controls_code <- batch_metadata[batch_metadata$SAMPLE_TYPE == "bed", "SAMPLE_CODE"]

if (length(bed_controls_code) > 0)
{
  if (length(bed_controls_code) == 1)
  {
    ms_area_count$BEDS_TOTAL <- ms_area_count[, paste0(bed_controls_code,"_area")][[1]]
    ms_spectra_count$BEDS_TOTAL <- ms_spectra_count[, paste0(bed_controls_code,"_spectra")][[1]]
  } else {
    ms_area_count$BEDS_TOTAL <- rowSums(ms_area_count[, paste0(bed_controls_code,"_area")])
    ms_spectra_count$BEDS_TOTAL <- rowSums(ms_spectra_count[, paste0(bed_controls_code,"_spectra")])
  }
  rm(bed_controls_code)
}

# remove duplicated entries due to merge
if ("DESREPLICATION" %in% names(ms_area_count))
{
  hflag_rows <- !is.na(ms_area_count$HFLAG) | !is.na(ms_area_count$DESREPLICATION)
  ms_spectra_count[hflag_rows, c("DESREPLICATION", "HFLAG")] <- 
  ms_area_count[hflag_rows, c("DESREPLICATION", "HFLAG")] <- Reduce(rbind, 
    lapply(which(hflag_rows), 
           function(i)
           {
             if (is.na(ms_area_count$HFLAG[[i]]))
               hflag <- NA
             else
               hflag <- paste0(unique(strsplit(ms_area_count$HFLAG[[i]], ";")[[1]]), collapse = ";")
             if (is.na(ms_area_count$DESREPLICATION[[i]]))
               desrep <- NA
             else
               desrep <- paste0(unique(strsplit(ms_area_count$DESREPLICATION[[i]], ";")[[1]]), collapse = ";")
             
             return(c(desrep, hflag))
           }), init = NULL)
}
tf <- Sys.time()
cat("\n  * Done in", round(tf-t1, 2), units(tf-t1), "*\n\n")

# if join_jobs, create the cleanClustID column here to store the current cleaned msclusterIDs
# and apply the integrative clustering heuristic to maintain the msclusterIDs from the
# first job selected as reference, using the msclusterID_integrative column created in step 4
# to extract the minimum original msclusterID; here the msclusterID_integrative
# will contain the original msclusterIDs that got joined in the clustering steps concatenated by ;
if ("msclusterID_integrative" %in% names(ms_spectra_count)) {
  cat("\n** Computing the integrative msclusterID **\n\n")
  t1 <- Sys.time()
  ms_spectra_count$cleanClustID <- ms_area_count$cleanClustID <- ms_spectra_count$msclusterID
  # convert to char
  ms_spectra_count$msclusterID_integrative <- ms_area_count$msclusterID_integrative <- as.character(ms_spectra_count$msclusterID_integrative)
  new_msclusterIDs <- sapply(ms_spectra_count$msclusterID_integrative, function(x) {
    min(as.integer(strsplit(x, ";")[[1]]))
  })
  ms_spectra_count$msclusterID <- ms_area_count$msclusterID <- new_msclusterIDs
  tf <- Sys.time()
  cat("\n  * Done in", round(tf-t1, 2), units(tf-t1), "*\n\n")
}

cat("\n** Saving the clean MGF file **\n\n")
t1 <- Sys.time()
# sort by msclusterID, this order will be applied to the scans in the MGF file
ms_area_count <- arrange(ms_area_count, msclusterID)
ms_spectra_count <- arrange(ms_spectra_count, msclusterID)
# write the peaks list and ints to the mgf file - the intensities will be 
# inversed scaled to be saved with no scaling
writeMgfDataFile_NP3_table(ms_count_table = ms_area_count,
                           file_MGF = file.path(output_path, "mgf",
                                                paste0(output_name, "_clean.mgf")),
                           output_name = output_name,
                           charge = ion_mode,
                           scale_factor = scale_factor)
tf <- Sys.time()
cat("\n  * Done in", round(tf-t1, 2), units(tf-t1), "*\n\n")

# round mzConsensus and rts
ms_spectra_count$mzConsensus <- ms_area_count$mzConsensus <- round(ms_area_count$mzConsensus, mz_rt_digits)
ms_spectra_count$rtMean <- ms_area_count$rtMean <- round(ms_area_count$rtMean, mz_rt_digits)
ms_spectra_count$rtMin <- ms_area_count$rtMin <- round(ms_area_count$rtMin, mz_rt_digits)
ms_spectra_count$rtMax <- ms_area_count$rtMax <- round(ms_area_count$rtMax, mz_rt_digits)

# create the retention time columns in minutes
ms_spectra_count$rtMean_minutes <- ms_area_count$rtMean_minutes <- round(ms_area_count$rtMean/60, mz_rt_digits)
ms_spectra_count$rtMin_minutes <- ms_area_count$rtMin_minutes <- round(ms_area_count$rtMin/60, mz_rt_digits)
ms_spectra_count$rtMax_minutes <- ms_area_count$rtMax_minutes <- round(ms_area_count$rtMax/60, mz_rt_digits)

# compute number of fragmented clusters
ms_spectra_count$fragmented_clusters <- ms_area_count$fragmented_clusters <- compute_fragmented_clusters(ms_area_count, rt_tol, mz_tol)

# write cleanned data without annotation
write_csv(ms_area_count, path = file.path(output_path, "count_tables", "clean", 
                                          paste0(output_name, "_peak_area_clean.csv")))
write_csv(ms_spectra_count, path = file.path(output_path, "count_tables", "clean", 
                                             paste0(output_name, "_spectra_clean.csv")))
rm(ms_spectra_count, ms_area_count)

tf <- Sys.time()
cat("\n** Done Clean in", round(tf-t0, 2), units(tf-t0), "**\n\n")