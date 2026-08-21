## ----load-libs
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


# function to compute similarity value between the members of defined groups in the provided df_mscluster_sets
# @ df_mscluster_sets data frame with the 'msclusterID' column and the grouping column called 'origin_cluster' with the msclusterID representative of each group
# @ path_mgh path to the mgf file to use in the pairwise comparison
# @ sim_output_path path to the similarity table to store the computed pairwise similarity values
# @ max_shift maximum mass difference allowed to search for shifted m/z fragmented ion in the cosine computation
pairwise_similarity_comparison_sets <- function(df_mscluster_sets, path_mgf, sim_output_path, 
                                                set_name="fragmented clusters from clustering Step 3",
                                                sim_cutoff=0.55,
                                                bin_size=0.05, 
                                                scale_factor=0.5,
                                                trim_mz=TRUE,
                                                max_shift=200,
                                                join_isotopic_peaks=1, # always on
                                                parallel_cores=2) {
  # function to returns a matrix with the sim values >= 0.1 from source msclusterID in
  # ms2_sample[[i]] and the targets in ms2_sample[set_j] where set_j is a vector of idxs
  # the matrix have 4 columns: "similarity_value","num_matched_peaks","msclusterID_source","msclusterID_target"
  compareSpectraNormDotProductSet <- function(i, set_j, set_size, sim_cutoff)
  {
    i_idx <- set_j[[i]]
    j_idx <- set_j[(i+1):set_size]
    # gets the list of valid cosine >= sim_cutoff, their number of matched peaks, and their idx in [set_j]
    # the cosine is already rounded in 3 decimal places
    np3_cos_matches_list <- normDotProductShiftList(peaks_A = ms2_sample$MZS[[i_idx]], 
                                                    ints_A = ms2_sample$INTS[[i_idx]], 
                                                    mz_A = ms2_sample$PREC_MZ[[i_idx]],
                                                    peaks_B = ms2_sample$MZS[j_idx],
                                                    ints_B = ms2_sample$INTS[j_idx], 
                                                    mzs_B = ms2_sample$PREC_MZ[j_idx],
                                                    bin_size = bin_size,
                                                    max_shift = max_shift,
                                                    sim_cutoff = sim_cutoff)
    # if at least one valid sim, also set source msclusterID - equals SCANS[[i]] - 
    # convert the target idx to their msclusterID and return as a matrix
    # otherwise return NULL
    number_valid_sim <- length(np3_cos_matches_list[[1]])
    if (number_valid_sim > 0) {
      np3_cos_matches_list <- matrix(unlist(np3_cos_matches_list), 
                                     nrow = number_valid_sim, ncol=3, byrow = FALSE)
      dimnames(np3_cos_matches_list) <- list(NULL, c("similarity_value", "num_matched_peaks","idx_target"))
      # fix the valid idx to be equal their idx in set_j
      np3_cos_matches_list[,"idx_target"] <- j_idx[np3_cos_matches_list[,"idx_target"]]
      np3_cos_matches_list <- cbind(np3_cos_matches_list, msclusterID_source = ms2_sample$SCANS[[i_idx]],
                                    msclusterID_target = ms2_sample$SCANS[np3_cos_matches_list[,"idx_target"]])
      return(np3_cos_matches_list[,c("similarity_value","num_matched_peaks","msclusterID_source","msclusterID_target")])
    } else {
      return(NULL)
    }
  }
  
  if (!all(c("msclusterID", "origin_cluster") %in% names(df_mscluster_sets))) 
  {
    stop("The provided data frame does not contain the expected columns with ",
         "the set of msclusters to be compared: 'msclusterID' and 'origin_cluster'.")
  }
 
  path_mgf <- normalizePath(path_mgf)
  if (!file.exists(path_mgf))
  {
    stop("The MGF file '", path_mgf,
         "' does not exists. Provide a valid path to the MGF file is located.")
  }
  
  if (!is.numeric(sim_cutoff))
    stop("The similarity cutoff must be a numeric value. Wrong value informed: ",sim_cutoff)
  
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
  
  # first remove any NA in the origin groups - which are the mz without a group
  df_mscluster_sets <- df_mscluster_sets[!is.na(df_mscluster_sets$origin_cluster),]

  # start pairwise comparison between the groups
  cluster_groups <- unique(df_mscluster_sets$origin_cluster)
  n_groups <- length(cluster_groups)
  n_comparison_within_groups <- (sum(table(mz_isomers_groups$origin_cluster)^2) - nrow(df_mscluster_sets))/2

  cat("\n  * Comparing", set_name, "pairwise by groups *\n")
  if (n_groups == 0) {
    cat("\n  - No comparison to make, there is no group in the provided set of putative clusters. Aborting similarity computation... ")
    return(0)
  }
  if (parallel_cores > 1 && !require(parallel)) {
    parallel_cores <- 1
    cat("    - Parallelization disabled. Could not load parallel package. ",
        "Serial pairwise comparison will be used.\n\n")
  }  else {
    cat("    - Number of parallel cores to use in the pairwise comparison:",parallel_cores,"\n\n")
  }
  
  # make progress have max length like in the clean step
  # add progress
  progress_comp <- unique(trunc(c(seq(from = 1, to = n_groups, by = n_groups/25), n_groups)))
  n_progress <- length(progress_comp)
  cat(paste0("        |", paste0(rep(" ", n_progress), collapse = ""), 
             "|\n        |", collapse = ""))
  
  ti <- Sys.time()
  
  # read the ms2 data and get the scan indexes
  ms2_sample <- readMgfPeaksList(path_mgf, bin_size = bin_size, 
                                 trim_mz = trim_mz, scale_factor = scale_factor,
                                 join_isotopic_peaks = join_isotopic_peaks,
                                 scans_to_keep=sort(df_mscluster_sets$msclusterID))
  
  if (!all(ms2_sample$SCANS %in% df_mscluster_sets$msclusterID)) 
  {
    stop("ERROR when reading the MGF '",path_mgf,
         "' with the provided list of o msclsuterIDs to keep. The list of SCANS ",
         "that were read is not contained in the set of provided msclusterIDs.")
    
  }
  # get scan order
  scan_index <- order(ms2_sample$SCANS) # get the indexes sorted by scan
  n_scans <- length(ms2_sample$SCANS)
  
  # order mgf data
  ms2_sample <- lapply(ms2_sample, function(x) x[scan_index])
  scan_index <- seq_along(ms2_sample$SCANS)
  scans_num <- ms2_sample$SCANS
  
  # create empty tibble to store the list of comparisons
  comp_groups_sim_matches <- dplyr::tibble(msclusterID_source=numeric(0), 
                                           msclusterID_target=numeric(0), 
                                           similarity_value=numeric(0), 
                                           num_matched_peaks=numeric(0))
  # create var to store number of comparison written during the processing
  n_stored_comps <- 0
  
  for (i in seq_along(cluster_groups)) {
    if (n_progress > 0 && i == progress_comp[[1]]) {
      progress_comp <- progress_comp[-1]
      n_progress <- n_progress - 1
      cat("=")
    }
    # retrieve origin cluster ID and members idx and size
    origin_cluster_ID <- cluster_groups[[i]]
    cluster_group_idxs <- which(df_mscluster_sets$origin_cluster == origin_cluster_ID)
    cluster_group_size <- length(cluster_group_idxs)
    
    if (parallel_cores > 1 && require(parallel))
    {
      # run paralellized verison o lapply, use forking
      comp_row_sim_matches <- dplyr::as_tibble(do.call(rbind, mclapply(seq_len(cluster_group_size-1), 
                                                                       compareSpectraNormDotProductSet, 
                                                                       set_j=cluster_group_idxs,
                                                                       set_size=cluster_group_size,
                                                                       sim_cutoff=sim_cutoff,
                                                                       mc.cores = parallel_cores)))
      #stopCluster(cl)
    } else {
      # sequential pairwise comparisions
      comp_row_sim_matches <- dplyr::as_tibble(do.call(rbind, lapply(seq_len(cluster_group_size-1), 
                                                                     compareSpectraNormDotProductSet,
                                                                     set_j=cluster_group_idxs,
                                                                     set_size=cluster_group_size,
                                                                     sim_cutoff=sim_cutoff,
                                                                     set_size=length(cluster_group_idxs),)))
    }
    
    # if at least one valid sim, concatenate it with previous similarities
    if (nrow(comp_row_sim_matches) > 0) {
      # check if all target ID is valid
      if (any(is.na(comp_row_sim_matches[,"msclusterID_target"])))
      {
        stop("ERROR: a wrong index mapping was detected in the similarity matrix for mgf '",
             path_mgf,"'. ")
      }
      comp_groups_sim_matches <- bind_rows(comp_groups_sim_matches, 
                                           comp_row_sim_matches)
      # make intermediary saves when reaching more than 20k rows
      if (nrow(comp_groups_sim_matches) > 20000) {
        readr::write_csv(x=comp_groups_sim_matches,
                         path = file.path(sim_output_path),
                         append=(n_stored_comps>0), 
                         col_names=(n_stored_comps==0))
        n_stored_comps <- n_stored_comps + nrow(comp_groups_sim_matches)
        # reset df
        comp_groups_sim_matches <- dplyr::tibble(msclusterID_source=numeric(0), 
                                                 msclusterID_target=numeric(0), 
                                                 similarity_value=numeric(0), 
                                                 num_matched_peaks=numeric(0))
      }
    }
    rm(comp_row_sim_matches)
  }
  
  # if no sim, save empty sim table
  if (nrow(comp_groups_sim_matches) == 0) {
    if (n_stored_comps == 0) {
    readr::write_csv(x=dplyr::tibble(msclusterID_source=numeric(0), 
                                     msclusterID_target=numeric(0), 
                                     similarity_value=numeric(0), 
                                     num_matched_peaks=numeric(0)),
                     path = file.path(sim_output_path))
    }
  } else {
    # TODO eval fwrite from data.table 
    # save cosine values
    readr::write_csv(x=comp_groups_sim_matches,
                     path = file.path(sim_output_path),
                     append=(n_stored_comps>0), 
                     col_names=(n_stored_comps==0))
    n_stored_comps <- n_stored_comps + nrow(comp_groups_sim_matches)
  }
  
  cat("|\n")
  tf <- Sys.time()
  cat("\n    - Total number of valid similarities stored = ", n_stored_comps,"\n")
  cat("\n    * Done making", n_comparison_within_groups,
      "pairwise comparisions in", round(tf-ti, 2), units(tf-ti), "*\n")
  
  return(n_stored_comps)
}