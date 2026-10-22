## ----load-libs

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

suppressPackageStartupMessages(library(data.table))
Rcpp::sourceCpp(file.path(script_path(), 'read_mgf_peak_list_R.cpp'))
Rcpp::sourceCpp(file.path(script_path(), 'dot_product_list.cpp'))


# function to compute similarity value between the members of defined groups in the provided df_mscluster_sets
# @ df_mscluster_sets data frame with the 'msclusterID' column and the grouping column called 'origin_cluster' with the msclusterID representative of each group
# @ path_mgh path to the mgf file to use in the pairwise comparison
# @ sim_output_path path to the similarity table to store the computed pairwise similarity values
# @ max_shift maximum mass difference allowed to search for shifted m/z fragmented ion in the cosine computation
# @ groups_per_batch number of groups processed (in parallel) before the results are collected and 
#   written to the output file; NULL (default) selects it automatically
#
# The parallelization is done over the groups (forking), each group is compared sequentially.
# The groups are processed in batches to keep the intermediary saves of the similarities
# and a low memory usage; the similarity table keeps the order of the groups.
pairwise_similarity_comparison_sets <- function(df_mscluster_sets, path_mgf, sim_output_path, 
                                                set_name="fragmented clusters from clustering Step 3",
                                                sim_cutoff=0.55,
                                                bin_size=0.05, 
                                                scale_factor=0.5,
                                                trim_mz=TRUE,
                                                max_shift=200,
                                                join_isotopic_peaks=1, # always on
                                                parallel_cores=2,
                                                groups_per_batch=NULL) {
  # function to returns a list with the columns (to be bound with rbindlist) with the sim values 
  # >= sim_cutoff from source msclusterID in ms2_sample[[i]] and the targets in ms2_sample[set_j] 
  # where set_j is a vector of idxs
  # the list have 4 columns: "msclusterID_source","msclusterID_target","similarity_value","num_matched_peaks"
  compareSpectraNormDotProductSet <- function(i, set_j, set_size, sim_cutoff)
  {
    i_idx <- set_j[[i]]
    j_idx <- set_j[(i+1):set_size]
    # gets the list of valid cosine >= sim_cutoff, their number of matched peaks, and their idx in [set_j]
    # the cosine is already rounded in 3 decimal places
    # returns a list with three elements: "similarity_value", "num_matched_peaks","idx_target"
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
    # convert the target idx to their msclusterID and return the columns as a list
    # (no matrix is created, the lists are bound with rbindlist)
    # otherwise return NULL
    number_valid_sim <- length(np3_cos_matches_list[[1]])
    if (number_valid_sim > 0) {
      # fix the valid idx to be equal their idx in ms2_sample
      idx_target <- j_idx[np3_cos_matches_list[[3]]]
      return(list(msclusterID_source = rep.int(ms2_sample$SCANS[[i_idx]], number_valid_sim),
                  msclusterID_target = ms2_sample$SCANS[idx_target],
                  similarity_value = np3_cos_matches_list[[1]],
                  num_matched_peaks = np3_cos_matches_list[[2]]))
    } else {
      return(NULL)
    }
  }
  
  # function to compare all the pairs of one group (sequential): g is the group index in
  # group_members, a list with the idxs (in ms2_sample) of the members of each group
  # returns a data.table with the valid similarities of the group, or NULL if there is none
  compareGroupPairs <- function(g, group_members, sim_cutoff)
  {
    set_j <- group_members[[g]]
    set_size <- length(set_j)
    if (set_size < 2) {
      return(NULL)
    }
    group_sims <- rbindlist(lapply(seq_len(set_size-1), 
                                   compareSpectraNormDotProductSet,
                                   set_j=set_j,
                                   set_size=set_size,
                                   sim_cutoff=sim_cutoff))
    if (nrow(group_sims) == 0) {
      return(NULL)
    }
    return(group_sims)
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
  # (only the two needed columns are extracted, the input table is not copied)
  keep_rows <- which(!is.na(df_mscluster_sets$origin_cluster))
  set_msclusterIDs <- df_mscluster_sets$msclusterID[keep_rows]
  set_origin_cluster <- df_mscluster_sets$origin_cluster[keep_rows]
  rm(keep_rows)
  
  # start pairwise comparison between the groups
  cluster_groups <- unique(set_origin_cluster)
  n_groups <- length(cluster_groups)
  # number of members of each group, in the order of cluster_groups
  group_sizes <- tabulate(match(set_origin_cluster, cluster_groups), nbins = n_groups)
  n_comparison_within_groups <- (sum(as.numeric(group_sizes)^2) - length(set_origin_cluster))/2
  
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
                                 scans_to_keep=sort(set_msclusterIDs))
  
  if (!all(ms2_sample$SCANS %in% set_msclusterIDs)) 
  {
    stop("ERROR when reading the MGF '",path_mgf,
         "' with the provided list of o msclsuterIDs to keep. The list of SCANS ",
         "that were read is not contained in the set of provided msclusterIDs.")
    
  }
  # get scan order
  scan_index <- order(ms2_sample$SCANS) # get the indexes sorted by scan
  
  # order mgf data
  ms2_sample <- lapply(ms2_sample, function(x) x[scan_index])
  
  # position of each msclusterID of the set in the (sorted by scan) ms2_sample: 
  # these are the indexes used in the comparisons. When the set is sorted by msclusterID,
  # as in the clean step, they are equal to the row positions in the set
  ms2_pos <- match(set_msclusterIDs, ms2_sample$SCANS)
  if (anyNA(ms2_pos))
  {
    stop("ERROR: ", sum(is.na(ms2_pos)), " msclusterIDs of the provided set were not found ",
         "in the MGF '", path_mgf, "'.")
  }
  # members (idxs in ms2_sample) of each group, in a single group by pass, 
  # keeping the groups order of cluster_groups and the members order of the set
  group_members <- data.table(pos = ms2_pos, origin = set_origin_cluster)[, .(pos = list(pos)), 
                                                                          by = origin]$pos
  rm(ms2_pos, set_msclusterIDs, set_origin_cluster)
  
  # column names of the similarity table, in the order they are written
  sim_columns <- c("msclusterID_source", "msclusterID_target", 
                   "similarity_value", "num_matched_peaks")
  
  # buffer to store the list of comparisons until the next save,
  # and var to store number of comparison written during the processing
  pending_sims <- list()
  n_pending_sims <- 0
  n_stored_comps <- 0
  
  # write the buffered similarities (intermediary saves when reaching more than 20k rows)
  write_pending_sims <- function() {
    sims_dt <- rbindlist(pending_sims)
    fwrite(sims_dt, file = file.path(sim_output_path), 
           append = (n_stored_comps > 0), 
           col.names = (n_stored_comps == 0), 
           showProgress = FALSE)
    n_stored_comps <<- n_stored_comps + nrow(sims_dt)
    pending_sims <<- list()
    n_pending_sims <<- 0
  }
  
  # the groups are processed in batches: the parallelization is done over the groups of the batch 
  # (mclapply, forking) and each group is compared sequentially; the results of each 
  # batch are collected, checked and buffered before the next one
  if (is.null(groups_per_batch)) {
    groups_per_batch <- max(parallel_cores * 500, ceiling(n_groups / 50))
  }
  groups_per_batch <- max(1, as.integer(groups_per_batch))
  batches <- split(seq_len(n_groups), ceiling(seq_len(n_groups) / groups_per_batch))
  
  for (batch in batches) {
    # process the biggest groups first to balance the work between the cores 
    # (the work grows with the squared size of the group), results are put back in the groups order
    process_order <- order(group_sizes[batch], decreasing = TRUE)
    if (parallel_cores > 1)
    {
      # run parallelized version over the groups, use forking
      batch_results <- mclapply(batch[process_order], 
                                compareGroupPairs, 
                                group_members=group_members,
                                sim_cutoff=sim_cutoff,
                                mc.cores = parallel_cores)
    } else {
      # sequential pairwise comparisons
      batch_results <- lapply(batch[process_order], 
                              compareGroupPairs, 
                              group_members=group_members,
                              sim_cutoff=sim_cutoff)
    }
    batch_results <- batch_results[order(process_order)]
    
    if (any(vapply(batch_results, inherits, logical(1), what = "try-error"))) {
      stop("ERROR: a failure happened in the pairwise comparison of the groups ",
           min(batch), " to ", max(batch), " for mgf '", path_mgf, "'.")
    }
    
    # progress, the marks of the groups of this batch
    n_marks <- sum(progress_comp <= max(batch))
    if (n_marks > 0) {
      progress_comp <- progress_comp[progress_comp > max(batch)]
      cat(strrep("=", n_marks))
    }
    
    # if at least one valid sim, buffer it to concatenate with previous similarities
    # (rbindlist ignores the NULL results of the groups without valid similarities)
    batch_sims <- rbindlist(batch_results)
    rm(batch_results)
    if (nrow(batch_sims) > 0) {
      # check if all target ID is valid
      if (anyNA(batch_sims$msclusterID_target))
      {
        stop("ERROR: a wrong index mapping was detected in the similarity matrix for mgf '",
             path_mgf,"'. ")
      }
      pending_sims[[length(pending_sims) + 1]] <- batch_sims
      n_pending_sims <- n_pending_sims + nrow(batch_sims)
      # make intermediary saves when reaching more than 50k rows
      if (n_pending_sims > 50000) {
        write_pending_sims()
      }
    }
    rm(batch_sims)
  }
  
  if (n_pending_sims > 0) {
    # save cosine values
    write_pending_sims()
  } else if (n_stored_comps == 0) {
    # if no sim, save empty sim table (only the header)
    writeLines(paste(sim_columns, collapse = ","), con = file.path(sim_output_path))
  }
  
  cat("|\n")
  tf <- Sys.time()
  cat("\n    - Total number of valid similarities stored = ", n_stored_comps,"\n")
  cat("\n    * Done making", n_comparison_within_groups,
      "pairwise comparisions in", round(tf-ti, 2), units(tf-ti), "*\n")
  
  return(n_stored_comps)
}
