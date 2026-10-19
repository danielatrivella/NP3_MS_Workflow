library(data.table)

# compute the number of fragmented clusters present in the provided count table and
# return a column with the putative fragmented clusters of each m/z by row (same order of the provided table)
# assigned as follow:
# 0 for a not fragmented cluster m/z; 1 or greater for a m/z that have these number of fragmented clusters;
# and -1 for a m/z signalized as a fragmented cluster
# where
# @putative_origin indicates if the column with the msclusterID of the m/z
# assigned as the putative origin cluster of the detected fragmented clusters should be returned
# if FALSE only returns the vector with the putative fragmented clusters assigned as described above
# @with_num_joins indicates to search for fragmented clusters only in the rows where num_joins > 0 (used in the clean step)
compute_fragmented_clusters <- function(clean_count_table,
                                        rt_tol = 2, mz_tol = 0.025,
                                        putative_origin = FALSE,
                                        table_type = "Step 5 clean",
                                        with_numJoins = FALSE,
                                        add_space_verbose = "") {
  cat("\n", add_space_verbose, "* Computing the number of fragmented cluster in the", table_type, "table *\n")
  t1 <- Sys.time()
  
  n <- nrow(clean_count_table)
  
  # extract the needed columns once (no copy: R shares the memory until modified)
  mz        <- clean_count_table[["mzConsensus"]]
  rtMean    <- clean_count_table[["rtMean"]]
  rtMin     <- clean_count_table[["rtMin"]]
  rtMax     <- clean_count_table[["rtMax"]]
  basePeakInt  <- clean_count_table[["basePeakInt"]]
  msclusterID      <- clean_count_table[["msclusterID"]]
  blanks_flag <- "BLANKS_TOTAL" %in% names(clean_count_table)
  bflag     <- if (blanks_flag) clean_count_table[["BFLAG"]] else NULL
  numJoins  <- if (with_numJoins) clean_count_table[["numJoins"]] else NULL
  
  frag   <- numeric(n)                 # fragmented_clusters (0 or >0 or -1)
  origin <- msclusterID[rep(NA_integer_, n)]  # origin_cluster: NA of the same type as msclusterID
  
  # ---- m/z index: sort once, find every row's window by binary search ----
  # data.table's fast radix ordering (NA last)
  ord <- data.table:::forder(mz) 
  #ord <- ord[!is.na(mz[ord])] # it should not have any NA
  mz_sorted <- mz[ord]
  eps <- 1e-9         # pad the window, exact |dmz| <= mz_tol test is done below
  lo <- findInterval(mz - mz_tol - eps, mz_sorted, left.open = TRUE) + 1L
  hi <- findInterval(mz + mz_tol + eps, mz_sorted)
  
  for (i in seq_len(n)) {
    if (frag[i] < 0) next
    if (with_numJoins && numJoins[i] == 0) next
    if (is.na(lo[i]) || lo[i] > hi[i]) next
    
    # clusters candidates in the m/z window
    cluster_peak_cand_idx <- ord[lo[i]:hi[i]]
    rt_min <- rtMin[i]; rt_max <- rtMax[i]; rt_mean <- rtMean[i]
    
    # get the cluster peak, with putative fragmented clusters
    keep_cluster_peak <- abs(mz[cluster_peak_cand_idx] - mz[i]) <= mz_tol &
      ((rtMean[cluster_peak_cand_idx] >= rt_min - rt_tol & rtMean[cluster_peak_cand_idx] <= rt_max + rt_tol) |
         (rt_mean >= rtMin[cluster_peak_cand_idx] - rt_tol & rt_mean <= rtMax[cluster_peak_cand_idx] + rt_tol))
    cluster_peak_cand_idx <- cluster_peak_cand_idx[which(keep_cluster_peak)]
    # if no fragmented cluster detected, go to next cluster
    if (length(cluster_peak_cand_idx) < 2L) next
    
    # keep the original table order (needed for ties in which.max and identical results)
    cluster_peak_cand_idx <- sort.int(cluster_peak_cand_idx)
    
    # if not bflag, check peak center and boundaries deviation, remove peaks not aligned
    # the spectra peak center deviation is <= 4*rt_tol or peak boundaries deviation <= 2*rt_tol
    if (!blanks_flag || !bflag[i]) {
      keep_cluster_peak <- abs(rt_mean - rtMean[cluster_peak_cand_idx]) <= 4 * rt_tol |
        sqrt(((rtMin[cluster_peak_cand_idx] - rt_min)^2 + (rtMax[cluster_peak_cand_idx] - rt_max)^2) / 2) <= 2 * rt_tol  # RMSE of (rtMin, rtMax)
      cluster_peak_cand_idx <- cluster_peak_cand_idx[which(keep_cluster_peak)]
    }
    
    # there if there is an fragmented cluster that was not merged and 
    # assign the origin of all clusters as the spectra in the current 
    # cluster_peak with the biggest base peak intensity
    if (length(cluster_peak_cand_idx) > 1L) { 
      j <- which.max(basePeakInt[cluster_peak_cand_idx])
      origin_cluster_idx <- cluster_peak_cand_idx[j]
      origin_ID <- msclusterID[origin_cluster_idx]
      frag[origin_cluster_idx]   <- length(cluster_peak_cand_idx) - 1L
      origin[origin_cluster_idx] <- origin_ID
      others_cluster_peak <- cluster_peak_cand_idx[-j]
      frag[others_cluster_peak]   <- -1
      origin[others_cluster_peak] <- origin_ID
    }
  }
  
  number_fragmented_clusters <- sum(frag != 0)
  number_unique_origin_cluster <- length(unique(origin[!is.na(origin)]))
  
  cat("\n", add_space_verbose, 
      "  - Total number of putative fragmented cluster and percentage overall m/zs:",
      number_fragmented_clusters, 
      "(", round(100 * number_fragmented_clusters / n, 2), "% ) \n")
  cat("", add_space_verbose, 
      "  - Total number of putative unique origin clusters and percentage overall fragmented clusters:",
      number_unique_origin_cluster, 
      "(", round(100 * number_unique_origin_cluster / number_fragmented_clusters, 2), "% ) \n")
  
  t2 <- Sys.time()
  cat("", add_space_verbose, "  * Done in", round(t2 - t1, 2), units(t2 - t1), "*\n")
  
  # return the resulting fragmented clusters column to be added to the clean tables
  if (putative_origin) {
    return(data.table(fragmented_clusters = frag, origin_cluster = origin))
  } else {
    return(frag)
  }
}
