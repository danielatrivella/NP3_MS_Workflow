library(readr)

RMSE <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  sqrt(mean((x - y)^2))
}

# compute the number of fragmented clusters present in the provided count table and
# return a column with the putative fragmented clusters of each m/z by row (same order of the provided table) 
# assigned as follow:
# 0 for a not fragmented cluster m/z; 1 or greater for a m/z that have these number of fragmented clusters;
# and -1 for a m/z signalized as a fragmented cluster
# where putative_origin indicates if the column with the msclusterID of the m/z 
# assigned as the putative origin cluster of the detected fragmented clusters should be returned
# if FALSE only returns the vector with the putative fragmented clusters assigned as described above
compute_fragmented_clusters <- function(clean_count_table, rt_tol=2, mz_tol=0.025, putative_origin = FALSE, table_type="Step 5 clean") {
  cat("\n  - Computing the number of fragmented cluster in the",table_type,"table\n")
    
  if ('BLANKS_TOTAL' %in% names(clean_count_table)) {
    blanks_flag <- TRUE
  } else {
    blanks_flag <- FALSE
  }
  clean_count_table$fragmented_clusters <- 0
  clean_count_table$origin_cluster <- NA
  clean_count_table$idx <- 1:nrow(clean_count_table)
  
  for (i in seq_len(nrow(clean_count_table))) {
    # print(i)
    if (clean_count_table$fragmented_clusters[[i]] < 0) {
      next()
    }
    
    cluster <- clean_count_table[i,]
    
    # get the cluster peak, with putative fragmented clusters
    cluster_peak <- clean_count_table[abs(clean_count_table$mzConsensus - cluster$mzConsensus) <= mz_tol &
                                    ((clean_count_table$rtMean >= cluster$rtMin - rt_tol & 
                                        clean_count_table$rtMean <= cluster$rtMax + rt_tol) |
                                       (cluster$rtMean >= clean_count_table$rtMin - rt_tol &
                                          cluster$rtMean <= clean_count_table$rtMax + rt_tol)),] 
    # if not bflag check peak center and boundaries deviation, remove peaks not aligned
    # the spectra peak center deviation is <= 4*rt_tol or peak boundaries deviation <= 2*rt_tol
    if (!blanks_flag || !cluster$BFLAG) {
      cluster_peak <- cluster_peak[abs(cluster$rtMean-cluster_peak$rtMean) <= 4*rt_tol |
                                     apply(cluster_peak[,c("rtMin", "rtMax")], 1,
                                           function(x,y)  RMSE(x, y),
                                           y = c(cluster$rtMin, cluster$rtMax)) <= 2*rt_tol,]
    }
    
    if (nrow(cluster_peak) > 1) # there is an adj cluster that was not merged
    {
      # assign the spectra with the bigger base peak intensity as the principal and the others as fragmented clusters
      j <- which.max(cluster_peak$basePeakInt)
      origin_ID <- cluster_peak$msclusterID[[j]]
      #cat("msclusterID ", origin_ID, " have fragmented clusters.\n")
      clean_count_table$fragmented_clusters[cluster_peak$idx[[j]]] <- nrow(cluster_peak) - 1
      clean_count_table$origin_cluster[cluster_peak$idx[[j]]] <- origin_ID
      cluster_peak <- cluster_peak[-j,]
      clean_count_table$fragmented_clusters[cluster_peak$idx] <- -1
      clean_count_table$origin_cluster[cluster_peak$idx] <- origin_ID
    }
  }
  
  n <- nrow(clean_count_table)
  number_fragmented_clusters <- sum(clean_count_table$fragmented_clusters<0)
  number_mzs_w_fragmented_clusters <- sum(clean_count_table$fragmented_clusters>0)
  cat("\n    - Total number of m/zs that had at least one putative fragmented clusters and percentage overall m/zs:", 
      number_mzs_w_fragmented_clusters, "(",round(100*number_mzs_w_fragmented_clusters/nrow(clean_count_table),2), "% ) \n")
  cat("\n    - Total number of putative fragmented clusters and percentage overall m/zs:", 
      number_fragmented_clusters, "(",round(100*number_fragmented_clusters/nrow(clean_count_table),2), "% ) \n")
  
  # return the resulting fragmented clusters column to be added to the clean tables
  if (putative_origin) {
    return(clean_count_table[,c("fragmented_clusters","origin_cluster")])
  } else {
    return(clean_count_table$fragmented_clusters)
  }
  #write_csv(clean_count_table, sub(".csv","_fragmentedClusters.csv", count_table_path))
}

# TODO make this function with data.table - benchmark time gain
# TODO make this function for only a selection of m/zs (the ones with a join (modified ones))

# TODO try clustering the origin cluster with the fragmented clusters, than 
# if more than 1 fragmented cluster remaining, 
# also try to cluster them between each other