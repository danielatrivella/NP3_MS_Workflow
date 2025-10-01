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
compute_fragmented_clusters <- function(clean_count_table, rt_tol=2, mz_tol=0.025) {
  cat("\n  - Computing the number of fragmented cluster in the clean table\n")
    
  if ('BLANKS_TOTAL' %in% names(clean_count_table)) {
    blanks_flag <- TRUE
  } else {
    blanks_flag <- FALSE
  }
  clean_count_table$fragmented_clusters <- 0
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
    # the spectra peak center deviation is <= rt_tol or peak boundaries deviation <= 2*rt_tol
    if (!blanks_flag || !cluster$BFLAG) {
      cluster_peak <- cluster_peak[abs(cluster$rtMean-cluster_peak$rtMean) <= rt_tol |
                                     apply(cluster_peak[,c("rtMin", "rtMax")], 1,
                                           function(x,y)  RMSE(x, y),
                                           y = c(cluster$rtMin, cluster$rtMax)) <= 2*rt_tol,]
    }
    
    if (nrow(cluster_peak) > 1) # there is an adj cluster that was not merged
    {
      # assign the spectra with the bigger base peak intensity as the principal and the others as fragmented clusters
      j <- which.max(cluster_peak$basePeakInt)
      #cat("msclusterID ", cluster_peak$msclusterID[[j]], " have fragmented clusters.\n")
      clean_count_table$fragmented_clusters[cluster_peak$idx[[j]]] <- nrow(cluster_peak) - 1
      cluster_peak <- cluster_peak[-j,]
      clean_count_table$fragmented_clusters[cluster_peak$idx] <- -1
      
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
  clean_count_table$idx <- NULL
  return(clean_count_table$fragmented_clusters)
  #write_csv(clean_count_table, sub(".csv","_fragmentedClusters.csv", count_table_path))
}

