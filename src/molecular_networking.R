##
# Step 10 - creates the SSMN from pairwise similarity tables in sparse format with 4 columns
#           Use awk to filter the complete pairwise similarity table
##

sim_min <- as.numeric(0.6)

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Two arguments must be supplied to create the molecular networking of similarities (MN):\n", 
       " 1 - Path to the output data folder, inside the 'outs' directory of the clustering results folder. The resulting MN edge table will be saved inside the molecular_network dir;\n", 
       " 2 - Minimum cosine score that must occur between a pair of consensus MS/MS spectra in order for an edge to be formed in the molecular network. Lower values will increase the size of the clusters by inducing the clustering of less related MS/MS spectra, higher values will limit to the opposite;\n", 
       call.=FALSE)
} else {
  output_path <- file.path(args[[1]])
  if (!dir.exists(output_path))
  {
    stop("The job output folder '", output_path, 
         "' do not exists. Provide a valid path to where the the job final result is located.")
  }
  output_name <- basename(output_path)
  
  path_sim_table <- file.path(output_path, "molecular_networking", "similarity_tables", 
                              paste0("similarity_table_", output_name, "_clean.csv"))
  if (!file.exists(path_sim_table)) 
  {
    stop("The final clean pairwise similarity table does not exists: ", path_sim_table,
         ". Look for errors in the clean Step 5 and retry.")
  }
  path_clean_table <- file.path(output_path, "count_tables", "clean", 
                              paste0(output_name, "_peak_area_clean_ann.csv"))
  if (!file.exists(path_clean_table)) 
  {
    # try not annotated counts
    path_clean_table <- file.path(output_path, "count_tables", "clean", 
                                  paste0(output_name, "_peak_area_clean.csv"))
    if (!file.exists(path_clean_table)) 
    {
      stop("The final clean count table does not exists: ", path_clean_table,
           ". Look for errors in the clean Step 5 and retry.")
    }
  }
  output_path <- file.path(output_path, "molecular_networking") 
  
  sim_min <- round(as.numeric(args[[2]]), 2)
  if (is.na(sim_min)) {
    stop("The provided similarity cutoff to connect a pair of spectra have ", 
         "an invalid value. Please correctly set this parameter value.")
  }
}

cat("Loading package readr...\n")
suppressPackageStartupMessages(library(readr))

build_mol_net_sim <- function(output_name, path_sim_table, path_clean_table,
                              output_path, sim_min)
{
  if (!file.exists(path_sim_table))
  {
    stop("The pairwise similarity table file '", path_sim_table,
         "' does not exists. Provide a valid path to where it is located.")
  }
  if (!file.exists(path_clean_table))
  {
    stop("The clean count table file '", path_sim_table,
         "' does not exists. Provide a valid path to where it is located.")
  }
  if (!dir.exists(output_path))
  {
    stop("The output folder '", output_path, 
         "' does not exists. Provide a valid path to where the MN should be saved.")
  }
  if (sim_min < 0 || sim_min > 1.0)
  {
    stop("The minimum similarity threshold must be a non negative numeric value less or equal than 1",
         " (x in [0,1.0]). Wrong sim_min value: ", sim_min)
  }
  
  # read the clean table to retrieve the complete list of msclusterIDs
  clean_counts_table <- suppressMessages(read_csv(path_clean_table, 
                                                  guess_max = 5000,
                                                  col_types = cols(.default="?", 
                                                                   msclusterID="i")))
  n_scans <- nrow(clean_counts_table)

  cat("\n** Creating the complete Spectra Similarity Molecular Networking (SSMN) with", n_scans, "nodes **\n")
  # create the complete SSMN by filtering the sim > cutoff from the clean sim table
  # use awk to make this filtering direct in the file without extra RAM usage
  # Start a timer to track performance
  t0 <- Sys.time()
  ssmn_file_path <- file.path(output_path, paste0(output_name,"_ssmn_w_",sub("\\.", "", sim_min),".selfloop"))
  # Run awk to filter the provided similarity table
  # stdout = sim_file streams directly to disk
  status <- system2(
    command = "awk",
    args = c(
      "-F,",
      "-v", paste0("sim_cutoff=", sim_min),
      "'NR == 1 || $3 >= sim_cutoff'",
      path_sim_table
    ),
    stdout = ssmn_file_path # Captures errors if something goes wrong
  )
  
  # check if status was 0 for success
  if (status != 0 || !file.exists(ssmn_file_path))
  {
    cat("    - ERROR. awk Status:",status,". Finished in:", round(tf - t0, 2), "seconds\n")
    stop("The final clean pairwise similarity table file '", path_sim_table,
         "' could not be filtered with awk to create the complete SSMN. ",
         "Check for warnings and the read/write permissions of the current user.")
  }
  
  # add self loops, check missing nodes in the complete ssmn
  # retrieve the list of nodes that appear in the SSMN
  msclusterIDs_in_ssmn <- sort(sapply(system2(
    command = "awk",
    args = c(
      "-F,",
      "'NR > 1 {if (!seen[$1]++) print $1; if (!seen[$2]++) print $2}'",
      ssmn_file_path
    ),
    stdout = TRUE), 
    as.numeric))
  if (msclusterIDs_in_ssmn[[1]] == 0) {
    stop("Error retrieving the msclusterIDs present in the complete SSMN, ",
         "something went wrong in the awk command:", msclusterIDs_in_ssmn)
  }
  # create isolated nodes check
  clean_counts_table$ssmn_complete_isolated <- !(clean_counts_table$msclusterID %in% msclusterIDs_in_ssmn)
  # retrieve the isolated nodes ids
  msclusterIDs_isolated <- clean_counts_table$msclusterID[clean_counts_table$ssmn_complete_isolated]
  
  # write the single nodes as self loops  for the isolated ones
  edges_isolated <- data.frame(msclusterID_source = msclusterIDs_isolated,
                         msclusterID_target = msclusterIDs_isolated,
                         similarity_value = rep(1.00, length(msclusterIDs_isolated)),
                         num_matched_peaks = rep(-1, length(msclusterIDs_isolated)),
                         stringsAsFactors = FALSE)
  
  readr::write_csv(edges_isolated, path = ssmn_file_path,
                   append = TRUE, col_names = FALSE)
  write_csv(clean_counts_table, path_clean_table)
  tf <- Sys.time()
  cat("  * Done in", round(tf-t0, 2), units(tf-t0), "*\n")
  
  cat("\n  * Retrieve the number of peaks and ionization annotations of the complete SSMN nodes *\n")
  
  # retrieve number of peaks
  msclusterIDs_num_peaks <- sapply(clean_counts_table$peaksList,
                                   stringi::stri_count_fixed, pattern=";", 
                                   USE.NAMES = FALSE) + 1
  names(msclusterIDs_num_peaks) <- clean_counts_table$msclusterID
  rm(clean_counts_table) # free memory
  
  # read the complete SSMN and add more information to the edges
  edges_ssmn <- suppressMessages(read_csv(ssmn_file_path, 
                                          guess_max = 5000,
                                          col_types = cols(.default="?", 
                                                           msclusterID_source="i",
                                                           msclusterID_target="i")))
  
  # add the number of peaks information of the source and target nodes to their edges, 
  # to be used when applying the min matched peaks parameter
  edges_ssmn$num_peaks_source <- msclusterIDs_num_peaks[match(edges_ssmn$msclusterID_source, names(msclusterIDs_num_peaks))]
  edges_ssmn$num_peaks_target <- msclusterIDs_num_peaks[match(edges_ssmn$msclusterID_target, names(msclusterIDs_num_peaks))]
  rm(msclusterIDs_num_peaks) # free memory
  
  # add the ionization variants annotations in the SSMN if the IVAMN is present
  ivamn_path <- file.path(output_path, paste0(output_name, "_ivamn.selfloop"))
  if (file.exists(ivamn_path))
  {
    edges_ivamn <- suppressMessages(readr::read_csv(ivamn_path))
    edges_ivamn <- edges_ivamn[!is.na(edges_ivamn$annotation),]
    
    edges_ssmn$annotation <- apply(edges_ssmn, 1, function(x, ann_edges)
    {
      match_undirect_edge <- (ann_edges$msclusterID_source == x[[1]] & 
                                  ann_edges$msclusterID_target == x[[2]]) | 
                                 (ann_edges$msclusterID_source == x[[2]] &
                                  ann_edges$msclusterID_target == x[[1]])
      if (any(match_undirect_edge)) {
        paste(unlist(ann_edges[match_undirect_edge, "annotation"]), collapse = ";") 
      } else {
        ""
      }
    }, edges_ivamn)
  }
  # write final complete ssmn
  readr::write_csv(edges_ssmn, path = ssmn_file_path)
 
  tend <- Sys.time()
  cat("    * Done in", round(tend-tf, 2), units(tend-tf), "*\n")
  cat("** Done creating the molecular network of similarities in", 
          round(tend-t0, 2), units(tend-t0), "**\n")
} 

build_mol_net_sim(output_name, path_sim_table, path_clean_table,
                  output_path, sim_min)
