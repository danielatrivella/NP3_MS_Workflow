##
# test the clean counts and the molecular networks consistency
##
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
Rcpp::sourceCpp(file.path(script_path(),
                          '../src/triangular_matrix_R.cpp'))
Rcpp::sourceCpp(file.path(script_path(),
                          '../src/read_mgf_peak_list_R.cpp'))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))

sim_tol <- 0.55
mn_tol <- 0.6
rt_tol <- 2
mz_tol <- 0.025
top_k <- 15

RMSE <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  sqrt(mean((x - y)^2))
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
  stop("Six arguments must be supplied to test the clean and annotated counts and the molecular networking consistency:\n", 
       " 1 - Path to the output folder where the molecular_networking directory is located;\n",
       " 2 - sim_tol;\n",
       " 3 - mn_tol;\n",
       " 4 - rt_tol;\n",
       " 5 - mz_tol;\n",
       " 6 - top_k;\n",
       " 7 - max_component_size;\n",
       " 8 - min_matched_peaks.",
       call.=FALSE)
} else {
  output_path <- args[[1]]
  sim_tol <- as.numeric(args[[2]])
  mn_tol <- as.numeric(args[[3]])
  rt_tol <- as.numeric(args[[4]])
  mz_tol <- as.numeric(args[[5]])
  top_k <- as.numeric(args[[6]])
  max_component_size <- as.numeric(args[[7]])
  min_matched_peaks <- as.numeric(args[[8]])
}
# print(args)
output_name <- basename(output_path)

path_area_count <- file.path(output_path, "count_tables", "clean", paste0(output_name,"_peak_area_clean_ann.csv"))
if (!file.exists(path_area_count))
{
  stop("The count file '", path_area_count,
       "' do not exists. Provide a valid output path to where csv file with the count ",
       "of peak area is located.")
}
ms_area_count <- suppressMessages(read_csv(path_area_count, guess_max = 10000))

path_sim_table <- file.path(output_path, "molecular_networking", "similarity_tables", 
                            paste0("similarity_table_", output_name, "_aggMax.csv"))

if (!file.exists(path_sim_table))
{
  stop("The pairwise similarity table file '", path_sim_table,
       "' do not exists. Provide a valid path to where it is located.")
}

path_sim_selfloops <- file.path(output_path, "molecular_networking", 
                                paste0(output_name,"_ssmn_w_",
                                        sub("\\.", "", mn_tol),
                                       "_mmp_",min_matched_peaks,
                                       "_k_", top_k, "_x_",
                                        max_component_size, ".selfloop"))
path_ann_selfloops <- file.path(output_path, "molecular_networking", 
                                paste0(output_name,
                                       "_ivamn.selfloop"))

scans_pairsim <- suppressMessages(read_csv(path_sim_table))
scans_order <- c(-1, scans_pairsim[[1]])

sim_selfloops <- suppressMessages(read_csv(path_sim_selfloops))
ann_selfloops <- suppressMessages(read_csv(path_ann_selfloops))
n_inconsistency <- 0

# check if all msclusterIDs are present in both networks
if (!all(ms_area_count$msclusterID %in% sim_selfloops$msclusterID_source | 
         ms_area_count$msclusterID %in% sim_selfloops$msclusterID_target)) {
  cat("The following msclusterIDs are missing in the molecular network of similarities: ",
      paste(ms_area_count$msclusterID[!(ms_area_count$msclusterID %in% sim_selfloops$msclusterID_source | 
                                          ms_area_count$msclusterID %in% sim_selfloops$msclusterID_target)], 
            collapse = ","), '\n')
  n_inconsistency <- n_inconsistency + length(ms_area_count$msclusterID[
    !(ms_area_count$msclusterID %in% sim_selfloops$msclusterID_source |
        ms_area_count$msclusterID %in% sim_selfloops$msclusterID_target)])
}
if (!all(ms_area_count$msclusterID %in% ann_selfloops$msclusterID_source | 
         ms_area_count$msclusterID %in% ann_selfloops$msclusterID_target)) {
  cat("The following msclusterIDs are missing in the molecular network of annotations: ",
      paste(ms_area_count$msclusterID[!(ms_area_count$msclusterID %in% ann_selfloops$msclusterID_source | 
                                          ms_area_count$msclusterID %in% ann_selfloops$msclusterID_target)], 
            collapse = ","), '\n')
  n_inconsistency <- n_inconsistency + length(ms_area_count$msclusterID[
    !(ms_area_count$msclusterID %in% ann_selfloops$msclusterID_source |
        ms_area_count$msclusterID %in% ann_selfloops$msclusterID_target)])
}
# check if all edges of the filtered SSMN have more peaks in common than the 
# min_matched_peaks, excluding selfloops (num_matched_peaks = -1) and edges
# between spectra with few peaks < min_matched_peaks, which should have at least 
# more than 2 matched peaks
if (any(sim_selfloops$num_matched_peaks < min_matched_peaks & 
        sim_selfloops$num_matched_peaks != -1 & 
        ((sim_selfloops$num_peaks_source > min_matched_peaks & 
          sim_selfloops$num_peaks_target > min_matched_peaks) | 
         (sim_selfloops$num_matched_peaks < 2 & 
          (sim_selfloops$num_peaks_source < min_matched_peaks | 
           sim_selfloops$num_peaks_target < min_matched_peaks))))) {
  cat("There are", sum(sim_selfloops$num_matched_peaks < min_matched_peaks & 
                         sim_selfloops$num_matched_peaks != -1), 
      "connections in the filtered SSMN with less common peaks than the minimum number",
      "of matched peaks =", min_matched_peaks, 
      ". Error in the SSMN filtering algorithm.")
  n_inconsistency <- n_inconsistency + sum(sim_selfloops$num_matched_peaks < min_matched_peaks & 
                                             sim_selfloops$num_matched_peaks != -1)
}
# check if all edges of the filtered SSMN have a cosine above the cutoff
if (any(sim_selfloops$cosine < mn_tol)) {
  cat("There are", sum(sim_selfloops$cosine < mn_tol), 
      "connections in the filtered SSMN with a cosine value smaller than the",
      "similarity cutoff =", mn_tol,".")
  n_inconsistency <- n_inconsistency + sum(sim_selfloops$cosine < mn_tol1)
}

# remove NA annotations of selfloops
ann_selfloops <- ann_selfloops[!is.na(ann_selfloops$annotation),] 

# check if similarity cutoff the fragments, isotopes and neutral loss annotations is being respected
if (any(grepl("fragment", ann_selfloops$annotation) & 
        ann_selfloops$cosine < 0.2))
{
  cat("There are", sum((grepl("fragment", ann_selfloops$annotation) & 
                          ann_selfloops$cosine < 0.2)),
      "wrong 'fragment' annotations with the cosine similarity bellow the cut-off of 0.2.\n")
  n_inconsistency <- n_inconsistency + sum((grepl("fragment", ann_selfloops$annotation) & 
                                              ann_selfloops$cosine < 0.2))
}
if (any(grepl("\\[M\\+1\\]\\+", ann_selfloops$annotation) & 
        ann_selfloops$cosine < 0.75))
{
  cat("There are", sum((grepl("\\[M\\+1\\]\\+", ann_selfloops$annotation) & 
                          ann_selfloops$cosine < 0.75)),
      "wrong '[M+1]+' isotope annotations with the cosine similarity bellow the cut-off of 0.75.\n")
  n_inconsistency <- n_inconsistency + sum((grepl("\\[M\\+1\\]\\+", ann_selfloops$annotation) & 
                                              ann_selfloops$cosine < 0.75))
}
if (any(grepl("\\[M\\+H-(NH3|H2O|NH3-H2O)]\\+", ann_selfloops$annotation) & 
        ann_selfloops$cosine < 0.2))
{
  cat("There are", sum((grepl("\\[M\\+H-(NH3|H2O|NH3-H2O)]\\+", 
                              ann_selfloops$annotation) & 
                          ann_selfloops$cosine < 0.2)),
      "wrong '[M+H-NH3]+' or '[M+H-H2O]+' or '[M+H-NH3-H2O]+' neutral loss",
      "annotations with the cosine similarity bellow the cut-off of 0.2.\n")
  n_inconsistency <- n_inconsistency + sum((grepl("\\[M\\+H-(NH3|H2O|NH3-H2O)]\\+", 
                                                  ann_selfloops$annotation) & 
                                              ann_selfloops$cosine < 0.2))
}

# check if readMgfPeaksList and the readMgfHeader function are consistent 
# and match the table attributes
path_mgf_clean <- file.path(output_path, "mgf", 
                            paste0( output_name, "_clean.mgf"))
mgf_clean <- readMgfPeaksList(path_mgf_clean, mz_tol, -1, 1, 0)
mgf_clean_header <- readMgfHeader(path_mgf_clean)
if (!all(mgf_clean$SCANS == mgf_clean_header$scans)) {
  cat('The readMgfPeaksList and the readMgfHeader functions have inconsistent',
      'results, some scans do not match. A total of', 
      sum(!(mgf_clean$SCANS == mgf_clean_header$scans)), 'scans were different.\n')
  n_inconsistency <- n_inconsistency + sum(!(mgf_clean$SCANS == mgf_clean_header$scans))
}
if (!all(mgf_clean$PREC_MZ == mgf_clean_header$mz)) {
  cat('The readMgfPeaksList and the readMgfHeader functions have inconsistent',
      'results, some precursor m/zs do not match. A total of', 
      sum(!(mgf_clean$PREC_MZ == mgf_clean_header$mz)), 
      'precursor m/zs were different.\n')
  n_inconsistency <- n_inconsistency + sum(!(mgf_clean$PREC_MZ == mgf_clean_header$mz))
}
# check if the mgfs headers match the table attributes' values
if (!(length(mgf_clean$SCANS) == nrow(ms_area_count) && 
      length(mgf_clean_header$scans) == nrow(ms_area_count))) {
  cat('The number of spectra read by the readMgfPeaksList and the readMgfHeader',
      'function do not match the number of rows in the clean count table.',
      '\nNumber of spectra read by the readMgf* functions:', length(mgf_clean$SCANS),
      '\nNumber of rows (spectra) in the clean count table:', nrow(ms_area_count),
      '\n')
  n_inconsistency <- n_inconsistency + 1
}
n_inconsistency_mgf <- 0
# dealing with floating point errors in the last decimal place of each field
for (i in seq_len(length(mgf_clean$SCANS))) {
  j <- which(ms_area_count$msclusterID == mgf_clean$SCANS[[i]])
  
  if (abs(mgf_clean$PREC_MZ[[i]] - ms_area_count$mzConsensus[[j]]) > 0.0005)
    n_inconsistency_mgf <- n_inconsistency_mgf + 1
  if (abs(mgf_clean_header$mz[[i]] - ms_area_count$mzConsensus[[j]]) > 0.0005)
    n_inconsistency_mgf <- n_inconsistency_mgf + 1
  if (abs(mgf_clean_header$rt[[i]] - ms_area_count$rtMean[[j]]) > 0.0005)
    n_inconsistency_mgf <- n_inconsistency_mgf + 1
  if (abs(mgf_clean_header$rtmin[[i]] - ms_area_count$rtMin[[j]]) > 0.0005)
    n_inconsistency_mgf <- n_inconsistency_mgf + 1
  if (abs(mgf_clean_header$rtmax[[i]] - ms_area_count$rtMax[[j]]) > 0.0005)
    n_inconsistency_mgf <- n_inconsistency_mgf + 1
  if (abs(mgf_clean_header$int[[i]] - ms_area_count$sumInts[[j]]) > 1)
    n_inconsistency_mgf <- n_inconsistency_mgf + 1
  # if (n_inconsistency_mgf >0)
  #   break()
}
if (n_inconsistency_mgf > 0)
{
  cat('A total of', n_inconsistency_mgf, 'inconsistencies were found between',
      'the mgf header values and the clean count table attributes.',
      'Something went wrong in the clean step or in the mgf parsing.\n')
  n_inconsistency <- n_inconsistency + n_inconsistency_mgf
}
remove(n_inconsistency_mgf, mgf_clean, mgf_clean_header)

if ('BLANKS_TOTAL' %in% names(ms_area_count)) {
  blanks_flag <- TRUE
} else {
  blanks_flag <- FALSE
}
# check clean and annoations - count
for (i in seq_len(nrow(ms_area_count))) {
  # print(i)
  cluster <- ms_area_count[i,]
  
  # check the similarity network
  # get the clusters ids that have a similarity with the current cluster including itself
  adj_clusters <- scans_order[which_ge(unlist(scans_pairsim[i,-1]), sim_tol, 1)] # add one to scape first column
  
  cluster_peak <- ms_area_count[abs(ms_area_count$mzConsensus - cluster$mzConsensus) <= mz_tol &
                                  ((ms_area_count$rtMean >= cluster$rtMin - rt_tol & 
                                      ms_area_count$rtMean <= cluster$rtMax + rt_tol) |
                                     (cluster$rtMean >= ms_area_count$rtMin - rt_tol &
                                        cluster$rtMean <= ms_area_count$rtMax + rt_tol)),] 
  # if not bflag check peak center and boundaries deviation, remove peaks not aligned
  # the spectra peak center deviation is <= 4 * rt_tol or peak boundaries deviation <= 2*rt_tol
  if (!blanks_flag || !cluster$BFLAG) {
    cluster_peak <- cluster_peak[abs(cluster$rtMean-cluster_peak$rtMean) <= 4*rt_tol |
                                   apply(cluster_peak[,c("rtMin", "rtMax")], 1,
                                         function(x,y)  RMSE(x, y),
                                         y = c(cluster$rtMin, cluster$rtMax)) <= 2*rt_tol,]
  }
  # if there is a not similar cluster in the current cluster peak check if they share any MS1 peak id 
  # and add them to the adj list to be joined
  # do not check clusters that have too much spectra
  non_adj_peak <- !(cluster_peak$msclusterID[cluster_peak$numSpectra < 5000] %in% adj_clusters)
  if (any(non_adj_peak) && cluster$numSpectra < 5000) {
    # do not consider fake peaks ids
    peakIds <- strsplit(cluster$peakIds,";")[[1]]
    # peaksIds <- paste0(peakIds,collapse = "|")
    peakIds <- peakIds[!startsWith(peakIds, "fake_")]
    if (length(peakIds) > 0 && length(peakIds) <= 500) {
      peaksIds <- paste0(peakIds, collapse = "|")
      non_adj_peak <- cluster_peak$msclusterID[cluster_peak$numSpectra < 5000][non_adj_peak][
        grepl(pattern = peaksIds, 
              cluster_peak$peakIds[cluster_peak$numSpectra < 5000][non_adj_peak])&
        (sapply(cluster_peak$peakIds[cluster_peak$numSpectra < 5000][non_adj_peak], 
                str_count, pattern=";", USE.NAMES = FALSE) <= 500)]
      if (length(non_adj_peak) > 0)
        adj_clusters <- c(adj_clusters, non_adj_peak)
    }
  }
  # filter only adj peaks
  cluster_peak <- cluster_peak[cluster_peak$msclusterID %in% adj_clusters,]
  if (nrow(cluster_peak) > 1) # there is an adj cluster that was not merged
  {
    n_inconsistency <- n_inconsistency + 1
    cat("msclusterID ", cluster[[1]], " was not joined correctly in the clean step.\n")
  }
  if (scans_pairsim[i,1] != cluster[[1]]) # there is an adj cluster that was not merged
  {
    n_inconsistency <- n_inconsistency + 1
    cat("msclusterID ", cluster[[1]], " wrong order with the sim table\n")
  }
  
  # TODO this should be done with the clean similarity table based on the clean mgf, which is now used in the mn creation
  # check adj clusters in the mn
  # adj_clusters <- length(c(scans_order[which_ge(unlist(scans_pairsim[i,-1]), mn_tol, 1)][-1], # add one to scape first column and rm selfloop
  #                          scans_order[rev(which_ge(unlist(scans_pairsim[,i+1]), mn_tol, 0))][-1])) # should be removing the last not the first, in the column the selfloop is the last not zero sim
  # if (adj_clusters == 0)
  #   adj_clusters <- 1 # self loop
  if (sum(sim_selfloops$msclusterID_source == cluster[[1]] | 
              sim_selfloops$msclusterID_target == cluster[[1]]) >
      # min(adj_clusters, top_k) ) {
      top_k ) {
    n_inconsistency <- n_inconsistency + 1
    cat("msclusterID", cluster[[1]], "wrong sim edges number. Number of edges ", 
            sum(sim_selfloops$msclusterID_source == cluster[[1]] | sim_selfloops$msclusterID_target == cluster[[1]]),
            # " Expected number of edges: ", min(adj_clusters, top_k),"\n" )
            " Expected number of edges: ", top_k,"\n" )
  }
  
  # CHECK if the count ann size matches the number of ann 
  cluster_anns <- paste(na.exclude(unlist(cluster[c("adducts", "isotopes", 
                                                    "dimers", "multiCharges", 
                                                    "fragments")])), collapse = ";")
  cluster_anns <- strsplit(cluster_anns, ";")[[1]]
  cluster_anns <- (as.numeric(regmatches(cluster_anns, 
                                               regexpr("(?<=\\[)[0-9]+(?=\\])", 
                                                       cluster_anns, perl = T))))                        
  
  ann_edges <- strsplit(paste(ann_selfloops[ann_selfloops$msclusterID_target == cluster$msclusterID, 
                                            "annotation"][[1]], collapse = ";"), ";")[[1]]
  if (length(cluster_anns) != length(ann_edges))
  {
    n_inconsistency <- n_inconsistency + 1
    cat("msclusterID", cluster[[1]], "annotations were not written to annotation network correctly.\n")
  }
}

if (n_inconsistency == 0) {
  cat("\nDone! :)\n")
} else {
  cat("\nERROR! A total of", n_inconsistency, "inconsistencies were found. :(\n")
}
