##
# test the consistency of the molecular networks
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
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))

mn_tol <- 0.6
top_k <- 15

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Three arguments must be supplied to test the molecular networking consistency:\n", 
       " 1 - Path to the output folder where the molecular_networking directory is located;\n",
       " 2 - mn_tol;\n",
       " 3 - top_k;\n",
       " 4 - max_component_size;\n",
       " 5 - min_matched_peaks.",
       call.=FALSE)
} else {
  output_path <- args[[1]]
  # if this is the testing of the mn of annotations, skip the checking of the similarity mn
  if (args[[2]] == "undefined" && args[[3]] == "undefined" && args[[4]] == "undefined") {
    mn_tol <- NA
  } else {
    mn_tol <- as.numeric(args[[2]])
    top_k <- as.numeric(args[[3]])
    max_component_size <- as.numeric(args[[4]])
    min_matched_peaks <- as.numeric(args[[5]])
  }
}

output_name <- basename(output_path)

path_area_count <- file.path(output_path, "count_tables", "clean", paste0(output_name,"_peak_area_clean_ann.csv"))
if (!file.exists(path_area_count))
{
  stop("The count file '", path_area_count,
       "' do not exists. Provide a valid output path to where csv file with the count ",
       "of peak area is located.")
}
ms_area_count <- suppressMessages(read_csv(path_area_count, guess_max = 10000))

# only test the mn of similarities if a similarity tolerance was passed
if (!is.na(mn_tol)) {
  path_sim_table <- file.path(output_path, "molecular_networking", "similarity_tables", 
                              paste0("similarity_table_", output_name, "_clean.csv"))
  if (!file.exists(path_sim_table))
  {
    stop("The pairwise similarity table file '", path_sim_table,
         "' do not exists. Provide a valid path to where it is located.")
  }
  scans_pairsim <- suppressMessages(read_csv(path_sim_table))
  scans_order <- c(-1, scans_pairsim[[1]])
  
  path_sim_selfloops <- file.path(output_path, "molecular_networking", 
                                  paste0(output_name,"_ssmn_w_",
                                         sub("\\.", "", mn_tol),
                                         "_mmp_",min_matched_peaks,
                                         "_k_", top_k, "_x_",
                                         max_component_size, ".selfloop"))
  sim_selfloops <- suppressMessages(read_csv(path_sim_selfloops))
}

path_ann_selfloops <- file.path(output_path, "molecular_networking", paste0(output_name,"_ivamn.selfloop"))
ann_selfloops <- suppressMessages(read_csv(path_ann_selfloops))
n_inconsistency <- 0

# check if all msclusterIDs are present in both networks
if (!is.na(mn_tol) && !all(ms_area_count$msclusterID %in% sim_selfloops$msclusterID_source | 
        ms_area_count$msclusterID %in% sim_selfloops$msclusterID_target)) {
  cat("The following msclusterIDs are missing in the molecular network of similarities: ",
      paste(ms_area_count$msclusterID[!(ms_area_count$msclusterID %in% sim_selfloops$msclusterID_source | 
                                         ms_area_count$msclusterID %in% sim_selfloops$msclusterID_target)], 
            collapse = ","))
  n_inconsistency <- n_inconsistency + length(ms_area_count$msclusterID[
    !(ms_area_count$msclusterID %in% sim_selfloops$msclusterID_source |
      ms_area_count$msclusterID %in% sim_selfloops$msclusterID_target)])
}
if (!all(ms_area_count$msclusterID %in% ann_selfloops$msclusterID_source | 
         ms_area_count$msclusterID %in% ann_selfloops$msclusterID_target)) {
  cat("The following msclusterIDs are missing in the molecular network of annotations: ",
      paste(ms_area_count$msclusterID[!(ms_area_count$msclusterID %in% ann_selfloops$msclusterID_source | 
                                          ms_area_count$msclusterID %in% ann_selfloops$msclusterID_target)], 
            collapse = ","))
  n_inconsistency <- n_inconsistency + length(ms_area_count$msclusterID[
    !(ms_area_count$msclusterID %in% ann_selfloops$msclusterID_source |
        ms_area_count$msclusterID %in% ann_selfloops$msclusterID_target)])
}

# check if all edges of the filtered SSMN have more peaks in common than the 
# min_matched_peaks, excluding selfloops (num_matched_peaks = -1) and edges
# between spectra with few peaks < min_matched_peaks, which should have at least 
# more than 2 matched peaks
if (!is.na(mn_tol) && any(sim_selfloops$num_matched_peaks < min_matched_peaks & 
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
if (!is.na(mn_tol) && any(sim_selfloops$cosine < mn_tol)) {
  cat("There are", sum(sim_selfloops$cosine < mn_tol), 
      "connections in the filtered SSMN with a cosine value smaller than the",
      "similarity cutoff =", mn_tol,".")
  n_inconsistency <- n_inconsistency + sum(sim_selfloops$cosine < mn_tol)
}

# remove NA annotations of selfloops
ann_selfloops <- ann_selfloops[!is.na(ann_selfloops$annotation),] 

# check if the similarity cutoff for fragments, isotopes and neutral loss 
# annotations is being respected
if (any(grepl("fragment", ann_selfloops$annotation) & 
        ann_selfloops$cosine < 0.2))
{
  cat("There are", sum((grepl("fragment", ann_selfloops$annotation) & 
                        ann_selfloops$cosine < 0.2)),
      "wrong 'fragment' annotations with the cosine similarity bellow the cut-off of 0.2.")
  n_inconsistency <- n_inconsistency + sum((grepl("fragment", ann_selfloops$annotation) & 
                                              ann_selfloops$cosine < 0.2))
}
if (any(grepl("\\[M\\+1\\]\\+", ann_selfloops$annotation) & 
        ann_selfloops$cosine < 0.75))
{
  cat("There are", sum((grepl("\\[M\\+1\\]\\+", ann_selfloops$annotation) & 
                        ann_selfloops$cosine < 0.75)),
      "wrong '[M+1]+' isotope annotations with the cosine similarity bellow the cut-off of 0.75.")
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
      "annotations with the cosine similarity bellow the cut-off of 0.2.")
  n_inconsistency <- n_inconsistency + sum((grepl("\\[M\\+H-(NH3|H2O|NH3-H2O)]\\+", 
                                                  ann_selfloops$annotation) & 
                                              ann_selfloops$cosine < 0.2))
}

# check clean and annotations - count
for (i in seq_len(nrow(ms_area_count))) {
  # print(i)
  cluster <- ms_area_count[i,]
  
  # check adj clusters in the mn
  if (!is.na(mn_tol)) {
    adj_clusters <- length(c(scans_order[which_ge(unlist(scans_pairsim[i,-1]), mn_tol, 1)][-1], # add one to scape first column and rm selfloop
                             scans_order[rev(which_ge(unlist(scans_pairsim[,i+1]), mn_tol, 0))][-1])) # should be removing the last not the first, in the column the selfloop is the last not zero sim
    if (adj_clusters == 0)
      adj_clusters <- 1 # self loop
    if (sum(sim_selfloops$msclusterID_source == cluster[[1]] | 
            sim_selfloops$msclusterID_target == cluster[[1]]) >
        min(adj_clusters, top_k) ) {
      n_inconsistency <- n_inconsistency + 1
      cat("msclusterID", cluster[[1]], "wrong sim edges number. Number of edges ", 
              sum(sim_selfloops$msclusterID_source == cluster[[1]] | sim_selfloops$msclusterID_target == cluster[[1]]),
              " Expected number of edges: ", min(adj_clusters, top_k) ,"\n")
    }
  }
  
  # CHECK if the count ann size matches the number of ann 
  cluster_anns <- paste(na.exclude(unlist(cluster[c("adducts", "isotopes", 
                                                    "dimers", "multiCharges", 
                                                    "fragments")])), collapse = ";")
  cluster_anns <- strsplit(cluster_anns, ";")[[1]]

  ann_edges <- strsplit(paste(ann_selfloops[ann_selfloops$msclusterID_target == cluster$msclusterID, 
                             "annotation"][[1]], collapse = ";"), ";")[[1]]
  # cluster_anns <- cluster_anns[as.numeric(regmatches(cluster_anns, 
  #                                                    regexpr("(?<=\\[)[0-9]+(?=\\])", 
  #                                                            cluster_anns, perl = T))) > 
  #                                cluster[[1]]]
  if (length(cluster_anns) != length(ann_edges))
  {
    n_inconsistency <- n_inconsistency + 1
    cat("msclusterID", cluster[[1]], "annotations were not written to annotation network correctly.\n")
  }
}

if (n_inconsistency == 0) {
  cat("Done! :)\n")
} else {
  cat("ERROR! A total of", n_inconsistency, "inconsistencies were found. :(\n")
}
