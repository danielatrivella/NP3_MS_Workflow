##
# return the top n hits above the cutoff, if no cadidate is selected 
# return the top k below the cutoff
# sort the table by the maximum biocorrelation and basePeakIntensity to better 
# deal with ties in the top n retrieval
# do not select blansk or bed m/zs
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

cat("Loading packages readr, dplyr...\n")
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
source(file.path(script_path(), "read_metadata_table.R"))

corr_cutoff <- 0.6
top_n_hits_gecutoff <- 5
top_k_hits_lecutoff <- 3
rm_flags <- TRUE

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("Two arguments must be supplied to retrieve the top hits:\n", 
       " 1 - Path to a to a NP3 count table;\n", 
       " 2 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks;\n",
       " 3 - The biocorrelation cutoff to select the top n hits;\n", 
       " 4 - The top n hits to be retrieve with the biocorrelation above the cutoff;\n", 
       #" 5 - The top k hist to be retrieve when there were less then top n retrieved, where the top k are below the biocorr cutoff;\n", 
       " 5 - TRUE to also remove bflag and bedflag, FALSE otherwise;\n", 
       call.=FALSE)
} else {
  path_count <- file.path(args[[1]])
  if (!file.exists(path_count))
  {
    stop("The count file '", path_count,
         "' do not exists. Provide a valid path to where csv file with the count ",
         "of spectra is located.")
  }
  
  path_metadata <- file.path(args[[2]])
  if (!file.exists(path_metadata))
  {
    stop("The CSV batch metadata file '", path_metadata, 
         "' do not exists. Provide a valid path to where the metadata is located.")
  }
  
  corr_cutoff <- as.numeric(args[[3]])
  top_n_hits_gecutoff <- as.numeric(args[[4]])
  #top_k_hits_lecutoff <- as.numeric(args[[5]])
  rm_flags <- as.logical(args[[5]])
}

# read the quantifications, skipping the rows with the bioactivities
ms_count <- suppressMessages(read_csv(path_count, guess_max = 10000, 
                                      skip = which(is.na(read_csv_chunked(path_count, 
                                                                          ListCallback$new(function(x, pos) 
                                                                            x[[1]]), 
                                                                          col_names =F)[[1]])), col_names = TRUE))
if (!('msclusterID' %in% names(ms_count))) {
  stop("ERROR. Missing the msclusterID column or skipt more rows than necessary (",
       which(is.na(read_csv_chunked(path_count, 
                                    ListCallback$new(function(x, pos) 
                                      x[[1]]), 
                                    col_names =F)[[1]])),")")
}
metadata <- readMetadataTable(path_metadata)

corr_cols <- names(ms_count)[startsWith(names(ms_count), 'COR_')]
# compute the max correlation
ms_count$max_corr <- do.call(pmax, c(ms_count[,corr_cols], na.rm=T))
ms_count$BCC_top_hits <- NA
ms_count$valid <- TRUE
# make blanks and bed mzs invalid for selection
if ('BLANKS_TOTAL' %in% names(ms_count)) {
  ms_count$valid[ms_count$BLANKS_TOTAL > 0] <- FALSE
  if (rm_flags) {
    ms_count$valid[ms_count$BFLAG] <- FALSE
  }
}
if ('BEDS_TOTAL' %in% names(ms_count)) {
  ms_count$valid[ms_count$BEDS_TOTAL > 0] <- FALSE
  if (rm_flags) {
    ms_count$valid[ms_count$BEDFLAG] <- FALSE
  }
}

for (corr_col in corr_cols) {
  # order ms_count table using the max correlation and the basePeakInt
  top_n_corr <- ms_count %>%
    filter(valid == TRUE,  !is.na(!!as.symbol(corr_col)), !!as.symbol(corr_col) >= corr_cutoff) %>%
    arrange(desc(!!as.symbol(corr_col)), desc(max_corr), desc(basePeakInt)) %>% 
    #slice(1:max(top_n_hits_gecutoff, top_k_hits_lecutoff)) %>% 
    select('msclusterID', !!as.symbol(corr_col))
  cat("\n** Top selected for ", corr_col, " **\n")
  print(top_n_corr)
  # select the top n above the cutoff
  consensus_spec_select <- top_n_corr[, 'msclusterID'][[1]]
  if (length(consensus_spec_select) >= top_n_hits_gecutoff) {
    consensus_spec_select <- consensus_spec_select[1:top_n_hits_gecutoff]
  } 
  #else {
  #  # not enought candidates above the cutoff, select more candidates to reach top_k
  #  if (top_k_hits_lecutoff > length(consensus_spec_select)) {
  #    consensus_spec_select <- top_n_corr[1:min(top_k_hits_lecutoff, nrow(top_n_corr)), 'msclusterID'][[1]]
  #  }
  #}
  
  # annotate the selected candidates
  #print(consensus_spec_select)
  select_index <- match(consensus_spec_select, ms_count$msclusterID)
  #print(select_index)
  ms_count$BCC_top_hits[select_index] <- sapply(seq_along(select_index), 
                                                function(i, sl_idxs, corr_name) {
    ifelse(is.na(ms_count$BCC_top_hits[sl_idxs[i]]), 
           paste0(corr_name,"_top_",i),
           paste0(ms_count$BCC_top_hits[sl_idxs[i]],";",corr_name,"_top_",i))
  }, sl_idxs = select_index, corr_name = corr_col)
}

ms_count$max_corr <- NULL
ms_count$valid <- NULL
#ms_count$BCC_top_hits[!is.na(ms_count$BCC_top_hits)]

write_csv(ms_count, sub(".csv", "_topHits_BCC.csv", path_count))
