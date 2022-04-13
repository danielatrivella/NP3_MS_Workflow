## ----load-libs, message = FALSE--------------------------------------------
cat("Loading functions...\n")
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
source(file.path(script_path(), "read_metadata_table.R"))
options(readr.show_progress = FALSE)

# default args values
mz_tol <- 0.025
rt_tol <- 2

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Two arguments must be supplied to clean and annotate the counts:\n", 
       " 1 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks;\n",
       " 2 - Path to the output file containing the csv file name;\n", 
       " 3 - Path to the quantification table of not fragmented MS1 peaks;\n",
       " 4 - mz tolerance;\n",
       " 5 - RT tolerance to enlarge the peaks before comparisions.\n",
       call.=FALSE)
} else {
  metadata_path <- file.path(args[[1]])
  if (!file.exists(metadata_path))
  {
    stop("The CSV batch metadata file '", metadata_path, 
         "' do not exists. Provide a valid path to where the metadata is located.")
  }
  
  output_path <- file.path(args[[2]])
  if (!dir.exists(dirname(output_path)))
  {
    stop("The output file '", output_path, 
         "' path do not exists. Provide a valid path to save the result.")
  }
  
  quantification_table_path <- file.path(args[[3]])
  if (!file.exists(quantification_table_path))
  {
    stop("The count file '", quantification_table_path,
         "' do not exists. Provide a valid path to where csv file with the MS1 not fragmented peaks count ",
         "of peak area is located.")
  }
  
  mz_tol <- as.numeric(args[[4]])
  rt_tol <- as.numeric(args[[5]])
}

# quantification_table_path <- "/home/crisfbazz/Documents/CNPEM/ms_entity_annotation/Data_Collections/control_helder/mzxml/processed_blankBaseline/mzs_no_MS2.csv"
# metadata_path <- "/home/crisfbazz/Documents/CNPEM/ms_entity_annotation/Data_Collections/control_helder/METADATA_control_helder.csv"
# output_path <- "/home/crisfbazz/Documents/CNPEM/ms_entity_annotation/Data_Collections/control_helder/control_test_anns_isoPattern/outs/control_test_anns_isoPattern/count_tables/control_test_anns_isoPattern_peak_area_MS1.csv"

# function to merge the count table by column
  merge_no_ms2_counts <- function(col_name, x)
{
  switch(col_name,
         peakIds = paste(x[[col_name]], collapse = ";"),
         numJoins =  sum(x$numJoins) + nrow(x) - 1,
         mz=,rtMean=
           round(weighted.mean(x[[col_name]], x$sumAreas), 4),
         rtMin=ifelse(any(x$rtMin == 0 & x$rtMax == 1000000), 0, round(weighted.mean(x[[col_name]], x$sumAreas), 4)),
         rtMax=ifelse(any(x$rtMin == 0 & x$rtMax == 1000000), 1000000, round(weighted.mean(x[[col_name]], x$sumAreas), 4)),
         sum(trunc(x[[col_name]]))) # sum unique area count columns and ms2_matches
}
 
# quantification_table path to the quantification table or a data.frame with it
# expected columns are: "peakIds", "mz", "rtMean", "rtMin", "rtMax", "ms2_matches" and the 
# metadata sample code names concatenated with _area suffix
# if in the output_path directory there exists the count of MS2 peak area, use the nrow of it to number the msclusterID
clean_quantification <- function(quantification_table, metadata, output_path, mz_tol, rt_tol) 
{
  if (is.character(metadata)) {
    if (!file.exists(metadata))
      stop("The CSV metadata file '", metadata, 
           "' do not exists. Provide a valid path to where the metadata is located.")
    metadata <- readMetadataTable(metadata)
  } else {
    metadata <- checkMetadataFormat(metadata)
  }
  metadata$SAMPLE_CODE <- paste0(metadata$SAMPLE_CODE, "_area")
  
  if (is.character(quantification_table)) {
    if (!file.exists(quantification_table))
      stop("The count file '", quantification_table,
           "' do not exists. Provide a valid path to where the csv ",
           "file with the count of peak area of not fragmented MS2 is located.")
    
    quantification_table <- suppressMessages(readr::read_csv(quantification_table, guess_max = 5000))
    if (!all(c("peakIds", "mz", "rtMean", "rtMin", "rtMax", 
               metadata$SAMPLE_CODE) %in% names(quantification_table))) {
      stop("Wrong count file format. It should contain the ",
           "following columns: 'peakIds', 'mz', 'rtMean', 'rtMin', 'rtMax' and all the ", 
           "sample code names concatenated with the suffix '_area'.")
    }
  } else if (is.data.frame(quantification_table)) {
    if (!all(c("peakIds", "mz", "rtMean", "rtMin", "rtMax", 
               metadata$SAMPLE_CODE) %in% names(quantification_table))) {
      stop("Wrong count file format. It should contain the ",
           "following columns: 'peakIds', 'mz', 'rtMean', 'rtMin', 'rtMax' and all the ", 
           "sample code names concatenated with the suffix '_area'.")
    }
  } else {
    stop("Wrong count file '", typeof(quantification_table),
         "' format. It should be a character with the path to the file with ",
         "the peak area count of not fragmented MS2 or a data.frame with the ",
         "following columns: 'peakIds', 'mz', 'rtMean', 'rtMin', 'rtMax' and all the ", 
         "sample code names concatenated with the suffix '_area'.")
  }
  
  if (is.na(mz_tol) || is.na(rt_tol) || is.null(mz_tol) || is.null(rt_tol) ||
      mz_tol < 0 || rt_tol < 0)
    stop("Wrong arguments 'mz_tol' and 'rt_tol' values. They should a numeric value greater or equal zero.")
  
  # filter only the samples that are in the metadata
  quantification_table <- quantification_table[,c("peakIds", "mz", "rtMean", 
                                                  "rtMin", "rtMax", metadata$SAMPLE_CODE)]
  if (any(duplicated(quantification_table$peakIds))) {
    warning("Duplicated peakIds's are present in the MS1_list_no_MS2.csv table and the first occurrences will be removed.")
    quantification_table <- quantification_table[!duplicated(quantification_table$peakIds),]  
  }
  
  # initialize indicators
  quantification_table$numJoins <- 0
  quantification_table$sumAreas <- rowSums(quantification_table[, metadata$SAMPLE_CODE])
  
  cat("\n** Cleanning counts and aggregating data of not fragmented MS1 m/z's **\n")
  t0 <- Sys.time()
  step_join <- 1
  num_joins_total <- 0
  # stores the number of joins by step and the total number of joins by cluster
  num_joins_last_step <- total_num_join_clusters <- quantification_table$numJoins
  nscans <- nrow(quantification_table)
  
  repeat
  {
    t1 <- Sys.time()
    cat("\n      * Step", step_join, "*\n")
    num_joins <- 0
    
    progress_joins <- unique(trunc(c(seq(from = 1, to = nscans, by = nscans/25), nscans)))
    n_progress <- length(progress_joins)
    cat(paste0("  |", paste0(rep(" ", n_progress), collapse = ""), "|\n  |", collapse = ""))
    # cat("  |", rep("", (n_progress-1)), "|\n  |")
    
    i <- 1
    while (i <= nscans) 
    {
      # print(i)
      # print(nscans)
      if (n_progress > 0 && i == progress_joins[[1]]) {
        progress_joins <- progress_joins[-1]
        n_progress <- n_progress - 1
        cat("=")
      }
      
      # get next cluster scan number and cluster info
      cluster <- quantification_table[i,]
      scan_num <- cluster[[1]]
      # print(cluster)
      
      
      # after step 1 only check clusters that changed
      if (step_join > 1 && cluster$numJoins == 0) 
      { 
        i <- i + 1
        next()
      }
    
      cluster_peak <- quantification_table[abs(quantification_table$mz - cluster$mz) <= mz_tol &
                                      ((quantification_table$rtMean >= cluster$rtMin - rt_tol & 
                                          quantification_table$rtMean <= cluster$rtMax + rt_tol) |
                                         (cluster$rtMean >= quantification_table$rtMin - rt_tol &
                                            cluster$rtMean <= quantification_table$rtMax + rt_tol)),] 
      # print(cluster_peak)
      if (nrow(cluster_peak) == 1) # just the current cluster, then go to next spectrum
      {
        i <- i + 1
        next()
      } 
      num_joins <- num_joins + nrow(cluster_peak) - 1
    
      # merge counts based on peak area then on spectra num
      cluster <- lapply(names(cluster_peak), merge_no_ms2_counts, cluster_peak)
      i_count <- match(cluster_peak$peakIds, quantification_table$peakIds) # get the cluster pos in the count table
      quantification_table[i,] <- cluster
      
      # remove merged row from count tables
      i_count <- i_count[i_count != i] # get the cluster_peak pos in the count table
      quantification_table <- quantification_table[-i_count,]
      
      num_joins_last_step <- num_joins_last_step[-i_count]
      total_num_join_clusters <- total_num_join_clusters[-i_count]
      
      nscans <- nscans - nrow(cluster_peak) + 1 # remove scans that were merged from the total count
      if (n_progress > 0 && progress_joins[n_progress] > nscans) # rm from progress the joined scans
      {
        cat("=")
        progress_joins <- progress_joins[-n_progress]
        n_progress <- n_progress - 1
      }
      
      # decrement number of rows preceding the current that were removed - prevent skipping rows
      i <- i + 1 - sum(i_count < i)
    }
    cat("|\n")
    
    # reset number of joins by m/z with the joined clusters of last step 
    quantification_table$numJoins <- quantification_table$numJoins - num_joins_last_step
    
    # reset number of joined clusters by m/z in this step
    num_joins_last_step <- quantification_table$numJoins
    # update number of total joins by cluster after this step
    total_num_join_clusters <- total_num_join_clusters + num_joins_last_step
    
    t2 <- Sys.time()
    cat("        * Joined", num_joins, "similar clusters in", 
            round(t2-t1, 2), units(t2-t1), "*\n")
    
    num_joins_total <- num_joins_total + num_joins
    
    # stop joining if no join was made in the last step
    if (num_joins == 0 || step_join == 10)
      break()
    
    step_join <- step_join + 1
  }
  
  # compute Blanks total and controls total and beds total
  blanks_code <- metadata[metadata$SAMPLE_TYPE == "blank", "SAMPLE_CODE"]
  if (length(blanks_code) > 0)
  {
    if (length(blanks_code) == 1) {
      quantification_table$BLANKS_TOTAL <- quantification_table[,  blanks_code][[1]]
    } else {
      quantification_table$BLANKS_TOTAL <- rowSums(quantification_table[, blanks_code])
    }
    mz_blank <- almost.unique(quantification_table[quantification_table$BLANKS_TOTAL > 0, 
                                                   "mz"][[1]], mz_tol)
    mz_blank_up <- mz_blank + mz_tol # upper bound
    mz_blank <- mz_blank - mz_tol # lower bound
  }
  controls_code <- metadata[metadata$SAMPLE_TYPE == "control", "SAMPLE_CODE"]
  if (length(controls_code) > 0)
  {
    if (length(controls_code) == 1)
    {
      quantification_table$CONTROLS_TOTAL <- quantification_table[, controls_code][[1]]
    } else {
      quantification_table$CONTROLS_TOTAL <- rowSums(quantification_table[, controls_code])
    }
    mz_control <- almost.unique(quantification_table[quantification_table$CONTROLS_TOTAL > 0, 
                                                     "mz"][[1]], mz_tol)
    mz_control_up <- mz_control + mz_tol # upper bound
    mz_control <- mz_control - mz_tol # lower bound
  }
  bed_controls_code <- metadata[metadata$SAMPLE_TYPE == "bed", "SAMPLE_CODE"]
  if (length(bed_controls_code) > 0)
  {
    if (length(bed_controls_code) == 1)
    {
      quantification_table$BEDS_TOTAL <- quantification_table[, bed_controls_code][[1]]
    } else {
      quantification_table$BEDS_TOTAL <- rowSums(quantification_table[, bed_controls_code])
    }
    mz_bed_control <- almost.unique(quantification_table[quantification_table$BEDS_TOTAL > 0, 
                                                         "mz"][[1]], mz_tol)
    mz_bed_control_up <- mz_bed_control + mz_tol # upper bound
    mz_bed_control <- mz_bed_control - mz_tol # lower bound
  }
  
  hits_code <- metadata[metadata$SAMPLE_TYPE == "hit", "SAMPLE_CODE"]
  if (length(hits_code) > 0) {
    mz_hits <- lapply(hits_code, function(y){
      almost.unique(unlist(quantification_table[quantification_table[,y] > 0, "mz"], use.names = FALSE))})
  }
  
  # add BFLAG and CFLAG and DESREPLICATION and HFLAG and BEDFLAG
  if (length(bed_controls_code) > 0 || length(controls_code) > 0 || length(blanks_code) > 0 || length(hits_code) > 0)
  {
    quantification_table <- dplyr::bind_cols(
      quantification_table, 
      dplyr::bind_rows(lapply(seq_along(quantification_table$mz), 
                       function(i)
                       {
                         x <- quantification_table$mz[[i]]
                         
                         flagColumns = list()
                         
                         if (length(hits_code) > 0) {
                           # get the hit samples in which the mass x appears
                           hflag <- unlist(sapply(seq_along(mz_hits), function(j){
                             if (any(x >= mz_hits[[j]] - mz_tol & x <= mz_hits[[j]] + mz_tol))
                               return(hits_code[[j]])
                             else
                               return(NULL)
                           }))
                           
                           if (is.null(hflag))
                             hflag <- NA
                           else
                             hflag <- paste(hflag, collapse = ";")
                           
                           if (length(hits_code) > 0)
                             desrep <- paste(hits_code[which(quantification_table[i, hits_code] > 0)], collapse = ";")
                           else
                             desrep <- NA
                           
                           flagColumns <- c(flagColumns, HFLAG = hflag, 
                                            DESREPLICATION = ifelse(desrep != "", 
                                                                    desrep, NA))
                         }
                         if (length(blanks_code) > 0)
                           flagColumns <- c(flagColumns, BFLAG = any(x >= mz_blank & x <= mz_blank_up))
                         if (length(controls_code) > 0)
                           flagColumns <- c(flagColumns, CFLAG = any(x >= mz_control & x <= mz_control_up))
                         if (length(bed_controls_code) > 0)
                           flagColumns <- c(flagColumns, BEDFLAG = any(x >= mz_bed_control & x <= mz_bed_control_up))
                         
                         flagColumns
                       })))
  }
  
  tf <- Sys.time()
  cat("\n * Done reducing", num_joins_total, "MS1 peaks in", 
          round(tf-t0, 2), units(tf-t0), " *\n\n")
  
  # order columns using the metadata SAMPLE_CODE order
  quantification_table <- quantification_table[, c("peakIds","mz", "rtMean", "rtMin", "rtMax", 
                                 "sumAreas", metadata$SAMPLE_CODE, 
                                 c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                   "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG")[
                                     c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                       "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG") %in% 
                                       names(quantification_table)])]
  # add ids
  quantification_table <- cbind(msclusterID = 1:nrow(quantification_table),
                                quantification_table)
  # set mz column to match mscluster output
  names(quantification_table)[names(quantification_table) == "mz"] <- "mzConsensus"
  
  # check if peak area count of MS2 file exists and use it to enumerate the mz's ids uniquely
  if (!missing(output_path)) {
    if (file.exists(sub("_MS1", "", output_path))) {
      # increment 1 to skip last msclusterID
      msclusterID_i <- length(readr::count_fields(sub("_MS1", "", output_path), 
                                           skip = 1, tokenizer = readr::tokenizer_csv()))+1
      quantification_table$msclusterID <- msclusterID_i:(msclusterID_i+nrow(quantification_table)-1)
    }
    
    readr::write_csv(quantification_table, path = output_path)
  }
  
  invisible(quantification_table)
}

clean_quantification(quantification_table_path, metadata_path, output_path, mz_tol, rt_tol)
