##
# Step 4 - quantifications from the pre-processed files
##

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
source(file.path(script_path(), "read_metadata_table.R"))

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Three arguments must be supplied to count number spectra by input file:\n", 
       " 1 - Path to the MSCluster output folder;\n",
       " 2 - MSCluster output prefix name;\n", 
       " 3 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks;\n",
       " 4 - Is a single batch without subclusters? 0 if there are sub batches and 1 if is a single run;\n", 
       " 5 - if (4) is 1, then processed_data_path should be provided\n",
       " 6 - if (4) is 1, then mz tolerance should be provided\n",
       call.=FALSE)
} else {
  outFilesDir <- file.path(args[[1]])
  # validate input
  if (!dir.exists(outFilesDir))
  {
    stop("The MSCluster output folder '", outFilesDir, 
         "' do not exists. Provide a valid path to where the clustered data is located.")
  }
  outFilesDir <- normalizePath(outFilesDir)
  
  outputName <- args[[2]]
  path_batch_metadata <- file.path(args[[3]])
  if (!file.exists(path_batch_metadata))
  {
    stop("The CSV batch metadata file '", path_batch_metadata, 
         "' do not exists. Provide a valid path to where the batch metadata is located.")
  }
  path_batch_metadata <- normalizePath(path_batch_metadata)
  
  singleBatch <- args[[4]] # is single batch??
}

clustFilesDir <- file.path(outFilesDir, "clust")
if (!dir.exists(clustFilesDir))
{
  stop("The MSCluster clust output folder '", clustFilesDir, 
       "' do not exists. Provide a valid path to where the clustered data is located.")
}
clustFiles    <- list.files(clustFilesDir, pattern = "\\.clust")
if (length(clustFiles) == 0)
{
  stop("The MSCluster clust output folder '", clustFilesDir, 
       "' is empty. Provide a valid path to where the clustered data is located.")
}

# get list of files used in the MSCluster run by index order 0..length(files)
filesIndex <- data.frame(FILENAME = basename(read.table(
  file.path(outFilesDir, paste0(outputName, "_0_spec_list.txt")), header = F, 
  stringsAsFactors = F, comment.char = "")[,1]), stringsAsFactors = F)

batch_metadata <- readMetadataTable(path_batch_metadata)

# using processed data, transform filenames
if (all(endsWith(filesIndex$FILENAME, "_peak_info.mgf")))
{ 
  batch_metadata$FILENAME <- paste0(batch_metadata$SAMPLE_CODE, 
                                              "_peak_info.mgf")
}
filesIndex <- left_join(x = filesIndex, y = batch_metadata, by = "FILENAME")

clusterList <- bind_rows(lapply(clustFiles, function(clustFileName)
{
  cat("- Start scanning file: ", clustFileName)
  clust_file <- read.csv(file = file.path(clustFilesDir, clustFileName), 
                         comment.char = "", stringsAsFactors = F)
  
  # get peak and area count by cluster
  clustersMembers <- bind_rows(lapply(unique(clust_file$clustId), function(id)
    {
      clust_members <- clust_file[clust_file$clustId == id,]
      
      areaSamples <- countSamples <- tabulate(clust_members$fileIndex+1, nbins = nrow(filesIndex))
      areaSamples[which(areaSamples > 0)] <- sapply((which(areaSamples > 0)-1), function(x) {
        sum(unique(trunc(clust_members$peakArea[clust_members$fileIndex == x])))
      })
      names(countSamples) <- paste0(filesIndex$SAMPLE_CODE, "_spectra")
      names(areaSamples) <- paste0(filesIndex$SAMPLE_CODE, "_area")
      
      if (all(is.na(clust_members$peakId)))
        peakIds <- NA
      else
        peakIds <- paste0(unique(clust_members$peakId[!is.na(clust_members$peakId)]), collapse = ";")
      
      scans <- paste0(clust_members$scanNumber, "_", 
                      filesIndex$SAMPLE_CODE[clust_members$fileIndex+1], collapse = ";")
      
      return(c(list(msclusterID = clust_members$clustId[[1]],
                    numSpectra = clust_members$clustSize[[1]],
                    mzConsensus = clust_members$mzConsensus[[1]],
                    rtMean = clust_members$rtMean[[1]],
                    rtMin = clust_members$rtMin[[1]],
                    rtMax = clust_members$rtMax[[1]],
                    sumInts = clust_members$precursorInt[[1]],
                    peakIds = peakIds,
                    scans = scans),
                    countSamples, areaSamples))
    }))
  
  rm(clust_file)

  cat(" DONE!\n")
  return(clustersMembers)
})) %>% arrange(as.numeric(msclusterID))

# check if any spectra is being counted twice
scans_list <- strsplit(paste(clusterList$scans, collapse = ";"), ";")[[1]]
if (any(duplicated(scans_list))) {
  stop("Error in the sample counts. The following SCANS are being counted more than once: ",
       paste(scans_list[duplicated(scans_list)], collapse = ","))
} 
rm(scans_list)

# separate the cluster count of spectra from the cluster peak area
clusterArea <- clusterList[, c("msclusterID","numSpectra", 
                               "mzConsensus", "rtMean", "rtMin", 
                               "rtMax", "peakIds", "scans", "sumInts", 
                               names(clusterList)[endsWith(names(clusterList), 
                                                             "_area")])]

clusterList <- clusterList[, c("msclusterID","numSpectra", 
                               "mzConsensus", "rtMean", "rtMin", 
                               "rtMax", "peakIds", "scans", "sumInts", 
                               names(clusterList)[endsWith(names(clusterList), 
                                                           "_spectra")])]

## check count sum
if (sum(clusterList$numSpectra) != 
    sum(clusterList[,!(names(clusterList) %in% c("msclusterID","numSpectra", 
                                                 "mzConsensus", "rtMean", "rtMin", 
                                                 "rtMax", "peakIds", "scans", "sumInts"))]))
  stop("Spectra count sum check failed. ",
       "The number of spectra in the column 'numSpectra' is different from the ",
       "sum of spectra counts by sample.")

if (singleBatch == "1") # only one batch, compute sum of blanks
{
  source(file.path(script_path(), "count_peak_area.R"))
  processed_data_path <- file.path(args[[5]])
  if (!dir.exists(processed_data_path))
  {
    stop("The processed data folder '", processed_data_path, 
         "' do not exists. Provide a valid path to where the original MGFs are. ",
         " The peak area count could not be computed")
  }
  
  # from package bazar
  almost.unique <- function(x, tol = sqrt(.Machine$double.eps))
  {
    y <- round(x/tol, 0)
    d <- duplicated(y)
    x[!d]
  }
  mz_tol <- as.numeric(args[[6]])
  
  blanksCode <- batch_metadata[batch_metadata$SAMPLE_TYPE == "blank", "SAMPLE_CODE"]
  controlsCode <- batch_metadata[batch_metadata$SAMPLE_TYPE =="control", "SAMPLE_CODE"]
  bedControlsCode <- batch_metadata[batch_metadata$SAMPLE_TYPE == "bed", "SAMPLE_CODE"]
  hitsCode <- batch_metadata[batch_metadata$SAMPLE_TYPE == "hit", "SAMPLE_CODE"]
  
  if (length(hitsCode) > 0) {
    hitsCode <- paste0(hitsCode, "_spectra")
    mz_hits <- lapply(hitsCode, function(y){
      almost.unique(unlist(clusterList[clusterList[,y] > 0, "mzConsensus"], use.names = FALSE))})
  }
  if (length(blanksCode) > 0) {
    clusterList <- mutate(clusterList, 
                          BLANKS_TOTAL = rowSums(select(clusterList, paste0(blanksCode,"_spectra"))))
    clusterArea <- mutate(clusterArea, 
                          BLANKS_TOTAL = rowSums(select(clusterArea, paste0(blanksCode,"_area"))))
    
    mz_blank <- almost.unique(clusterList[clusterList$BLANKS_TOTAL > 0, "mzConsensus"][[1]], mz_tol)
    mz_blank_up <- mz_blank + mz_tol # upper bound
    mz_blank <- mz_blank - mz_tol # lower bound
  }
  if (length(controlsCode) > 0) {
    clusterList <- mutate(clusterList, 
                          CONTROLS_TOTAL = rowSums(select(clusterList, paste0(controlsCode,"_spectra"))))
    clusterArea <- mutate(clusterArea, 
                          CONTROLS_TOTAL = rowSums(select(clusterArea, paste0(controlsCode,"_area"))))
    
    mz_control <- almost.unique(clusterList[clusterList$CONTROLS_TOTAL > 0, "mzConsensus"][[1]], mz_tol)
    mz_control_up <- mz_control + mz_tol # upper bound
    mz_control <- mz_control - mz_tol # lower bound
  }
  if (length(bedControlsCode) > 0) {
    clusterList <- mutate(clusterList, 
                          BEDS_TOTAL = rowSums(select(clusterList, paste0(bedControlsCode,"_spectra"))))
    clusterArea <- mutate(clusterArea, 
                          BEDS_TOTAL = rowSums(select(clusterArea, paste0(bedControlsCode,"_area"))))
    
    mz_bed_control <- almost.unique(clusterList[clusterList$BEDS_TOTAL > 0, "mzConsensus"][[1]], mz_tol)
    mz_bed_control_up <- mz_bed_control + mz_tol # upper bound
    mz_bed_control <- mz_bed_control - mz_tol # lower bound
  }
  
  # add BFLAG and CFLAG and DESREPLICATION and HFLAG and BEDFLAG
  if (length(bedControlsCode) > 0 || length(controlsCode) > 0 || length(blanksCode) > 0 || length(hitsCode) > 0)
  {
    clusterList <- bind_cols(clusterList, bind_rows(lapply(seq_along(clusterList$mzConsensus), 
      function(i)
      {
       x <- clusterList$mzConsensus[[i]]
       
       flagColumns = list()
       
       if (length(hitsCode) > 0) {
         # get the hit samples in which the mass x appears
         hflag <- unlist(sapply(seq_along(mz_hits), function(j){
           if (any(x >= mz_hits[[j]] - mz_tol & x <= mz_hits[[j]] + mz_tol))
             return(hitsCode[[j]])
           else
             return(NULL)
         }))
         
         if (is.null(hflag))
           hflag <- NA
         else
           hflag <- paste(hflag, collapse = ";")
         
         if (length(hitsCode) > 0)
           desrep <- paste(hitsCode[which(clusterList[i, hitsCode] > 0)], collapse = ";")
         else
           desrep <- NA
         
         flagColumns <- c(flagColumns, HFLAG = hflag, 
                          DESREPLICATION = ifelse(desrep != "", 
                                                  desrep, NA))
       }
       if (length(blanksCode) > 0)
         flagColumns <- c(flagColumns, BFLAG = any(x >= mz_blank & x <= mz_blank_up))
       if (length(controlsCode) > 0)
         flagColumns <- c(flagColumns, CFLAG = any(x >= mz_control & x <= mz_control_up))
       if (length(bedControlsCode) > 0)
         flagColumns <- c(flagColumns, BEDFLAG = any(x >= mz_bed_control & x <= mz_bed_control_up))
       
       flagColumns
      })))
    clusterArea <- bind_cols(clusterArea, clusterList[, 
      names(clusterList)[names(clusterList) %in% 
                           c("BFLAG", "CFLAG", "BEDFLAG", "HFLAG", "DESREPLICATION")]])
  }
  
  # compute peak area count and the base peak int
  clusterArea_basePeakInt <- compute_peak_area(processed_data_path, 
                                               clusterList$msclusterID,
                                               lapply(clusterList$scans, function(x) strsplit(x, ";")[[1]]),
                                               lapply(clusterList$peakIds, function(x) strsplit(x, ";")[[1]]),
                                               batch_metadata)
  clusterArea[,names(clusterArea_basePeakInt)] <- clusterArea_basePeakInt
  clusterList$basePeakInt <- clusterArea$basePeakInt
  
  # order columns using the metadata SAMPLE_CODE order
  clusterList <- clusterList[, c("msclusterID","numSpectra", 
                                 "mzConsensus", "rtMean", 
                                 "rtMin", "rtMax", "peakIds", "scans", "sumInts",
                                 "basePeakInt",
                                 paste0(batch_metadata$SAMPLE_CODE,"_spectra"), 
                                 c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                   "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG")[
                                     c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                       "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG") %in% 
                                       names(clusterList)])]
  clusterArea <- clusterArea[, c("msclusterID","numSpectra", 
                                 "mzConsensus", "rtMean", 
                                 "rtMin", "rtMax", "peakIds", "scans", "sumInts",
                                 "basePeakInt",
                                 paste0(batch_metadata$SAMPLE_CODE,"_area"), 
                                 c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                   "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG")[
                                     c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                       "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG") %in% 
                                       names(clusterArea)])]
}

dir.create(file.path(outFilesDir, "count_tables"), showWarnings = FALSE)
write.csv(clusterList, file = file.path(outFilesDir, "count_tables", paste0(outputName, "_spectra.csv")), 
          row.names = FALSE)
write.csv(clusterArea, file = file.path(outFilesDir, "count_tables", paste0(outputName, "_peak_area.csv")), 
          row.names = FALSE)

# rm spec from root folder
invisible(file.copy(from = file.path(outFilesDir, paste0(outputName, "_0_spec_list.txt")), 
          to = file.path(outFilesDir, 'clust'), 
          recursive = TRUE, copy.mode = TRUE, copy.date = FALSE))
unlink(file.path(outFilesDir, paste0(outputName, "_0_spec_list.txt")))
