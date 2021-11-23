suppressPackageStartupMessages(library(dplyr))
##
# Step 4 - quantifications from the subclustering results
##

library(purrr)
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

# from package bazar
almost.unique <- function(x, tol = sqrt(.Machine$double.eps))
{
  y <- round(x/tol, 0)
  d <- duplicated(y)
  x[!d]
}

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Three arguments must be supplied to count number spectra by input file:\n",
       " 1 - Path to the batches output folder;\n",
       " 2 - Batch name;\n",
       " 3 - Is final step? 0 if is a batch and 1 if is the final integration step.\n",
       " 4 - if (3) is 1, then metadata path should be provided\n",
       " 5 - if (3) is 1, then processed_data_path should be provided\n",
       " 6 - if (3) is 1, then mz tolerance should be provided\n",
       call.=FALSE)
} else {
  outFilesDir <- file.path(args[[1]])
  # validate input
  if (!dir.exists(outFilesDir))
  {
    stop("The path to the batches output folder '", outFilesDir,
         "' do not exists. Provide a valid path to where the batches output are located.")
  }
  outFilesDir <- normalizePath(outFilesDir)
  outputName <- args[[2]]
  isLast <- args[[3]]
  
  # count spectra of each subcluster
  if (isLast == "0") # it is not the final step
  {
    batchesNum <- sort(as.numeric(sub(pattern = paste0(outputName,"_"), replacement = "", 
                                      x= list.files(outFilesDir)[startsWith(list.files(outFilesDir), paste0(outputName,"_"))])))
    
    # get the count of spectra of the subbatches
    countSubclusters <- lapply(paste0(outputName,"_", batchesNum), function(x)
    {
      countFile <- read.csv(file.path(outFilesDir, x, "count_tables", paste0(x, "_spectra.csv")), 
                            stringsAsFactors = F, comment.char = "")
      countFile <- countFile[, c("msclusterID", "numSpectra", "peakIds", "scans",
                                 names(countFile)[!(names(countFile) %in% 
                                                      c("mzConsensus", "msclusterID", "peakIds", "scans", "numSpectra",
                                                        "rtMean", "rtMin", "rtMax", "sumInts"))])]
      countFile
    })
    names(countSubclusters) <- batchesNum
  } else { # its the final step, integrating batches
    path_batch_metadata <- file.path(args[[4]])
    if (!file.exists(path_batch_metadata))
    {
      stop("The CSV batch metadata file '", path_batch_metadata, 
           "' do not exists. Provide a valid path to where the metadata file is located.")
    }
    path_batch_metadata <- normalizePath(path_batch_metadata)
    
    batch_metadata <- readMetadataTable(path_batch_metadata)
    n_batches <- unique(batch_metadata$DATA_COLLECTION_BATCH)
    
    if (length(n_batches) > 1) # more than one batch, use only the final batches output
    {
      finalBatchesNum <- sort(as.numeric(sub(pattern = "B_", replacement = "", 
                                             x = list.files(outFilesDir)[grepl("B_[0-9]+$", list.files(outFilesDir))])))
      
      # get the count of peak area of the subbatches
      countSubclusters <- lapply(paste0("B_", finalBatchesNum), function(x)
      {
        countFile <- read.csv(file.path(outFilesDir, x, "count_tables", paste0(x, "_spectra.csv")), 
                              stringsAsFactors = F, comment.char = "")
        countFile <- countFile[, c("msclusterID", "numSpectra", "peakIds", "scans",
                                   names(countFile)[!(names(countFile) %in% 
                                                        c("mzConsensus", "msclusterID", "peakIds", "scans", "numSpectra",
                                                          "rtMean", "rtMin", "rtMax", "sumInts"))])]
        countFile
      })
      names(countSubclusters) <- finalBatchesNum
    } else { # just one batch, use all sub batches output
      batchesNum <- sort(as.numeric(sub(pattern = "B_1_", replacement = "", 
                                        x= list.files(outFilesDir)[grepl("^B_[0-9]+", list.files(outFilesDir))])))
      # get the subbatch count of spectra
      countSubclusters <- lapply(paste0("B_1_", batchesNum), function(x)
      {
        countFile <- read.csv(file.path(outFilesDir, x, "count_tables", paste0(x, "_spectra.csv")), 
                              stringsAsFactors = F, comment.char = "")
        countFile <- countFile[, c("msclusterID", "numSpectra", "peakIds", "scans",
                                   names(countFile)[!(names(countFile) %in% 
                                                        c("mzConsensus", "msclusterID", "peakIds", "scans", "numSpectra",
                                                          "rtMean", "rtMin", "rtMax", "sumInts"))])]
        countFile
      })
      names(countSubclusters) <- batchesNum
    }
    if (length(countSubclusters) == 0) ## something went wrong when locating the subbatches
    {
      stop("No sub batches for '", outputName,
           "' were found. Maybe calling script with wrong parameters.")
    }
  }
} 


clustFilesDir <- file.path(outFilesDir, outputName, "clust")
if (!dir.exists(clustFilesDir))
{
  stop("The path to the batch .clust folder '", clustFilesDir,
       "' do not exists. Provide a valid path to where the batches output are located.")
}
clustFiles  <- list.files(clustFilesDir, pattern = "\\.clust")

clusterList <- bind_rows(lapply(clustFiles, function(clustFileName)
{
  cat("- Start scanning file: ", clustFileName)
  clust_file <- read.csv(file = file.path(clustFilesDir, clustFileName), 
                         comment.char = "", stringsAsFactors = F)
  
  # get peak and area count by cluster
  clustersMembers <- bind_rows(lapply(unique(clust_file$clustId), function(id)
  {
    # print(id)
    clust_members <- clust_file[clust_file$clustId == id,]
    
    # aggregate counts for the same fileIndex
    countSamples <- bind_rows(lapply(unique(clust_members$fileIndex), function(subCluster)
    {
      ids <- clust_members[clust_members$fileIndex == subCluster, "scanNumber"]
      return(c(list(peakIds = paste0(countSubclusters[[subCluster+1]][
        countSubclusters[[subCluster+1]]$msclusterID %in% ids, 3], 
        collapse = ";"), # collapse peaksIds
        scans = paste0(countSubclusters[[subCluster+1]][
          countSubclusters[[subCluster+1]]$msclusterID %in% ids, 4], collapse = ";")),
        colSums(countSubclusters[[subCluster+1]][   
          countSubclusters[[subCluster+1]]$msclusterID %in% ids, c(-3,-4)], na.rm = TRUE)))  # collapse scans
    }))
    
    # if peakId is not missing, aggregate them
    if (any(!is.na(countSamples$peakIds))) { 
      peakIds <- paste0(unique(unlist(strsplit(countSamples$peakIds[!is.na(countSamples$peakIds)], 
                                             split = ";"))), collapse = ";")
    } else {
      peakIds <- NA
    }
    # get scans
    scans <- paste0(countSamples$scans[!is.na(countSamples$scans)], collapse = ";")
    
    countSamples <- countSamples[, !(names(countSamples) %in% 
                                       c("msclusterID", "peakIds", "scans"))]
    # sum counts from all file indexes
    if (NROW(countSamples) > 1)
    {
      countSamples <- colSums(countSamples, na.rm = TRUE)
    }
    
    numSpectra <- gregexpr(";", scans)[[1]] # count number of ;
    if (numSpectra[[1]] == -1) {
      numSpectra <- 1 # no ; means only one spectra
    } else {
      numSpectra <- length(numSpectra)+1 # number of scans = number of ; + 1
    }
    if (numSpectra != countSamples[["numSpectra"]]) {
      stop("Something went wrong in the counting. Wrong number of spectra for ",
           "msclusterId ", clust_members$clustId[[1]])
    }
    countSamples <- countSamples[!(names(countSamples) %in% c("numSpectra"))]
    countSamples[is.na(countSamples)] <- 0
    
    return(c(list(msclusterID = clust_members$clustId[[1]],
                  numSpectra = numSpectra,
                  mzConsensus = clust_members$mzConsensus[[1]],
                  rtMean = clust_members$rtMean[[1]],
                  rtMin = clust_members$rtMin[[1]],
                  rtMax = clust_members$rtMax[[1]],
                  peakIds = peakIds,
                  scans = scans,
                  sumInts = clust_members$precursorInt[[1]]),
             countSamples))
  }))
  
  rm(clust_file)
  
  cat(" DONE!\n")
  return(clustersMembers)
})) %>% arrange(as.numeric(msclusterID))
# set the missing counts to zero
peakIds <- clusterList$peakIds
scans <- clusterList$scans
clusterList[is.na(clusterList)] <- 0
clusterList$peakIds <- peakIds
clusterList$scans <- scans
rm(peakIds, scans)

# check if any spectra is being counted twice
scans_list <- strsplit(paste(clusterList$scans, collapse = ";"), ";")[[1]]
if (any(duplicated(scans_list))) {
  stop("Error in the sample counts. The following SCANS are being counted more than once: ",
       paste(scans_list[duplicated(scans_list)], collapse = ","))
} 
rm(scans_list)

# order columns
clusterList <- clusterList[, c("msclusterID","numSpectra", 
                               "mzConsensus", "rtMean", "rtMin",
                               "rtMax", "peakIds", "scans", "sumInts",
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_spectra")])]

## check count sum
if (sum(clusterList$numSpectra) != 
    sum(clusterList[,!(names(clusterList) %in% c("msclusterID","numSpectra", 
                                                 "mzConsensus", "rtMean", "rtMin", 
                                                 "rtMax", "peakIds", "scans", "sumInts"))])) {
  stop("Spectra count sum check failed. Could not count number of spectra.")
}

# if it is the final step, compute indicators and peak area
if (isLast == "1") 
{
  source(file.path(script_path(), "count_peak_area.R"))
  processed_data_path <- file.path(args[[5]])
  if (!dir.exists(processed_data_path))
  {
    stop("The processed data folder '", processed_data_path, 
         "' do not exists. Provide a valid path to where the original MGFs are. ",
         " The peak area count could not be computed")
  }
  # compute peak area count and the base peak int
  clusterArea <- compute_peak_area(processed_data_path, 
                                   clusterList$msclusterID,
                                   lapply(clusterList$scans, function(x) strsplit(x, ";")[[1]]),
                                   lapply(clusterList$peakIds, function(x) strsplit(x, ";")[[1]]),
                                   batch_metadata)
  clusterArea <- bind_cols(clusterList[, c("msclusterID","numSpectra", 
                                           "mzConsensus", "rtMean", "rtMin",
                                           "rtMax", "peakIds", "scans", "sumInts")],
                           clusterArea)
  clusterList$basePeakInt <- clusterArea$basePeakInt
  
  # compute indicators
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
    clusterList <- bind_cols(clusterList, 
                             bind_rows(lapply(seq_along(clusterList$mzConsensus), 
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
  
  # order columns using the metadata SAMPLE_CODE order
  clusterList <- clusterList[, c("msclusterID","numSpectra", 
                                 "mzConsensus", "rtMean", 
                                 "rtMin", "rtMax", "peakIds", "scans", 
                                 "sumInts", "basePeakInt",
                                 paste0(batch_metadata$SAMPLE_CODE,"_spectra"), 
                                 c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                   "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG")[
                                     c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                       "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG") %in% 
                                       names(clusterList)])]
  clusterArea <- clusterArea[, c("msclusterID","numSpectra",
                                 "mzConsensus", "rtMean",
                                 "rtMin", "rtMax", "peakIds", "scans", 
                                 "sumInts", "basePeakInt",
                                 paste0(batch_metadata$SAMPLE_CODE,"_area"),
                                 c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                   "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG")[
                                     c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                       "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG") %in%
                                       names(clusterArea)])]
  
  dir.create(file.path(outFilesDir, outputName, "count_tables"), showWarnings = FALSE)
  write.csv(clusterList, 
            file = file.path(outFilesDir, outputName, "count_tables", paste0(outputName, "_spectra.csv")), 
            row.names = FALSE)
  write.csv(clusterArea,
            file = file.path(outFilesDir, outputName, "count_tables", paste0(outputName, "_peak_area.csv")),
            row.names = FALSE)
} else {
  dir.create(file.path(outFilesDir, outputName, "count_tables"), showWarnings = FALSE)
  write.csv(clusterList, 
            file = file.path(outFilesDir, outputName, "count_tables", paste0(outputName, "_spectra.csv")), 
            row.names = FALSE)
}

# rm spec file
# rm spec from root folder
invisible(file.copy(from = file.path(outFilesDir, outputName, paste0(outputName, "_0_spec_list.txt")), 
          to = file.path(outFilesDir, outputName, 'clust'),
          recursive = TRUE, copy.mode = TRUE, copy.date = FALSE))
unlink(file.path(outFilesDir, outputName, paste0(outputName, "_0_spec_list.txt")))