suppressPackageStartupMessages(library(dplyr))
##
# Step 4 - quantifications from the subjobs clean results of joining jobs
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
       " 1 - Path to the joined jobs final output folder, inside the outs dir, named with the output name;\n",
       " 2 - Output name;\n",
       " 3 - Path to the jobs data folder where the jobs results are located.\n",
       " 4 - The metadata path with the joining jobs\n",
       " 5 - The mz tolerance - to compute flags indicators \n",
       call.=FALSE)
} else {
  output_path <- file.path(args[[1]])
  # validate input
  if (!dir.exists(output_path))
  {
    stop("The path to the joined jobs output folder '", output_path,
         "' do not exists. Provide a valid path to where the joined jobs result is located.")
  }
  output_path <- normalizePath(output_path)
  output_name <- args[[2]]
  
  path_jobs_data <- file.path(args[[3]])
  if (!dir.exists(path_jobs_data))
  {
    stop("The path to the jobs data folder '", path_jobs_data,
         "' do not exists. Provide a valid path to where the jobs results are located.")
  }
  path_jobs_data <- normalizePath(path_jobs_data)
  
  path_jobs_metadata <- file.path(args[[4]])
  if (!file.exists(path_jobs_metadata))
  {
    stop("The CSV join jobs metadata file '", path_jobs_metadata, 
         "' do not exists. Provide a valid path to where the metadata file is located.")
  }
  path_jobs_metadata <- normalizePath(path_jobs_metadata)
  
  metadata <- readMetadataTableJoinJobs(path_jobs_metadata, path_jobs_data)
  
  mz_tol <- as.numeric(args[[5]])
} 

# get the count of spectra of the subjobs
# remove columns with mz, rt, sumInts  or flags columns or
# peaks list and int in string or ionization variants annotation cols 
# or tremolo result
# TODO see if these columns are necessary - or maxarea, numjoins, meanInt, blankdist and protonated cols
count_subjobs <- lapply(metadata$JOBNAME, function(x)
{
  jobcode <- metadata[metadata$JOBNAME == x,"JOB_CODE"]
  # retrieve the spectra columns
  countFile <- read.csv(file.path(path_jobs_data, x, "outs", x, "count_tables", 
                                  "clean", paste0(x, "_spectra_clean_ann.csv")), 
                        stringsAsFactors = F, comment.char = "")
  # create the joinedJobsIDs column equals the msclusterID_JobCode if it does not exists
  if (is.null(countFile$joinedJobsIDs)) 
  {
    countFile$joinedJobsIDs <- paste0(countFile$msclusterID, "_", jobcode)
  }
  # TODO filter only the columns that will be used, ignore the rest - safer than removing the not wanted ones, 
  # to exclude any other column that the user may create
  countFile <- countFile[, c("msclusterID", "numSpectra", "peakIds", "scans", "joinedJobsIDs", "basePeakInt",
                             names(countFile)[
                               !((names(countFile) %in%
                                    c("mzConsensus", "msclusterID", "peakIds", "scans", "joinedJobsIDs", 
                                      "basePeakInt", "joinedIDs", 
                                      "numSpectra", "rtMean", "rtMin", "rtMax", "sumInts",
                                      "BFLAG", "CFLAG", "BEDFLAG", "HFLAG", "DESREPLICATION",
                                      "BLANK_DIST","maxArea","numJoins", "meanInt",
                                      "protonated_representative", "protonated_mzError_sum",
                                      "protonated_rtError_sum",
                                      "peaksList", "peaksInt", "adducts", "isotopes",
                                      "dimers", "multiCharges", "fragments", "analogs")) | 
                                   startsWith(names(countFile), "tremolo_") | 
                                   startsWith(names(countFile), "gnps_"))])]
  # read area table and concatenate _area columns
  countFile_area <- read.csv(file.path(path_jobs_data, x, "outs", x, "count_tables", 
                                  "clean", paste0(x, "_peak_area_clean_ann.csv")), 
                        stringsAsFactors = F, comment.char = "")
  # retrieve the peak area and total columns
  countFile_area <- countFile_area[,c("msclusterID",
                                      names(countFile_area)[endsWith(names(countFile_area), "_area") | 
                                                              endsWith(names(countFile_area), "_TOTAL")])]
  # rename the area TOTAL columns with a suffix equals _area, to separate them at the end
  names(countFile_area)[endsWith(names(countFile_area), "_TOTAL")] <- paste0(names(countFile_area)[endsWith(names(countFile_area), "_TOTAL")], "_area")
  countFile <- left_join(countFile, countFile_area, by="msclusterID")
  rm(countFile_area)
  countFile
})
names(count_subjobs) <- metadata$JOB_CODE

if (length(count_subjobs) == 0) ## something went wrong when locating the subjobs
{
  stop("No sub job for '", output_name,
       "' were found. Maybe calling script with wrong parameters.")
}

clustFilesDir <- file.path(output_path, "clust")
if (!dir.exists(clustFilesDir))
{
  stop("The path to the joined jobs 'clust' folder '", clustFilesDir,
       "' do not exists. Provide a valid path to where the joined jobs output are located.")
}
clustFiles  <- list.files(clustFilesDir, pattern = "\\.clust")

# aggregate the spectra from the cluster lists
clusterList <- bind_rows(lapply(clustFiles, function(clustFileName)
{
  cat("- Start scanning file: ", clustFileName)
  clust_file <- read.csv(file = file.path(clustFilesDir, clustFileName), 
                         comment.char = "", stringsAsFactors = F)
  
  # get peak and area count by cluster
  clustersMembers <- bind_rows(lapply(unique(clust_file$clustId), function(id)
  {
    #cat("\n- Start counting clustId: ", id)
    clust_members <- clust_file[clust_file$clustId == id,]
    
    # aggregate counts for the same fileIndex
    countJobs <- bind_rows(lapply(unique(clust_members$fileIndex), function(subJob)
    {
      #cat("\n- Start counting subJob: ", subJob)
      ids <- clust_members[clust_members$fileIndex == subJob, "scanNumber"]
      return(c(list(peakIds = paste0(count_subjobs[[subJob+1]][
        count_subjobs[[subJob+1]]$msclusterID %in% ids, 3], 
        collapse = ";"), # collapse peaksIds
        scans = paste0(count_subjobs[[subJob+1]][
          count_subjobs[[subJob+1]]$msclusterID %in% ids, 4], collapse = ";"), # collapse scans
        joinedJobsIDs = paste0(count_subjobs[[subJob+1]][
          count_subjobs[[subJob+1]]$msclusterID %in% ids, 5], collapse = ";"), # collapse joinedJobsIDs
        basePeakInt = max(count_subjobs[[subJob+1]][
          count_subjobs[[subJob+1]]$msclusterID %in% ids, 6])),  # max basePeakInt
        # sum the remaining columns - spectra counts
        colSums(count_subjobs[[subJob+1]][   
          count_subjobs[[subJob+1]]$msclusterID %in% ids, c(-3,-4,-5,-6)], na.rm = TRUE)))  
    }))
    
    # if peakId is not missing, aggregate them
    if (any(!is.na(countJobs$peakIds))) { 
      peakIds <- paste0(unique(unlist(strsplit(countJobs$peakIds[!is.na(countJobs$peakIds)], 
                                             split = ";"))), collapse = ";")
    } else {
      peakIds <- NA
    }
    # get scans and aggregate them
    scans <- paste0(countJobs$scans[!is.na(countJobs$scans)], collapse = ";")
    # aggregate the unique joinedJobsIDs
    joinedJobsIDs <- paste0(unique(unlist(strsplit(countJobs$joinedJobsIDs[!is.na(countJobs$joinedJobsIDs)], 
                                             split = ";"))), collapse = ";")
    
    countJobs <- countJobs[, !(names(countJobs) %in% 
                               c("msclusterID", "peakIds", "scans", "joinedJobsIDs"))]
    # sum counts from all file indexes
    if (NROW(countJobs) > 1)
    {
      countJobs <- colSums(countJobs, na.rm = TRUE)
    }
    
    numSpectra <- gregexpr(";", scans)[[1]] # count number of ;
    if (numSpectra[[1]] == -1) {
      numSpectra <- 1 # no ; means only one spectra
    } else {
      numSpectra <- length(numSpectra)+1 # number of scans = number of ; + 1
    }
    if (numSpectra != countJobs[["numSpectra"]]) {
      stop("Something went wrong in the counting. Wrong number of spectra for ",
           "msclusterId ", clust_members$clustId[[1]])
    }
    countJobs <- countJobs[!(names(countJobs) %in% c("numSpectra"))]
    countJobs[is.na(countJobs)] <- 0
    
    return(c(list(msclusterID = clust_members$clustId[[1]],
                  numSpectra = numSpectra,
                  mzConsensus = clust_members$mzConsensus[[1]],
                  rtMean = clust_members$rtMean[[1]],
                  rtMin = clust_members$rtMin[[1]],
                  rtMax = clust_members$rtMax[[1]],
                  peakIds = peakIds,
                  scans = scans,
                  joinedJobsIDs = joinedJobsIDs,
                  sumInts = clust_members$precursorInt[[1]]),
             countJobs))
  }))
  
  rm(clust_file)
  
  cat(" DONE!\n")
  return(clustersMembers)
})) %>% arrange(as.numeric(msclusterID))
# set the missing spectra counts to zero
clusterList[is.na(clusterList)] <- 0 

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
                               "rtMax", "peakIds", "scans", 
                               "joinedJobsIDs", "sumInts", "basePeakInt",
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_spectra")],
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_TOTAL")],
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_ion")],
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_area")])]

## check count sum
if (sum(clusterList$numSpectra) != 
    sum(clusterList[,names(clusterList)[endsWith(names(clusterList), "_spectra")]])) {
  stop("Spectra count sum check failed. Could not count number of spectra.")
}

# add flag columns based on the total columns present in the final joined counts
# compute indicators
blanksCode <- "BLANKS_TOTAL" %in% names(clusterList)
controlsCode <- "CONTROLS_TOTAL" %in% names(clusterList)
bedControlsCode <- "BEDS_TOTAL" %in% names(clusterList)
# DESREPLICATION AND HITS flag will be ignored here

if (blanksCode) {
  mz_blank <- almost.unique(clusterList[clusterList$BLANKS_TOTAL > 0, "mzConsensus"][[1]], mz_tol)
  mz_blank_up <- mz_blank + mz_tol # upper bound
  mz_blank <- mz_blank - mz_tol # lower bound
}
if (controlsCode) {
  mz_control <- almost.unique(clusterList[clusterList$CONTROLS_TOTAL > 0, "mzConsensus"][[1]], mz_tol)
  mz_control_up <- mz_control + mz_tol # upper bound
  mz_control <- mz_control - mz_tol # lower bound
}
if (bedControlsCode) {
  mz_bed_control <- almost.unique(clusterList[clusterList$BEDS_TOTAL > 0, "mzConsensus"][[1]], mz_tol)
  mz_bed_control_up <- mz_bed_control + mz_tol # upper bound
  mz_bed_control <- mz_bed_control - mz_tol # lower bound
}

# add BFLAG and CFLAG and and BEDFLAG
if (bedControlsCode || controlsCode || blanksCode)
{
  clusterList <- bind_cols(clusterList, 
                           bind_rows(lapply(seq_along(clusterList$mzConsensus), 
    function(i)
    {
    x <- clusterList$mzConsensus[[i]]
    
    flagColumns = list()
    
    if (blanksCode)
    flagColumns <- c(flagColumns, BFLAG = any(x >= mz_blank & x <= mz_blank_up))
    if (controlsCode)
    flagColumns <- c(flagColumns, CFLAG = any(x >= mz_control & x <= mz_control_up))
    if (bedControlsCode)
    flagColumns <- c(flagColumns, BEDFLAG = any(x >= mz_bed_control & x <= mz_bed_control_up))
    
    flagColumns
    })))
}

# TODO compute the peak are here correctly - using the original pre process data

# order columns to separate spectra and area count 
clusterArea <- clusterList[, c("msclusterID","numSpectra",
                               "mzConsensus", "rtMean",
                               "rtMin", "rtMax", "peakIds", "scans",
                               "joinedJobsIDs", "sumInts", "basePeakInt",
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_area")],
                               c("BFLAG", "CFLAG", "BEDFLAG")[
                                   c("BFLAG", "CFLAG", "BEDFLAG") %in%
                                     names(clusterList)],
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_ion")])]
# rename the area TOTAL columns removing the suffix equals _area
names(clusterArea)[endsWith(names(clusterArea), "_TOTAL_area")] <- gsub("_area$", "", names(clusterArea)[endsWith(names(clusterArea), "_TOTAL_area")])
# reorder columns to put the total at the end
clusterArea <- clusterArea[, c("msclusterID","numSpectra",
                               "mzConsensus", "rtMean",
                               "rtMin", "rtMax", "peakIds", "scans",
                               "joinedJobsIDs", "sumInts", "basePeakInt",
                               names(clusterArea)[endsWith(names(clusterArea),
                                                           "_area")],
                               c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                 "BFLAG", "CFLAG", "BEDFLAG")[
                                   c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                     "BFLAG", "CFLAG", "BEDFLAG") %in%
                                     names(clusterArea)],
                               names(clusterArea)[endsWith(names(clusterArea),
                                                           "_ion")])]
# filter the spectra count columns
clusterList <- clusterList[, c("msclusterID","numSpectra",
                               "mzConsensus", "rtMean",
                               "rtMin", "rtMax", "peakIds", "scans",
                               "joinedJobsIDs", "sumInts", "basePeakInt",
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_spectra")],
                               c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                 "BFLAG", "CFLAG", "BEDFLAG")[
                                   c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                     "BFLAG", "CFLAG", "BEDFLAG") %in%
                                     names(clusterList)],
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_ion")])]


dir.create(file.path(output_path, "count_tables"), showWarnings = FALSE)
write.csv(clusterList, 
          file = file.path(output_path, "count_tables", paste0(output_name, "_spectra.csv")), 
          row.names = FALSE)
# save count by area
write.csv(clusterArea,
          file = file.path(output_path, "count_tables", paste0(output_name, "_peak_area.csv")),
          row.names = FALSE)

# rm spec file from root folder, paste it in the clust folder
invisible(file.copy(from = file.path(output_path, paste0(output_name, "_0_spec_list.txt")), 
          to = file.path(output_path, 'clust'),
          recursive = TRUE, copy.mode = TRUE, copy.date = FALSE))
unlink(file.path(output_path, paste0(output_name, "_0_spec_list.txt")))