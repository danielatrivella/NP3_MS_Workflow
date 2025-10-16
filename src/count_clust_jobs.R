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
source(file.path(script_path(), "count_peak_area.R"))

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
  
  mz_tol <- as.numeric(args[[5]])
} 

# read the join jobs metadata
metadata <- readMetadataTableJoinJobs(path_jobs_metadata, path_jobs_data)

# To accept original sample codes from different original jobs with the same value, 
# here when reading the original clean data, the column names must be updated 
#to the new naming present in the original samples metadata, 
# column SAMPLES_CODE_ORIGINAL will store the original sample codes from the original jobs
# and the column SAMPLE_CODE_JOINED_JOBS will store the sample codes used in the last run, 
# this last code will be used to match the quantification columns names and updated them
# with new sample codes present in the SAMPLE_CODE column without duplicates

# read original samples metadata
metadata_samples <- readMetadataTable(file.path(output_path, "../..", 
                                                "original_samples_METADATA.csv"))
# get the samples that changed the code
diff_sample_code <- which(metadata_samples$SAMPLE_CODE != metadata_samples$SAMPLE_CODE_JOINED_JOBS)

# get the count of spectra of the subjobs
# and select only necessary columns that will be kept in the process
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
  # filter only the columns that will be used, ignore the rest - safer than removing the not wanted ones, 
  # this excludes any other column that the user may create
  countFile <- countFile[, c("msclusterID", "numSpectra", "peakIds", "scans", 
                             "joinedJobsIDs",  
                             names(countFile)[endsWith(names(countFile), "_spectra")])]
  
  # rename quantification columns using previous used name defined in 
  # the original samples metadata when there was a changed sample code
  # also replace new sample code naming in the original scans and peakIds references
  # only rename the changed samples codes (difference between SAMPLE_CODE and SAMPLE_CODE_JOINED_JOB)
  diff_sample_code_job <- diff_sample_code[diff_sample_code %in% which(metadata_samples$JOB_CODE == jobcode)]
  if (length(diff_sample_code_job) > 0)
  {
    # rename quantification columns of the new sample codes
    old_sample_code <- metadata_samples$SAMPLE_CODE_JOINED_JOBS[diff_sample_code_job]
    new_sample_code <- metadata_samples$SAMPLE_CODE[diff_sample_code_job]
    new_cols_names_idx <- match(paste0(old_sample_code, "_spectra"), names(countFile)) 
    names(countFile)[new_cols_names_idx] <- paste0(new_sample_code, "_spectra")
    # get rows where these samples appears and update its code in the 
    # respective peakIds and scans usgin regex
    for (i in seq_len(length(new_sample_code))) {
      rows_new_sample_code_i <- which(countFile[,new_cols_names_idx[i]]>0)
      pattern_peakIds_old_code <- paste0("(?<=;|^)([0-9]*|fake)(_",old_sample_code[i],")(?=;|$)")
      pattern_scans_old_code <- paste0("(?<=;|^)([0-9]*)(_",old_sample_code[i],")(?=;|$)")
      replace_new_code <- paste0("\\1_",new_sample_code[i])
      
      countFile$peakIds[rows_new_sample_code_i] <- sapply(countFile$peakIds[rows_new_sample_code_i],
        function(x) {
          gsub(pattern_peakIds_old_code, replace_new_code, x, perl = TRUE)
        }, USE.NAMES = FALSE)
      countFile$scans[rows_new_sample_code_i] <- sapply(countFile$scans[rows_new_sample_code_i],
        function(x) {
          gsub(pattern_scans_old_code, replace_new_code, x, perl = TRUE)
        }, USE.NAMES = FALSE)
    }
  }
  
  # not retrieving the peak area count, it will be correctly computed at the end
  # return the count of spectra
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
          count_subjobs[[subJob+1]]$msclusterID %in% ids, 5], collapse = ";")), # collapse joinedJobsIDs
        # sum the remaining columns - spectra counts
        colSums(count_subjobs[[subJob+1]][   
          count_subjobs[[subJob+1]]$msclusterID %in% ids, c(-3,-4,-5)], na.rm = TRUE)))  
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

# create the msclusterID_integrative columns here, it will be used in the clean step for 
# replacing the final msclusterIDs and here will serve as reference for original IDs from
# reference job; the msclusterID continues to be the result of MSCluster_NP3, this will
# be used for the pairwise similarity correct computation and in the cleaning

# integrative clustering heuristic
# maintain the first job ids as reference in a new column named msclusterID_integrative
# - in an integrative clustering manner - the rest of the clusters that do not
# appear in the first job should receive an incremental ID starting from the
# last ID of the first job that is present - first extract the IDs
# from the first job, them set the remaining IDs - 
# use the column joinedJobsIDs to retrieve the cluster IDs in the original 
# jobs using the first job in the metadata as the reference
reference_job_code <- metadata$JOB_CODE[[1]]
pattern_job_code_ref <- paste0("(?<=;|^)([0-9]*)(_",reference_job_code,")(?=;|$)")
max_original_ID <-  1
clusterList$msclusterID_integrative <- sapply(clusterList$joinedJobsIDs, function(x) {
  matched_job_code <- gregexpr(pattern_job_code_ref, x, perl = TRUE)
  # set the match length to the first captured group length using [,1]
  attr(matched_job_code[[1]], 'match.length') <- as.vector(attr(matched_job_code[[1]], 'capture.length')[,1])
  # extract first group matches (the original numerical IDs) - without the reference_job_code
  ref_ids <- as.integer(regmatches(x, matched_job_code)[[1]])
  # if any reference original id, return all of them concatenated with ; for reference - the minimum ID will be extracted in the clean step
  # keep the maximum ID from the reference that appears
  if (length(ref_ids) > 0) {
    max_original_ID <<- max(c(max_original_ID, ref_ids))
    paste(ref_ids, collapse = ";")
  } else {
    NA
  }
})
new_msclusterID_integrative <- ((max_original_ID+1):(max_original_ID+sum(is.na(clusterList$msclusterID_integrative))))
clusterList$msclusterID_integrative[is.na(clusterList$msclusterID_integrative)] <- new_msclusterID_integrative

if (any(duplicated(clusterList$msclusterID_integrative))) {
  stop("ERROR. An inconsistency was found in the msclusterIDs after applying the ",
       "integrative heuristic to maintain the original IDs of the first job ",
       reference_job_code," selected as reference, ",
       "and renumbering the m/zs that do not appear in this first reference job. ",
       "Something went wrong when renumbering the msclusterIDs. Please contact the dev team.")
}

# order columns
cluster_columns_order <- c("msclusterID", "numSpectra", 
                           "mzConsensus", "rtMean", "rtMin",
                           "rtMax", "peakIds", "scans", 
                           "joinedJobsIDs", "msclusterID_integrative", "sumInts")
clusterList <- clusterList[, c(cluster_columns_order,
                               names(clusterList)[endsWith(names(clusterList),
                                                           "_spectra")])]

# create the counts tables output folder
dir.create(file.path(output_path, "count_tables"), showWarnings = FALSE)
# save count by number of spectra - partial result - it will be overwritten at the end
# help in debugging
write.csv(clusterList, 
          file = file.path(output_path, "count_tables", paste0(output_name, "_spectra.csv")), 
          row.names = FALSE)

## check count sum
if (sum(clusterList$numSpectra) != 
    sum(clusterList[,names(clusterList)[endsWith(names(clusterList), "_spectra")]])) {
  stop("Spectra count sum check failed. Could not correctly count the total number of spectra. ",
       "Something went wrong in the quantification step. Please contact the dev team.")
}

# compute the peak area here correctly - using the original pre process data 
# from the original samples
# compute peak area count and the base peak int
clusterArea <- compute_peak_areas_joined_jobs(output_path, 
                                              clusterList$msclusterID,
                                              clusterList$scans,
                                              clusterList$peakIds, 
                                              clusterList$joinedJobsIDs)
clusterArea <- bind_cols(clusterList[, cluster_columns_order], clusterArea)
clusterList$basePeakInt <- clusterArea$basePeakInt

# use the original jobs samples metadata to recompute the sample types indicators, 
# for both spectra and area count

# compute indicators
blanksCode <- metadata_samples[metadata_samples$SAMPLE_TYPE == "blank", "SAMPLE_CODE"]
controlsCode <- metadata_samples[metadata_samples$SAMPLE_TYPE =="control", "SAMPLE_CODE"]
bedControlsCode <- metadata_samples[metadata_samples$SAMPLE_TYPE == "bed", "SAMPLE_CODE"]
hitsCode <- metadata_samples[metadata_samples$SAMPLE_TYPE == "hit", "SAMPLE_CODE"]

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
clusterList <- clusterList[, c(cluster_columns_order, "basePeakInt",
                               paste0(metadata_samples$SAMPLE_CODE,"_spectra"), 
                               c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                 "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG")[
                                   c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                     "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG") %in% 
                                     names(clusterList)])]
clusterArea <- clusterArea[, c(cluster_columns_order, "basePeakInt",
                               paste0(metadata_samples$SAMPLE_CODE,"_area"),
                               c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                 "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG")[
                                   c("BLANKS_TOTAL", "CONTROLS_TOTAL", "BEDS_TOTAL",
                                     "DESREPLICATION", "BFLAG", "CFLAG", "BEDFLAG", "HFLAG") %in%
                                     names(clusterArea)])]

# save the quantification table
# save count by number of spectra
write.csv(clusterList, 
          file = file.path(output_path, "count_tables", paste0(output_name, "_spectra.csv")), 
          row.names = FALSE)
# save count by peak area
write.csv(clusterArea,
          file = file.path(output_path, "count_tables", paste0(output_name, "_peak_area.csv")),
          row.names = FALSE)

# rm spec file from root folder, paste it in the clust folder
invisible(file.copy(from = file.path(output_path, paste0(output_name, "_0_spec_list.txt")), 
          to = file.path(output_path, 'clust'),
          recursive = TRUE, copy.mode = TRUE, copy.date = FALSE))
unlink(file.path(output_path, paste0(output_name, "_0_spec_list.txt")))