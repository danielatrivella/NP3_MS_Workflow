##
# test the consistency of the count table - check the headers and MS2 completeness
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
cat("Loading packages readr...\n")
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
Rcpp::sourceCpp(file.path(script_path(),
                          '../src/read_mgf_peak_list_R.cpp'))
source(file.path(script_path(),"../src/read_metadata_table.R"))

min_num_peaks <- 5  # minimum number of fragmented peaks that a ms2 spectra must have to be kept in the final counts table

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Five arguments must be supplied to test the count tables consistency:\n", 
       " 1 - Path to the processed data folder (where the original MGFs are);\n",
       " 2 - Path to the spectra count table;\n", 
       " 3 - Path to the peak area count table;\n", 
       " 4 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks; For join_jobs this is the original samples metadata;\n",
       " 5 - Number of minimum peaks used to output a MS2 spectra - will be used to validate the MS2 completeness in the job",
       call.=FALSE)
} else {
  processed_data_path <- file.path(args[[1]])
  if (!dir.exists(processed_data_path))
  {
    stop("The processed data folder '", processed_data_path, 
         "' do not exists. Provide a valid path to where the original MGFs are.")
  }
  
  path_spectra_count <- file.path(args[[2]])
  if (!file.exists(path_spectra_count))
  {
    stop("The count file '", path_spectra_count,
         "' do not exists. Provide a valid path to where csv file with the count ",
         "of spectra is located.")
  }
  
  path_area_count <- file.path(args[[3]])
  if (!file.exists(path_area_count))
  {
    stop("The count file '", path_area_count,
         "' do not exists. Provide a valid path to where csv file with the count ",
         "of peak area is located.")
  }
  
  path_metadata <- file.path(args[[4]])
  if (!file.exists(path_metadata))
  {
    stop("The CSV batch metadata file '", path_metadata, 
         "' do not exists. Provide a valid path to where the metadata is located.")
  }
  
  min_num_peaks <- as.numeric(args[[5]])
}

ms_spectra_count <- suppressMessages(read_csv(path_spectra_count, guess_max = 10000))
ms_area_count <- suppressMessages(read_csv(path_area_count, guess_max = 10000))
metadata <- readMetadataTable(path_metadata)
n_inconsistency <- 0

convert_sample_processed_data_path <- function(is_join_jobs, processed_data_path) {
  if (!is_join_jobs) {
    convert_processed_data_path <- function(x) {
      return(processed_data_path)
    }
  } else {
    # join_jobs
    metadata_jobs <- readMetadataTableJoinJobs(file.path(dirname(path_metadata),
                                                         "original_jobs_METADATA_JOIN.csv"))
    convert_processed_data_path <- function(x) {
      return(metadata_jobs$PRE_PROCESSED_DATA_PATH[metadata_jobs$JOB_CODE == x][1])
    }
  }
  return(convert_processed_data_path)
}

retrieve_processed_data_path <- convert_sample_processed_data_path(is_join_jobs=("joinedJobsIDs" %in% names(ms_spectra_count)),
                                                                   processed_data_path = processed_data_path)

# check if the count tables are consistent
if (!isTRUE(all.equal(unlist(ms_spectra_count[,c("mzConsensus", "rtMean", "rtMin", "rtMax", "sumInts", "basePeakInt",
                                   "numSpectra", "peakIds", "scans", "peaksInt", "peaksList")]),
               unlist(ms_area_count[,c("mzConsensus", "rtMean", "rtMin", "rtMax", "sumInts", "basePeakInt",
                                "numSpectra", "peakIds", "scans", "peaksInt", "peaksList")]))))
{
  n_inconsistency <- n_inconsistency + 1
  warning("Some of the count tables columns - not the samples ",
          "quantification columns - are not consistent. ",
          "Some informations differs between the clean peak area and the clean spectra count tables.",
          "Differences:", all.equal(unlist(ms_spectra_count[,c("mzConsensus", "rtMean", "rtMin", "rtMax", "sumInts", "basePeakInt",
                                                               "numSpectra", "peakIds", "scans", "peaksInt", "peaksList")]),
                                    unlist(ms_area_count[,c("mzConsensus", "rtMean", "rtMin", "rtMax", "sumInts", "basePeakInt",
                                                            "numSpectra", "peakIds", "scans", "peaksInt", "peaksList")])))
}

# removed not used columns
ms_spectra_count <- ms_spectra_count[, !(startsWith(names(ms_spectra_count), "tremolo_") | 
                                         endsWith(names(ms_spectra_count), "_metfrag"))]
ms_spectra_count <- bind_cols(ms_spectra_count,  
                              ms_area_count[,endsWith(names(ms_area_count), "_area")])
rm(ms_area_count)

# get the count columns position
spectra_count_columns <- match(paste0(metadata$SAMPLE_CODE, "_spectra"), 
                               names(ms_spectra_count))

# check if the number of spectra by sample is consistent with the number of spectra
# in the scan list by sample and if the total #spectra is consistent with the 
# numSpectra column
scans_peakIds <- Reduce(rbind, (lapply(seq_len(nrow(ms_spectra_count)), function(i)
{
  scans <- strsplit(ms_spectra_count$scans[[i]], ";")[[1]]
  peakIds <- strsplit(ms_spectra_count$peakIds[[i]], ";")[[1]]
  
  num_spectra <- sum(sapply(seq_len(nrow(metadata)), function(j)
  {
    n_spectra_count <- ms_spectra_count[[i,spectra_count_columns[j]]]
    n_spectra_scans <- sum(grepl(pattern = paste0("$", metadata$SAMPLE_CODE[[j]], "$"), 
                                 fixed = TRUE, x = paste0(sub(pattern = "_", 
                                                              replacement = "$", scans), "$")))
    if (n_spectra_count != n_spectra_scans)
    {
      n_inconsistency <<- n_inconsistency + 1
      warning("The number of spectra of the sample ", metadata$SAMPLE_CODE[[j]],
              " for the msclusterID ", ms_spectra_count[[i,1]], 
              " does not match with the respectful number of scans.")
    }
    n_spectra_count
  }))
  
  if (num_spectra != ms_spectra_count$numSpectra[[i]]) {
    n_inconsistency <<- n_inconsistency + 1
    warning("Wrong number of spectra for the msclusterID ", ms_spectra_count[[i,1]], 
            ", it does not match with the respectful number of scans.")
  }
  
  list(scans = scans, peakIds = peakIds)
})))

min_basePeakInt <- min(ms_spectra_count$basePeakInt)
real_headers_total <- tibble(msclusterId = 0,
                             mz_w = 0,
                             rt_w = 0,
                             rtmin_w = 0,
                             rtmax_w = 0,
                             total_int = 0, .rows = nrow(ms_spectra_count))
# read each mgf from the samples and check if peakId, peak areas and 
# scan header is consistent for each msclusterId - i = which(ms_spectra_count$msclusterID == 193)
wrong_scans <- lapply(metadata$SAMPLE_CODE, function(x)
{
  preprocessed_mgf_path <- file.path(retrieve_processed_data_path(metadata$JOB_CODE[metadata$SAMPLE_CODE == x]), 
                                     paste0(x, '_peak_info.mgf'))
  if (!file.exists(preprocessed_mgf_path)) {
    stop("Error checking peakId, peak areas and scan header consistency of sample code ",x,
         ". Its pre processed MGF file does not exists: ", preprocessed_mgf_path, 
         ". Please provide a valid path to where the pre processed data is located.")
  }
  mgf_data <- readMgfHeader(preprocessed_mgf_path)
  mgf_data$scans <- paste0(mgf_data$scans, "_", x)
  
  # check consistency for each msclusterId
  msclusterId_wrong_area <- c()
  msclusterId_wrong_peakId <- c()
  real_headers <- bind_rows(lapply(seq_len(nrow(scans_peakIds)), function(i)
  {
    scans <- scans_peakIds[[i,1]]
    scans_x <- grepl(pattern = paste0("$", x, "$"), fixed = TRUE, 
                     x = paste0(sub(pattern = "_", replacement = "$", scans), "$")) # add $ at the end to prevent matching with a sample with the same prefix
    # at least one scan appears in this sample, check consistency
    if (any(scans_x)) {
      scans_header <- mgf_data[match(scans[scans_x], mgf_data$scans), ]
      peakIds <- scans_peakIds[[i,2]][grepl(pattern = paste0("$", x, "$"), fixed = T, 
                                            x = paste0(sub(pattern = "_", 
                                                           replacement = "$", 
                                                           scans_peakIds[[i,2]]), "$"))]
      
      # test if all peakIds are correct for this sample
      if (!all(peakIds %in% scans_header$peak_id) & 
          !all(scans_header$peak_id  %in% peakIds)) {
        n_inconsistency <<- n_inconsistency + 1
        warning("Wrong peakIds for msclusterId ", ms_spectra_count[[i,1]], 
                " in sample ", x, 
                ". The following peaks are not present in the original MGF: ",
                paste(peakIds[!(peakIds %in% scans_header$peak_id)],collapse = ","))
        msclusterId_wrong_peakId <<- c(msclusterId_wrong_peakId, ms_spectra_count[[i,1]])
      }
      
      # check if peak area of this sample is correct, with 1% of error due to rounding
      real_peak_area_x <- sum(unique(trunc(scans_header$peak_area)))
      if (ms_spectra_count[[i,paste0(x, "_area")]]/real_peak_area_x > 1.01 |
          real_peak_area_x/ms_spectra_count[[i,paste0(x, "_area")]] > 1.01) {
        n_inconsistency <<- n_inconsistency + 1
        warning("Wrong peak area for msclusterId ", ms_spectra_count[[i,1]], 
                " in sample ", x, 
                ". Correct area = ",real_peak_area_x, "; computed area = ", 
                ms_spectra_count[[i,paste0(x, "_area")]])
        msclusterId_wrong_area <<- c(msclusterId_wrong_area, ms_spectra_count[[i,1]])
      }
      list(msclusterId = ms_spectra_count[[i,1]],
           mz_w = sum(scans_header$mz*scans_header$int),
           rt_w = sum(scans_header$rt*scans_header$int),
           rtmin_w = sum(scans_header$rtmin*scans_header$int),
           rtmax_w = sum(scans_header$rtmax*scans_header$int),
           total_int = sum(scans_header$int))
    } else {
      list(msclusterId = ms_spectra_count[[i,1]],
           mz_w = 0,
           rt_w = 0,
           rtmin_w = 0,
           rtmax_w = 0,
           total_int = 0)
    }
  }))
  real_headers_total <<- real_headers_total + real_headers
  
  # check the completeness of the ms2 spectra
  # valid ms2 will receive 1 if the scans is present in the final count table
  # 2 if the scans could be removed by the filters (min num peaks and min base peak int)
  # 0 otherwise, meaning it is missing from the final count table - something went wrong
  mgf_data$valid <- 0 
  valid_scans <- sapply(mgf_data$scans, function(scan_num, x) {
    any(grepl(pattern = scan_num, fixed = TRUE, x = x))
  }, x = unlist(scans_peakIds[,1]))
  mgf_data$valid[valid_scans] <- 1
  mgf_data$valid[!valid_scans & 
                 (mgf_data$num_peaks < min_num_peaks | 
                  mgf_data$base_peak_int < min_basePeakInt)] <- 2
  if (any(mgf_data$valid == 0)) {
    n_inconsistency <<- n_inconsistency + sum(mgf_data$valid == 0)
    warning("Missing ", round(sum(mgf_data$valid == 0)/length(mgf_data$scans)*100,2) ,
            "% of valid MS2 spectra from sample ", 
            x, ". Raw scans that are missing: ", 
            paste(mgf_data$scans[mgf_data$valid==0], collapse = ","), "\n")
  }
  cat("Sample", x, "\t- # wrong peakIds ", length(msclusterId_wrong_peakId),
      "- # wrong peak areas ", length(msclusterId_wrong_area),
      "- # missing valid MS2 spectra ", sum(mgf_data$valid == 0),"\n")
  
  list(msclusterId_wrong_peakId = list(msclusterId_wrong_peakId),
       msclusterId_wrong_area = list(msclusterId_wrong_area),
       scans_missing = paste(mgf_data$scans[mgf_data$valid==0], collapse = ","))
})
warnings()

# only test consistency in the header of not merged counts
if (!grepl(pattern = "_merged_ann.csv", fixed = TRUE, x = path_area_count))
{
  # norm header fields
  real_headers_total$msclusterId <- real_headers_total$msclusterId/nrow(metadata)
  real_headers_total$mz_w <- real_headers_total$mz_w/real_headers_total$total_int
  real_headers_total$rt_w <- real_headers_total$rt_w/real_headers_total$total_int
  real_headers_total$rtmin_w <- real_headers_total$rtmin_w/real_headers_total$total_int
  real_headers_total$rtmax_w <- real_headers_total$rtmax_w/real_headers_total$total_int
  # set blanks rtmin to 0 and rtmax to 1000000
  if ('BLANKS_TOTAL' %in% names(ms_spectra_count)) {
    real_headers_total$rtmin_w[ms_spectra_count$BLANKS_TOTAL > 0] <- 0
    real_headers_total$rtmax_w[ms_spectra_count$BLANKS_TOTAL > 0] <- 1000000
  }
  
  diff_header <- real_headers_total[,-1]/ms_spectra_count[,c("mzConsensus", "rtMean", "rtMin", "rtMax", "sumInts")]
  if (any(diff_header > 1.05 | diff_header < 0.95, na.rm = TRUE))
  {
    warning("Some attributes of m/z and/or retention time are inconsistent")
    n_inconsistency <- n_inconsistency + sum(diff_header > 1.05 | 
                                             diff_header < 0.95, na.rm = TRUE)
  }
}

if (n_inconsistency == 0) {
  cat("Done! :)\n")
} else {
  cat("ERROR! A total of", n_inconsistency, "inconsistencies were found. :(\n")
}

#TODO Write test script for groups method.