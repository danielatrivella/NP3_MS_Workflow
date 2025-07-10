suppressPackageStartupMessages(library(dplyr))

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
Rcpp::sourceCpp(file.path(script_path(), 'read_mgf_peak_list_R.cpp'))
# compute_peak_area(processed_data_path, 
#                   ms_spectra_count$msclusterID,
#                   lapply(ms_spectra_count$scans, function(x) strsplit(x, ";")[[1]]),
#                   lapply(ms_spectra_count$peakIds, function(x) strsplit(x, ";")[[1]]),
#                   batch_metadata)

# retrieve the peak areas and max base peak int for each msclusterID based on
# the scans_count and check if the peakIds_count are consistent
# scans_count and peakIds_count: list of lists of scans splitted by ';'
# scans_count_format: list of lists of scans already formated, splitted by ';', with the '_' replaced 
# by $ and a $ added at the end of each scan - the pattern $sample_code$ will 
# be used to match with the sample codes of each job - the $ added at the end 
# prevents matching with a sample with the same prefix
compute_peak_area <- function(processed_data_path, msclusterIDs, scans_count,
                              peakIds_count, metadata, scans_count_format=NULL,
                              peakIds_count_format=NULL, job_code_match=NULL)
{
  # if the peakIds_count_format is not informed, format it from the peakIds_counts lists
  if (is.null(peakIds_count_format)) {
    peakIds_count_format <- lapply(peakIds_count, function(x) paste0(sub(pattern = "_", 
                                                                         replacement = "$", 
                                                                         x), "$"))
  }
  # the same for the scans_count
  if (is.null(scans_count_format)) {
    scans_count_format <- lapply(scans_count, function(x) paste0(sub(pattern = "_", 
                                                                         replacement = "$", 
                                                                         x), "$"))
  }
  # compute the peak area for each sample code
  peak_areas <- Reduce(bind_cols, lapply(metadata$SAMPLE_CODE, function(x)
  {
    preprocessed_mgf_path <- file.path(processed_data_path, paste0(x, '_peak_info.mgf'))
    if (!file.exists(preprocessed_mgf_path)) {
      stop("Error computing the peak area of sample code ",x,
           ". Its pre processed MGF file does not exists: ", preprocessed_mgf_path, 
           ". Please provide a valid path to where the pre processed data is located.")
    }
    mgf_data <- readMgfHeader(preprocessed_mgf_path)
    mgf_data$scans <- paste0(mgf_data$scans, "_", x)
    
    peak_area_x <- bind_rows(lapply(seq_len(length(scans_count)), function(i)
    {
      # when joining jobs, if the current id do not appear in this job, 
      # skip scans matching with the jobs' samples codes and return no matching directly
      if (!is.null(job_code_match) && !job_code_match[[i]]) {
        return(list(x = 0, base_peak_int = 0))
      }
      # match the scans within the current sample code, between $ symbols
      scans_f <- scans_count_format[[i]]
      scans <- scans_count[[i]]
      scans_x <- grepl(pattern = paste0("$", x, "$"), fixed = TRUE, 
                       x = scans_f)
      # at least one scan appears in this sample, check consistency 
      # and retrieve the unique peak area, summing with fake peaks, and the 
      # max base peak intensity for this sample
      if (any(scans_x)) {
        scans_header <- mgf_data[match(scans[scans_x], mgf_data$scans), ]
        peakIds <- peakIds_count[[i]][grepl(pattern = paste0("$", x, "$"), fixed = T, 
                                              x = peakIds_count_format[[i]])]
        
        # test if all peakIds are correct for this sample
        if (!all(peakIds %in% scans_header$peak_id) & 
            !all(scans_header$peak_id  %in% peakIds)) {
          stop("Wrong peakIds for msclusterId ", msclusterIDs[[i]], 
                  " in sample ", x, 
                  ". The following peaks are not present in the original MGF: ",
                  paste(peakIds[!(peakIds %in% scans_header$peak_id)],collapse = ","))
        }
        
        # sum the unique peak areas (fake peaks got integrated here because they have different peak_areas)
        real_peak_area_x <- sum(unique(trunc(scans_header$peak_area)))
        # return the maximum base peak found
        base_peak_int = max(scans_header$base_peak_int)
        
        list(x = real_peak_area_x, base_peak_int = base_peak_int)
      } else {
        list(x = 0, base_peak_int = 0)
      }
    }))
    names(peak_area_x) <- c(x, paste0('base_peak_int_',x))
    
    list(peak_area_x)
  }))
  if (length(metadata$SAMPLE_CODE) == 1) {
    peak_areas <- peak_areas[[1]]
  }
  # retrieve the maximum base peak int among all samples
  base_peak_int <- apply(peak_areas[,startsWith(names(peak_areas), "base_peak_int")], 1, max)
  # test if all base peak ints were retrieved correctly
  # when joining jobs, ignore the scans that do not appear in this job - they will have 0 base peak int
  if (is.null(job_code_match)) {
    if (any(base_peak_int == 0)) {
      stop("Wrong base peak intensities retrieved for msclusterId's ", 
           paste(msclusterIDs[(base_peak_int == 0)],collapse = ","))
    }
  } else {
    if (any(base_peak_int[job_code_match] == 0)) {
      stop("Wrong base peak intensities retrieved for msclusterId's ", 
           paste(msclusterIDs[(base_peak_int == 0) & job_code_match],
                 collapse = ","))
    }
  }
  # removes the base peak int from individual samples and set the final maximum
  # base peak int among samples
  peak_areas <- peak_areas[,!startsWith(names(peak_areas), "base_peak_int")]
  names(peak_areas) <- paste0(names(peak_areas), "_area")
  peak_areas$basePeakInt <- base_peak_int
  return(peak_areas)
}

# Compute peak area for different joined jobs
# call compute peak area for each original unique job_code and then bind the cols of the different jobs
# get the max base peak among jobs
# return areas concatenated by col
# How: uses the two default metadatas here: the join original samples to select the job codes
# and the original samples to filter the samples from this job codes
# use the output_path to retrieve the default metadatas path created in the setup step (create_batch_lists_join_jobs)
# output_path: path to the final output folder, inside the outs dir, named with the output_name
# msclusterIDs, scans, peakIds, joinedJobsIDs: the columns of the count tables - where
# the last three columns contain values separated by ;
compute_peak_areas_joined_jobs <- function(output_path, 
                                          msclusterIDs, scans, peakIds, 
                                          joinedJobsIDs)
{
  # split the concatenated scans, peakIds and joinedJobsIDs - and also
  # replace all _ by $, add a $ at the end of each value for matching
  scans <- lapply(scans, function(x) strsplit(x, ";")[[1]])
  scans_fomart <- lapply(scans, function(x) paste0(sub(pattern = "_", 
                                                       replacement = "$", 
                                                       x), "$"))
  peakIds <- lapply(peakIds, function(x) strsplit(x, ";")[[1]])
  peakIds_format <- lapply(peakIds, function(x) paste0(sub(pattern = "_", 
                                                                 replacement = "$",
                                                                 x), "$"))
  joinedJobsIDs <- lapply(joinedJobsIDs, function(x) {
    x <- strsplit(x, ";")[[1]]
    paste0(sub(pattern = "_", replacement = "$", x), "$")
  })
  
  # read the default metadata
  metadata_joined_jobs <- readMetadataTableJoinJobs(file.path(output_path, "../..", 
                                                              "original_jobs_METADATA_JOIN.csv"))
  metadata_samples <- readMetadataTable(file.path(output_path, "../..", 
                                                  "original_samples_METADATA.csv"))
  
  
  joined_jobs_peak_areas <- Reduce(bind_cols, lapply(metadata_joined_jobs$JOB_CODE, function(job_code) 
  {
    metadata_samples_job <- metadata_samples[metadata_samples$JOB_CODE == job_code,]
    pre_processed_data_path_job <- metadata_joined_jobs[metadata_joined_jobs$JOB_CODE == job_code,
                                                        "PRE_PROCESSED_DATA_PATH"]
    
    job_code_match <- grepl(pattern = paste0("$", job_code, "$"), fixed = TRUE, 
                     x = joinedJobsIDs)
    compute_peak_area(pre_processed_data_path_job,
                      msclusterIDs,
                      scans,
                      peakIds,
                      metadata_samples_job,
                      scans_fomart,
                      peakIds_format,
                      job_code_match)
  }))
  
  if (!(length(msclusterIDs) == NROW(joined_jobs_peak_areas))) {
    stop("Something went wrong when computing the peak areas of the joined jobs. ", 
         "The final number of m/zs in the peak area table of the joined jobs (", NROW(joined_jobs_peak_areas),
         ") is different than the number of provided msclusterIDs (", length(msclusterIDs),
         "). Please contact the dev team, an error probably occuried in the peak area counting process.")
  }
  
  # extract the maximum values from the basePeakInt columns from the different jobs
  max_basePeakInt <- apply(joined_jobs_peak_areas[,startsWith(names(joined_jobs_peak_areas), 
                                                              "basePeakInt")], 
                           1, max)
  # remove the basePeakInt columns from the different jobs and set a new one with the maximum among jobs
  joined_jobs_peak_areas <- joined_jobs_peak_areas[,!startsWith(names(joined_jobs_peak_areas), 
                                                                "basePeakInt")]
  joined_jobs_peak_areas$basePeakInt <- max_basePeakInt
  
  return(joined_jobs_peak_areas)
}
