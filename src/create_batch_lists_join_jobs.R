##
# Script to create the specification lists to be passed to the MSCluster when running 
# the joining jobs command to cluster the final clean mgfs from different jobs
# Also join the metadata of joined jobs recursively and 
# concatenate the original jobs metadata containing the original samples - 
# check if they are all present - and fix duplicates among jobs
##

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Five arguments must be supplied to create the specification files for MSCLuster for joining jobs and concatenating the metadatas:\n", 
       " 1 - path_batch_metadata: Path to the CSV batch metadata file containing jobnames and job codes;\n",
       " 2 - path_jobs_data: Path to the jobs data folder;\n", 
       " 3 - path_pre_processed_dir: Path to the directory containing the original np3 jobs pre processed folder",
       " 4 - output_path: Path to the desired join jobs output folder location;\n", 
       " 5 - output_name: Output name of the joined jobs.",call.=FALSE)
} else {
  path_batch_metadata <- file.path(args[[1]])
  # validate input
  if (!file.exists(path_batch_metadata))
  {
    stop("The CSV join jobs metadata file '", path_batch_metadata, 
         "' do not exists. Provide a valid path to where the join jobs metadata table is located.")
  }
  path_batch_metadata <- normalizePath(path_batch_metadata)
  
  path_jobs_data <- file.path(args[[2]])
  if (!dir.exists(path_jobs_data))
  {
    stop("The jobs data file folder '", path_jobs_data, 
         "' do not exists. Provide a valid path to where the jobs data is located.")
  }
  path_jobs_data <- normalizePath(path_jobs_data)
  
  path_pre_processed_dir <- file.path(args[[3]])
  if (!dir.exists(path_pre_processed_dir))
  {
    stop("The jobs pre processed data folder '", path_pre_processed_dir, 
         "' do not exists. Provide a valid path to where the original np3 jobs pre processed data is located.")
  }
  path_pre_processed_dir <- normalizePath(path_pre_processed_dir)
  
  output_path <- file.path(args[[4]])
  output_name <- args[[5]]
}

library(purrr)
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

create_batch_lists_join_jobs_metadata <- function(path_batch_metadata,
                                                  path_jobs_data, path_pre_processed_dir,
                                                  output_path, output_name) 
{
  spec_lists_path <- file.path(output_path, "spec_lists")
  
  metadata <- readMetadataTableJoinJobs(path_batch_metadata, path_jobs_data, path_pre_processed_dir)
  
  # check if the jobs final clean mgfs exists
  if (!all(file.exists(file.path(metadata$JOB_PATH, "mgf", 
                     paste0(metadata$JOBNAME,"_clean.mgf")))))
  {
    stop("Could not create the specification lists for '", output_name, 
         "'. The following jobs mgf clean data do not exists:\n", 
         paste(metadata$JOBNAME[!file.exists(file.path(metadata$JOB_PATH, "mgf", 
                                                       paste0(metadata$JOBNAME,"_clean.mgf")))], 
               collapse = "\n"),
         "\nPlease check if the jobs names are correctly defined in the join metadata or run the jobs again and retry.")
  }
  
  if (nrow(metadata) < 2) {
    stop("At least two jobs must be included in the join jobs metadata to allow",
         " executing the join_jobs command. The current job named '", output_name, 
         "' has less than two jobs for joining and failed.")
  }
  
  # create output dir for storing spec lists and the MScluster output
  dir.create(file.path(spec_lists_path), showWarnings = TRUE, recursive = TRUE)
  if (!dir.exists(file.path(spec_lists_path))) {
    stop("Could not create the directory to store the specification lists. Directory '",
         file.path(spec_lists_path), "' not found.")
  }
  
  # create the specification list with final output clean mgf from the jobs to be joined
  write.table(file.path(metadata$JOB_PATH, "mgf", 
                        paste0(metadata$JOBNAME,"_clean.mgf")), 
              file = file.path(spec_lists_path, "out_spec_lists.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # create the joined metadatas for join jobs and original samples
  # if there is any joined job to be joined with other jobs, recursively extract 
  # the original jobs from the joined metadata_join and the original samples from the joined samples metadata
  # otherwise, save the join jobs metadata directly
  # # for new original jobs: saving original sample codes in a column called SAMPLE_CODE_ORIGINAL for reference when computing the peak area
  # and the last used sample codes in a column called SAMPLE_CODE_JOINED_JOBS for reference when updated the clean counts column names
  if (any(metadata$JOINED_JOB == 1)) {
    # there is a joined job in the list of jobs to be joined,
    # extract the original joined jobs metadata (original_jobs_METADATA_JOIN.csv) and the original samples metadata (original_samples_METADATA.csv)
    # for each joined job in the list and concatenate the joined jobs metadata with the 
    # remaining original samples in the list, if any,
    # and concatenate the original samples metadata with the metadata from original samples, if any
    
    # read all the original joined jobs metadata
    if (sum(metadata$JOINED_JOB == 1) == 1)
    {
      # only one joined job to join
      i <- which(metadata$JOINED_JOB == 1)
      orig_joined_jobs_metadata <- readMetadataTableJoinJobs(file.path(metadata$JOB_PATH[i], "../..", 
                                                    "original_jobs_METADATA_JOIN.csv"),
                                                    path_jobs_data=NULL, path_pre_processed_dir)
      orig_joined_jobs_metadata$JOINED_JOB_CODE <- metadata$JOB_CODE[i]
      
      # read the joined sample metadata, update current JOINED_JOB_CODE and the SAMPLE_CODE_JOINED_JOBS to the last used SAMPLE_CODE
      orig_samples_metadata <- readMetadataTable(file.path(metadata$JOB_PATH[i], "../..", 
                                                           "original_samples_METADATA.csv"))
      orig_samples_metadata$JOINED_JOB_CODE <- metadata$JOB_CODE[i]
      orig_samples_metadata$SAMPLE_CODE_JOINED_JOBS <- orig_samples_metadata$SAMPLE_CODE
      # if this is an join_jobs result before allowing repeated SAMPLE_CODE among different jobs, set the original sample code
      if (!("SAMPLE_CODE_ORIGINAL" %in% names(orig_samples_metadata))) {
        orig_samples_metadata$SAMPLE_CODE_ORIGINAL <- orig_samples_metadata$SAMPLE_CODE
      }
    } else {
      # more than one joined job to join
      orig_joined_jobs_metadata <- bind_rows(lapply(which(metadata$JOINED_JOB == 1), 
        function(i) {
          m_join <- readMetadataTableJoinJobs(file.path(metadata$JOB_PATH[i], "../..", 
                                                        "original_jobs_METADATA_JOIN.csv"),
                                              path_jobs_data=NULL, path_pre_processed_dir)
          m_join$JOINED_JOB_CODE <- metadata$JOB_CODE[i]
          m_join
        }))
      
      # read the joined sample metadata, update current JOINED_JOB_CODE and the SAMPLE_CODE_JOINED_JOBS to the last used SAMPLE_CODE
      orig_samples_metadata <- bind_rows(lapply(which(metadata$JOINED_JOB == 1), 
        function(i) {
          m_samples <- readMetadataTable(file.path(metadata$JOB_PATH[i], "../..", 
                                                        "original_samples_METADATA.csv"))
          m_samples$JOINED_JOB_CODE <- metadata$JOB_CODE[i]
          m_samples$SAMPLE_CODE_JOINED_JOBS <- m_samples$SAMPLE_CODE
          # if this is an join_jobs result before allowing repeated SAMPLE_CODE among different jobs, set the original sample code
          if (!("SAMPLE_CODE_ORIGINAL" %in% names(m_samples))) {
            m_samples$SAMPLE_CODE_ORIGINAL <- m_samples$SAMPLE_CODE
          }
          m_samples
        }))
    }
    
    # if there is any original job in the list to join, 
    # concatenate the joined jobs metadata
    # and the original samples metadata
    if (any(metadata$JOINED_JOB == 0)) {
      orig_joined_jobs_metadata <- bind_rows(orig_joined_jobs_metadata,
                                             metadata[metadata$JOINED_JOB==0,])
      # updated the SAMPLE_CODE_ORIGINAL and SAMPLE_CODE_JOINED_JOBS using the original SAMPLE_CODE
      orig_samples_metadata <- bind_rows(orig_samples_metadata,
                                         lapply(which(metadata$JOINED_JOB==0), 
                                                function(i) {
        m_samples_path <- file.path(metadata[i,"JOB_PATH"], "../..", 
                                    metadata[i,"METADATA_NAME"])
        m <- readMetadataTable(m_samples_path)
        m$JOB_CODE <- metadata$JOB_CODE[i]
        m$SAMPLE_CODE_ORIGINAL <- m$SAMPLE_CODE_JOINED_JOBS <- m$SAMPLE_CODE
        m
        }))
    }
    metadata_join <- metadata
    metadata <- orig_joined_jobs_metadata
  } else {
    # no joined job, only orignal jobs then concatenate their metadata
    orig_samples_metadata <- bind_rows(lapply(which(metadata$JOINED_JOB==0), 
            function(i) {
              m_samples_path <- file.path(metadata[i,"JOB_PATH"], "../..", 
                                          metadata[i,"METADATA_NAME"])
              m <- readMetadataTable(m_samples_path)
              m$JOB_CODE <- metadata$JOB_CODE[i]
              m$SAMPLE_CODE_ORIGINAL <- m$SAMPLE_CODE_JOINED_JOBS <- m$SAMPLE_CODE
              m
            }))
  }
  
  # if there is no joined job in the list of jobs to be joined, 
  # save the join jobs metadata directly
  write.csv(metadata, 
            file = file.path(output_path, "original_jobs_METADATA_JOIN.csv"),
            row.names = FALSE)
  
  # check for duplicated sample_codes among jobs, and resolve them 
  # duplicated codes are concatenate with _<i> where i is the number of duplicates
  # check if the sample_codes are unique names among jobs - and resolve duplicates
  # these will be necessary in the following workflow steps, for referencing the quantification columns
  if (any(duplicated(orig_samples_metadata$SAMPLE_CODE)))
  {
    duplicated_orig_sample_codes <- orig_samples_metadata$SAMPLE_CODE[duplicated(orig_samples_metadata$SAMPLE_CODE)]
    cat("\n  - WARNING: There are original samples with duplicated SAMPLE_CODE among different jobs:",
        paste(unique(duplicated_orig_sample_codes), collapse = ", "), 
        ". The duplicated values will receive a numerical suffix to differentiate",
        "them in the following steps.\n")
    
    # repeat while there is a duplicated SAMPLE_CODE
    while(any(duplicated(orig_samples_metadata$SAMPLE_CODE))) {
      duplicated_orig_sample_codes <- orig_samples_metadata$SAMPLE_CODE[duplicated(orig_samples_metadata$SAMPLE_CODE)]
      # remove remaining duplicated sample code
      for (duplicated_sample_code in unique(duplicated_orig_sample_codes)) {
        # create the new sample codes
        new_sample_codes <- sapply(seq_len(sum(duplicated_orig_sample_codes == duplicated_sample_code)), 
               function(i, x=duplicated_sample_code) { paste0(x,"_",i)})
        # if any of them also exists, increment the last digits of all until there is no repetition
        # these way the new codes keep an ascending suffix value
        while(any(new_sample_codes %in% orig_samples_metadata$SAMPLE_CODE)) {
          new_sample_codes <- sapply(new_sample_codes, function(x) {
              y <- strsplit(x,"_")[[1]]
              n <- length(y)
              y[n] <- as.integer(y[n])+1
              return(paste0(y, collapse="_"))
            }, USE.NAMES = FALSE)
        }
        # reset current sample code
        orig_samples_metadata$SAMPLE_CODE[
          which(orig_samples_metadata$SAMPLE_CODE == duplicated_sample_code)[-1]] <- new_sample_codes
      }
    } 
  }
  
  # reorder the columns and save the original samples metadata for the joining jobs list
  metadata_cols_order <- c("FILENAME","SAMPLE_CODE", "SAMPLE_CODE_JOINED_JOBS", 
                           "SAMPLE_CODE_ORIGINAL", "DATA_COLLECTION_BATCH",
                           "SAMPLE_TYPE","JOB_CODE")
  if (!is.null(metadata$JOINED_JOB_CODE)) {
    metadata_cols_order <- c(metadata_cols_order, "JOINED_JOB_CODE")
    orig_samples_metadata <- orig_samples_metadata[, c(metadata_cols_order,
                                                       colnames(orig_samples_metadata)[
                                                         !(colnames(orig_samples_metadata) %in% metadata_cols_order)])]
    # if there was a joined job, set the metadata as the metadata_join for the last check
    # based on the joining jobs clean result
    metadata <- metadata_join
  } else {
    # do not include the joined_job_code columns
    orig_samples_metadata <- orig_samples_metadata[, c(metadata_cols_order,
                                                       colnames(orig_samples_metadata)[
                                                         !(colnames(orig_samples_metadata) %in% metadata_cols_order)])]
  }
  
  # save the original jobs samples metadata
  write.csv(orig_samples_metadata, 
            file = file.path(output_path, "original_samples_METADATA.csv"),
            row.names = FALSE)
  
  # check if the columns with the samples quantifications from the original jobs clean count tables match the
  # list of samples present in the orig_samples_metadata - samples consistency test
  count_samples_list <- unlist(lapply(metadata$JOBNAME, function(x)
  {
    jobcode <- metadata[metadata$JOBNAME == x,"JOB_CODE"]
    # check if the clean spectra table exists
    spectraCountFilePath <- file.path(path_jobs_data, x, "outs", x, "count_tables", 
                                      "clean", paste0(x, "_spectra_clean_ann.csv"))
    if (!file.exists(spectraCountFilePath)) 
    {
      stop("Could not find the original clean count table from JOB_CODE = ", 
           jobcode,
           ". Check if this job output is consistent and no original file is missing or was modified. ",
           "Path to the missing count table file: ", spectraCountFilePath)
    }
    # retrieve the spectra quantification columns
    countFileHeader <- read.csv(spectraCountFilePath, 
                          stringsAsFactors = F, comment.char = "", nrows = 1)
    # get the sample codes present in the quantification columns
    jobSamples <- sub("_spectra$", "",  names(countFileHeader)[endsWith(names(countFileHeader), "_spectra")])
    # check if all samples are present in the original samples metadata
    if (!all(jobSamples %in% orig_samples_metadata$SAMPLE_CODE_JOINED_JOBS)) {
      stop("There are original samples present in the original clean counts of JOB_CODE = ", jobcode,
           " that are missing from the provided original samples metadata. ",
           " The following SAMPLE_CODEs are missing: ", 
           paste(jobSamples[!(jobSamples %in% orig_samples_metadata$SAMPLE_CODE_JOINED_JOBS)], collapse= ","), 
           " Please make sure all original data is not modified and fully present in the provided paths ",
           "to prevent inconsistency errors during the join_jobs integrative clustering")
    }
    return(jobSamples)
  }))
  if (!all(orig_samples_metadata$SAMPLE_CODE_JOINED_JOBS %in% count_samples_list)) {
    stop("There are original samples present in the original samples metadata that are missing from ", 
         " the original samples clean counts. The following SAMPLE_CODEs are missing: ", 
         paste(orig_samples_metadata$SAMPLE_CODE_JOINED_JOBS[!(orig_samples_metadata$SAMPLE_CODE_JOINED_JOBS %in% count_samples_list)], collapse= ","), 
         " Please make sure all original data is not modified and fully present in the provided paths ",
         "to prevent inconsistency errors during the join_jobs integrative clustering")
  }
}


# call function to create the batch lists and to join the original metadatas recursively
create_batch_lists_join_jobs_metadata(path_batch_metadata, path_jobs_data,
                                      path_pre_processed_dir,
                                      output_path, output_name)

