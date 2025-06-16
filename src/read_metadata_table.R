# from package bazar
almost.unique <- function(x, tol = sqrt(.Machine$double.eps))
{
  y <- round(x/tol, 0)
  d <- duplicated(y)
  x[!d]
}

#
# functions to read the metadata for running the np3 workflow major steps
readMetadataTable <- function(path_metadata, path_raw_data = NULL) 
{
  metadata <- read.csv(path_metadata, stringsAsFactors = FALSE, 
                             comment.char = "", strip.white = TRUE, colClasses = "character")
  
  return(checkMetadataFormat(metadata, path_raw_data))
}


checkMetadataFormat <- function(metadata, path_raw_data) {
  names(metadata) <- toupper(names(metadata))
  
  # check columns
  mandatory_columns <- c("FILENAME", "SAMPLE_CODE", "DATA_COLLECTION_BATCH", "SAMPLE_TYPE")
  
  if (!all(mandatory_columns %in% names(metadata)))
  {
    stop("Wrong batch metadata file format. It should have at least 4 columns named as follow:
  - FILENAME: must contain the raw MS/MS file name, without the path.
  - SAMPLE_CODE: must contain a unique syntactically valid sample code identifying each file
  - DATA_COLLECTION_BATCH: must contain a numeric index, enumerated sequentially starting with 1, indicading which files should be grouped in a same batch
  - SAMPLE_TYPE: must be 'blank' if the file is a blank, 'bed' or 'control' if the file is a control, 'hit' if the sample was a hit that need desreplication and 'sample' otherwise", call. = F)
  }
  
  # set sample type in lower case
  metadata$SAMPLE_TYPE <- tolower(metadata$SAMPLE_TYPE)
  # convert to numeric
  metadata$DATA_COLLECTION_BATCH <- as.numeric(metadata$DATA_COLLECTION_BATCH)
  # add raw data path to filename if present
  if (!is.null(path_raw_data)) {
    metadata$FILENAME <- file.path(path_raw_data, metadata$FILENAME)
  }
  
  # check if the sample codes are syntatic valid names
  if (!all(make.names(metadata$SAMPLE_CODE) == metadata$SAMPLE_CODE))
  {
    stop("Wrong metadata file format. The following samples do not have a syntactically valid sample code:\n", 
         paste(metadata$SAMPLE_CODE[which(make.names(metadata$SAMPLE_CODE) != metadata$SAMPLE_CODE)], collapse = ", "),
         "\nThe sample code must be a unique *syntactically valid* name consisting of ", 
         "letters, numbers, and underscore characters, starting with a letter. ",
         "R reserved words are not syntactically valid names. Spaces are ignored.")
  }
  
  # check if the sample codes are unique names
  if (any(duplicated(metadata$SAMPLE_CODE)))
  {
    stop("Wrong metadata file format. The following samples have duplicated sample codes:\n", 
         paste(metadata$SAMPLE_CODE[duplicated(metadata$SAMPLE_CODE)], 
               collapse = ", "),
         "\nThe sample code must be a *unique* syntactically valid name consisting of ", 
         "letters, numbers, and underscore characters, starting with a letter. ",
         "R reserved words are not syntactically valid names.")
  }
  
  # check if the data collection batch number are enumarated correctly
  if (!all(sort(unique(metadata$DATA_COLLECTION_BATCH)) == 
           1:length(unique(metadata$DATA_COLLECTION_BATCH))))
  {
    stop("Wrong metadata file format. The data collection batch number must be enumerated sequentially ", 
         "starting with 1. Refactor the following batch numbers sequence in the metadata file and retry: ", 
         paste(sort(unique(metadata$DATA_COLLECTION_BATCH)), collapse = ", "),
         ".")
  }
  
  # check if all samples types are valid
  if (!all(metadata$SAMPLE_TYPE %in% c("blank", "sample", "control", "bed", "hit")))
  {
    stop("Wrong metadata file format. The following sample types are not valid: ", 
         paste(metadata$SAMPLE_TYPE[which(!(metadata$SAMPLE_TYPE %in% 
                                                    c("blank", "sample", "control", "bed", "hit")))], collapse = ", "),
         "Accepted values are 'sample', 'hit', 'blank', 'control' or 'bed'.")
  }
  
  return(metadata)
}

#
# functions to read the metadata for joining jobs
readMetadataTableJoinJobs <- function(path_metadata, path_jobs_data = NULL, 
                                      path_pre_processed_dir = NULL) 
{
  metadata <- read.csv(path_metadata, stringsAsFactors = FALSE, 
                       comment.char = "", strip.white = TRUE, 
                       colClasses = "character")
  
  return(checkJoinJobsMetadataFormat(metadata, path_jobs_data, 
                                     path_pre_processed_dir))
}


checkJoinJobsMetadataFormat <- function(metadata, path_jobs_data, path_pre_processed_dir) {
  names(metadata) <- toupper(names(metadata))
  
  # check columns
  mandatory_columns <- c("JOBNAME", "JOB_CODE", "METADATA_NAME", 
                         "PRE_PROCESSED_DATA_NAME", "JOINED_JOB")
  
  if (!all(mandatory_columns %in% names(metadata)))
  {
    stop("Wrong join jobs metadata file format. It should have at least 2 columns named as follow:
  - JOBNAME: must contain the jobs final name, without the path;
  - JOB_CODE: must contain a unique syntactically valid job code identifying each job to be joined;
  - METADATA_NAME: must contain the name of the metadata table of the original np3 job if this is not a joined job, or the name of the metadata_join table of the joined job;
  - PRE_PROCESSED_DATA_NAME: must contain the nome of the pre process directory of this job, if this is not a joined job. Otherwise, empty;
  - JOINED_JOB: must contain 0 or 1 defining if this is a original np3 job (0) or a joined job result (1) from a join_jobs process, from different original np3 jobs.", 
         call. = F)
  }
  
  # create the jobs path concatenating the jobs data path to the jobname
  if (!is.null(path_jobs_data)) {
    metadata$JOB_PATH <- file.path(path_jobs_data, metadata$JOBNAME, 'outs', metadata$JOBNAME)
    
    if (!all(file.exists(metadata$JOB_PATH))) 
    {
      stop("Checking the join jobs metadata failed for '", output_name, 
           "'. The jobs paths of the following jobs code do not exists:\n", 
           paste(metadata$JOB_CODE[!file.exists(metadata$JOB_PATH)], 
                 collapse = ","),
           "\nPlease check if the jobs names and path are correctly defined in the join metadata and parameter and retry.")
    }
  }
  
  # check if the sample codes are syntatic valid names
  if (!all(make.names(metadata$JOB_CODE) == metadata$JOB_CODE))
  {
    stop("Wrong join jobs metadata file format. The following jobs do not have a syntactically valid job code:\n", 
         paste(metadata$JOB_CODE[which(make.names(metadata$JOB_CODE) != metadata$JOB_CODE)], collapse = ", "),
         "\nThe job code must be a unique *syntactically valid* name consisting of ", 
         "letters, numbers, and underscore characters, starting with a letter. ",
         "R reserved words are not syntactically valid names. Spaces are ignored.")
  }
  
  # check if the job codes are unique names
  if (any(duplicated(metadata$JOB_CODE)))
  {
    stop("Wrong join jobs metadata file format. The following jobs have duplicated job codes:\n", 
         paste(metadata$JOB_CODE[duplicated(metadata$JOB_CODE)], 
               collapse = ", "),
         "\nThe job code must be a *unique* syntactically valid name consisting of ", 
         "letters, numbers, and underscore characters, starting with a letter. ",
         "R reserved words are not syntactically valid names.")
  }
  
  # check if the JOINED_JOB column contains only 1 or 0 values
  if (!all(metadata$JOINED_JOB %in% c(0,1))) 
  {
    stop("Wrong join jobs metadata file format. The following jobs do not have a valid JOINED_JOB value:\n", 
         paste(metadata$JOB_CODE[!(metadata$JOINED_JOB %in% c(0,1))], collapse=", "),
         "\nThe JOINED_JOB must be a 0 or 1 value indicating if the respective ", 
         "job is an original np3 job or a np3 joined job from different original np3 jobs. ")
  }
  
  # if the directory pointing to the path where the pre processed results are located
  # is informed, also check if all the pre_processing are present - 
  # remove the joined jobs from the checking if any
  if (!is.null(path_pre_processed_dir) && any(metadata$JOINED_JOB == 0))
  {
    metadata$PRE_PROCESSED_DATA_PATH[metadata$JOINED_JOB == 0] <- file.path(path_pre_processed_dir, 
                                                                            metadata$PRE_PROCESSED_DATA_NAME[metadata$JOINED_JOB == 0])
    
    if (!all(file.exists(metadata$PRE_PROCESSED_DATA_PATH[metadata$JOINED_JOB == 0]))) 
    {
      stop("Checking the join jobs metadata failed for '", output_name, 
           "'. The pre processed data paths of the following jobs code do not exists:\n", 
           paste(metadata$JOB_CODE[metadata$JOINED_JOB == 0 & !file.exists(metadata$PRE_PROCESSED_DATA_PATH)], 
                 collapse = "\n"),
           "\nPlease check if the jobs names and path are correctly defined in the join metadata and parameter and retry.")
    }
  }
  
  return(metadata)
}