##
# Script to create the specification lists to be passed to the MSCluster when running 
# the joining jobs command to cluster the final clean mgfs
# Also join the metadata of joined jobs recursively and 
# concatenate the original jobs metadata containing the original samples - 
# check if they are all present
##

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Five arguments must be supplied to create the specification files for MSCLuster for joining jobs and concatenating the metadatas:\n", 
       " 1 - path_batch_metadata: Path to the CSV batch metadata file containing jobnames and job codes;\n",
       " 2 - path_jobs_data: Path to the jobs data folder;\n", 
       " 3 - path_pre_processed_dir: Path to the directory containing the original np3 jobs pre processed folder",
       " 4 - output_path: Path to the desired output folder location;\n", 
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
  # the original jobs from the joined metadata
  # otherwise, save the join jobs metadata directly
  if (any(metadata$JOINED_JOB == 1)) {
    # there is a joined job in the list of jobs to be joined,
    # extract the original joined jobs metadata (original_jobs_METADATA_JOIN.csv)
    # for each joined job in the list and concatenate this metadata with the 
    # remaining original samples in the list if any
    
    # read all the original joined jobs metadata
    if (sum(metadata$JOINED_JOB == 1) == 1)
    {
      # only one joined job to join
      i <- which(metadata$JOINED_JOB == 1)
      orig_joined_jobs_metadata <- readMetadataTableJoinJobs(file.path(metadata$JOB_PATH[i], "../..", 
                                                    "original_jobs_METADATA_JOIN.csv"),
                                          path_jobs_data, path_pre_processed_dir)
      orig_joined_jobs_metadata$JOINED_JOB_CODE <- metadata$JOB_CODE[i]
    } else {
      # more than one joined job to join
      orig_joined_jobs_metadata <- bind_rows(lapply(which(metadata$JOINED_JOB == 1), 
        function(i) {
          m_join <- readMetadataTableJoinJobs(file.path(metadata$JOB_PATH[i], "../..", 
                                                        "original_jobs_METADATA_JOIN.csv"),
                                              path_jobs_data, path_pre_processed_dir)
          m_join$JOINED_JOB_CODE <- metadata$JOB_CODE[i]
          m_join
        }))
    }
    
    # if there is any original job in the list to join, concatenate the joined jobs metadata
    if (any(metadata$JOINED_JOB == 0)) {
      orig_joined_jobs_metadata <- bind_rows(orig_joined_jobs_metadata,
                                             metadata[metadata$JOINED_JOB==0,])
    }
    metadata <- orig_joined_jobs_metadata
  } 
  
  # if there is no joined job in the list of jobs to be joined, 
  # save the join jobs metadata directly
  write.csv(metadata, 
            file = file.path(output_path, "original_jobs_METADATA_JOIN.csv"),
            row.names = FALSE)
  
  # Then concatenate the original samples metadata 
  # retrieve and concatenate the original jobs samples metadata
  orig_samples_metadata_path <- file.path(metadata$JOB_PATH, "../..", metadata$METADATA_NAME)
  if (!all(file.exists(orig_samples_metadata_path))) 
  {
    stop("Could not find the original samples metadata from the original np3 jobs for the joining job '", output_name, 
         "'. The samples' metadata table of the following jobs code do not exists or are misspelled in the join metadata:\n", 
         paste(metadata$JOB_CODE[!file.exists(orig_samples_metadata_path)], 
               collapse = "\n"),
         "\nPlease check the original samples' metadata name in the jobs output folder,",
         " match them with the METADATA_NAME in the provided join metadata and retry to run the np3 command.")
  }
  
  # read all the original samples metadata tables
  orig_samples_metadata <- bind_rows(sapply(seq_along(orig_samples_metadata_path), function(i) {
    m <- readMetadataTable(orig_samples_metadata_path[i])
    m$JOB_CODE <- metadata$JOB_CODE[i]
    m$JOINED_JOB_CODE <- metadata$JOINED_JOB_CODE[i]
    m
  }))
  
  # reorder the columns and save the original samples metadata for the joining jobs list
  if (!is.null(metadata$JOINED_JOB_CODE)) {
    orig_samples_metadata <- orig_samples_metadata[, c("FILENAME","SAMPLE_CODE","DATA_COLLECTION_BATCH","SAMPLE_TYPE","JOB_CODE","JOINED_JOB_CODE",
                                                       colnames(orig_samples_metadata)[!(colnames(orig_samples_metadata) %in% c("FILENAME","SAMPLE_CODE","DATA_COLLECTION_BATCH","SAMPLE_TYPE","JOB_CODE","JOINED_JOB_CODE"))])]
  } else {
    # do not include the joined_job_code columns
    orig_samples_metadata <- orig_samples_metadata[, c("FILENAME","SAMPLE_CODE","DATA_COLLECTION_BATCH","SAMPLE_TYPE","JOB_CODE",
                                                       colnames(orig_samples_metadata)[!(colnames(orig_samples_metadata) %in% c("FILENAME","SAMPLE_CODE","DATA_COLLECTION_BATCH","SAMPLE_TYPE","JOB_CODE"))])]
  }
  
  # check if the sample_codes are unique names among jobs - 
  # these will be necessary in the following workflow steps, for referencing the quantification columns
  if (any(duplicated(orig_samples_metadata$SAMPLE_CODE)))
  {
    stop("Wrong original samples metadata file format from the joining jobs. The following samples among jobs have duplicated sample codes:\n", 
         paste(orig_samples_metadata$SAMPLE_CODE[duplicated(orig_samples_metadata$SAMPLE_CODE)], 
               collapse = ", "),
         "\nThe sample code must be a *unique* syntactically valid name consisting of ", 
         "letters, numbers, and underscore characters, starting with a letter. ",
         "R reserved words are not syntactically valid names. When joining jobs, ",
         "the sample codes must also be *unique* among the different jobs.")
  }
  
  # save the original jobs samples metadata
  write.csv(orig_samples_metadata, 
            file = file.path(output_path, "original_samples_METADATA.csv"),
            row.names = FALSE)
}


# call function to create the batch lists and to join the original metadatas recursively
create_batch_lists_join_jobs_metadata(path_batch_metadata, path_jobs_data,
                                      path_pre_processed_dir,
                                      output_path, output_name)

