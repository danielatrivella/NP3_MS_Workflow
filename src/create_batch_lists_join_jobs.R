##
# Script to create the specification lists to be passed to the MSCluster when running 
# the joining jobs command to cluster the final clean mgfs
##

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Four arguments must be supplied to create the specification files for MSCLuster for joining jobs:\n", 
       " 1 - Path to the CSV batch metadata file containing jobnames and job codes;\n",
       " 2 - Path to the jobs data folder;\n", 
       " 3 - Path to the desired output folder location;\n", 
       " 4 - Output name of the joined jobs.",call.=FALSE)
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
  
  output_path <- file.path(args[[3]])
  output_name <- args[[4]]
}

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
spec_lists_path <- file.path(output_path, "spec_lists")

metadata <- readMetadataTableJoinJobs(path_batch_metadata, path_jobs_data)

# check if the jobs final clean mgfs exists
if (!all(file.exists(file.path(metadata$JOBPATH, "mgf", 
                   paste0(metadata$JOBNAME,"_clean.mgf")))))
{
  stop("Could not create the specification lists for '", output_name, 
       "'. The following jobs mgf clean data do not exists:\n", 
       paste(metadata$JOBNAME[!file.exists(file.path(metadata$JOBPATH, "mgf", 
                                                     paste0(metadata$JOBNAME,"_clean.mgf")))], 
             collapse = "\n"),
       "\nPlease run the jobs again and retry.")
}

# create output dir for storing spec lists and the MScluster output
dir.create(file.path(spec_lists_path), showWarnings = TRUE, recursive = TRUE)
if (!dir.exists(file.path(spec_lists_path)))
  stop("Could not create the directory to store the specification lists. Directory '",
       file.path(spec_lists_path), "' not found.")

# create the specification list with final output clean mgf from the jobs to be joined
write.table(file.path(metadata$JOBPATH, "mgf", 
                      paste0(metadata$JOBNAME,"_clean.mgf")), 
            file = file.path(spec_lists_path, "out_spec_lists.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

