##
# Script to create the specification lists to be passed to the MSCluster when running 
# in batches of experiments with blank files indicated
##

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Five arguments must be supplied to create the specification files for MSCLuster:\n", 
       " 1 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks;\n",
       " 2 - Path to the raw data folder;\n", 
       " 3 - Path to the desired output folder location;\n", 
       " 4 - Output name of the data collection batches.",
       " 5 - The processed data dir or FALSE to not use this data.",call.=FALSE)
} else {
  path_batch_metadata <- file.path(args[[1]])
  # validate input
  if (!file.exists(path_batch_metadata))
  {
    stop("The CSV batch metadata file '", path_batch_metadata, 
         "' do not exists. Provide a valid path to where the metadata is located.")
  }
  path_batch_metadata <- normalizePath(path_batch_metadata)
  
  path_raw_data <- file.path(args[[2]])
  if (!dir.exists(path_raw_data))
  {
    stop("The raw data file folder '", path_raw_data, 
         "' do not exists. Provide a valid path to where the raw data is located.")
  }
  path_raw_data <- normalizePath(path_raw_data)
  
  output_path <- file.path(args[[3]])
  output_name <- args[[4]]
  use_processed <- args[[5]] # the processed data dir
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
spec_lists_path <- file.path(output_path, output_name, "spec_lists")
outs_path <- file.path(output_path, output_name, "outs")

metadata <- readMetadataTable(path_batch_metadata)

# check if the processed data exists, if not run with the raw data instead
if (dir.exists(file.path(path_raw_data, use_processed)))
{
  metadata$FILENAME <- file.path(path_raw_data, use_processed,
                                  paste0(metadata$SAMPLE_CODE, 
                                         "_peak_info.mgf"))
  # check if all files exists
  if (!all(file.exists(metadata$FILENAME)))
  {
    stop("Could not create the specification lists for '", output_name, 
         "'. The following pre processed data files do not exists:\n", 
         paste(metadata$FILENAME[!file.exists(metadata$FILENAME)], collapse = "\n"),
         "\nPlease run the pre_process command again.")
  }
} else { # no processed data
  stop("Could not create the specification lists for '", output_name, 
       "'. The pre processed folder '", file.path(path_raw_data, use_processed), 
       "' do not exists. Please run the pre_process command again.")
}

n_batches <- sort(unique(metadata$DATA_COLLECTION_BATCH))

# create output dir for storing spec lists and the MScluster output
dir.create(file.path(spec_lists_path), showWarnings = TRUE, recursive = TRUE)
if (!dir.exists(file.path(spec_lists_path)))
  stop("Could not create the directory to store the specification lists. Directory '",
       file.path(spec_lists_path), "' not found.")

blanks_total <- sum(metadata$SAMPLE_TYPE == "blank")

# for each batch create the specification list with the files to be run
# if just one batch with all files blanks or not blanks, just create the final integration step
if (length(n_batches) > 1 || (blanks_total > 0 && blanks_total < nrow(metadata))) 
{
  # if more than one batch
  if (length(n_batches) > 1)
  {
    dir.create(file.path(spec_lists_path, "batch_lists"),
               showWarnings = TRUE, recursive = TRUE)
  }
  
  # if there are blanks create dir for the blank lists
  if (blanks_total > 0) 
  {
    dir.create(file.path(spec_lists_path, "blank_lists"),
               showWarnings = TRUE, recursive = TRUE)
    
    # if there are also not blanks, create dir for the collections lists
    if (blanks_total < nrow(metadata))
    {
      dir.create(file.path(spec_lists_path, "data_lists"),
                 showWarnings = FALSE, recursive = TRUE)
    }
  }
  
  # create spec list for each batch
  invisible(lapply(n_batches, function(bat)
  {
    batch_data <- metadata[metadata$DATA_COLLECTION_BATCH == bat,]
    batch_blanks <- sum(batch_data$SAMPLE_TYPE == "blank")
    
    # there are blanks and more than one file and at least one data file in the data collection batch
    if (batch_blanks > 0 && nrow(batch_data) > 1 && batch_blanks < nrow(batch_data))
    {
      # create data lists
      if (length(batch_data[batch_data$SAMPLE_TYPE != "blank", "FILENAME"]) > 0)
      {
        write.table(batch_data[batch_data$SAMPLE_TYPE != "blank", "FILENAME"], 
                    file = file.path(spec_lists_path, "data_lists", paste0("B_", bat, "_0.txt")), 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
        subbatches <- c(0)
      } else {
        subbatches <- c()
      }
      
      # create blank lists
      write.table(batch_data[batch_data$SAMPLE_TYPE == "blank", "FILENAME"], 
                  file = file.path(spec_lists_path, "blank_lists", paste0("B_", bat, "_1.txt")), 
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      # invisible(lapply(seq_len(batch_blanks), 
      #                  function(n_blank, blank_name) 
      # {
      #   write.table(blank_name[[n_blank]], 
      #               file = file.path(spec_lists_path, "blank_lists", paste0("B_", bat, "_", n_blank, ".txt")), 
      #               row.names = FALSE, col.names = FALSE, quote = FALSE)
      # }, batch_data$FILENAME[batch_data$SAMPLE_TYPE == "blank"]))
      
      outs_dir <- paste0("B_", bat, "_", c(subbatches, 1))
      if (length(n_batches) > 1)
      {
      # create batch lists
        write.table(file.path(outs_path, outs_dir, "mgf", 
                              paste0(outs_dir, "_all.mgf")), 
                    file = file.path(spec_lists_path, "batch_lists", paste0("B_", bat, ".txt")), 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
      } else {
        # create the final output list directly
        write.table(file.path(outs_path, outs_dir, "mgf", 
                               paste0(outs_dir, "_all.mgf")), 
                    file = file.path(spec_lists_path, "out_spec_lists.txt"), 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
      }
    } else {
      # create batch list directly
      write.table(batch_data$FILENAME, 
                  file = file.path(spec_lists_path, "batch_lists", paste0("B_", bat, ".txt")), 
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }))
  
  if (length(n_batches) > 1)
  {
    # create the final output lists with all batches
    write.table(file.path(outs_path, paste0("B_", n_batches), "mgf", 
                          paste0(paste0("B_", n_batches), "_all.mgf")), 
                file = file.path(spec_lists_path, "out_spec_lists.txt"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
} else {
  # create the final output lists with all batches
  write.table(file.path(metadata$FILENAME), 
              file = file.path(spec_lists_path, "out_spec_lists.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
