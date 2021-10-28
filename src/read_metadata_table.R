# from package bazar
almost.unique <- function(x, tol = sqrt(.Machine$double.eps))
{
  y <- round(x/tol, 0)
  d <- duplicated(y)
  x[!d]
}


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