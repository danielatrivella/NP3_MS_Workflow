cat("Loading packages dplyr...\n")
suppressPackageStartupMessages(library(dplyr))

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Two arguments must be supplied to update the joined IDs in the tremolo result:\n", 
       " 1 - Path to the tremollo result;\n",
       " 2 - Path to a clean count CSV file.;\n", 
       call.=FALSE)
} else {
  path_tremolo_result <- file.path(args[[1]])
  if (!file.exists(path_tremolo_result))
  {
    stop("The tremolo result file '", path_tremolo_result,
         "' do not exists. Provide a valid output path to where csv file with the ",
         "tremolo result is located.")
  }
  
  path_clean_count <- file.path(args[[2]])
  if (!file.exists(path_clean_count))
  {
    stop("The count file '", path_clean_count,
         "' do not exists. Provide a valid output path to where csv file with the ",
         "clean count is located.")
  }
}

# read count files
ms_clean_count <- read.csv(path_clean_count, stringsAsFactors = FALSE,
                          comment.char = "", strip.white = TRUE)[, c("msclusterID", "joinedIDs")]
# get the joined ids
ms_clean_count <- ms_clean_count[!is.na(ms_clean_count$joinedIDs),]

if (nrow(ms_clean_count) == 0)
  q('no')

# read tremolo result
tremolo_result <- read.csv(path_tremolo_result, stringsAsFactors = FALSE,
                           comment.char = "", strip.white = TRUE)
# save original msclusterIDs
if (is.null(tremolo_result$scans))
  tremolo_result <- cbind(scans = tremolo_result$msclusterID, tremolo_result)

# ms_clean_count$joinedIDs <- lapply(ms_clean_count$joinedIDs, function(x) strsplit(x, ";")[[1]])
# overwrite joined scans ids by the parent id
invisible(apply(ms_clean_count, 1, function(x) 
  {
    joinedIDs <- as.numeric(strsplit(x[[2]], ";")[[1]])
    
    tremolo_result$msclusterID[tremolo_result$scans %in% joinedIDs] <<- as.numeric(x[[1]])
  }))

tremolo_result <- arrange(tremolo_result, msclusterID, desc(MQScore))
write.csv(tremolo_result, file = path_tremolo_result, row.names = FALSE)
