##
# step 9 - biocorrelation
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
source(file.path(script_path(), "read_metadata_table.R"))
suppressPackageStartupMessages(library(dplyr))

# default args
method <- "pearson"
bio_cutoff <- 0

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Two arguments must be supplied to compute the correlation between the samples number of spectra and bioactivity:\n",
       " 1 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches, blanks and correlation groups;\n",
       " 2 - Path to the CSV file containing the count of spectra or of peak area by cluster ID;\n",
       " 3 - Correlation method, one of: 'pearson', 'kendall', or 'spearman' (default);\n",
       " 4 - bioactivity cutoff - set to zero bioactivities values bellow this arg;\n",
       call.=FALSE)
} else {
  path_metadata <- file.path(args[[1]])
  if (!file.exists(path_metadata))
  {
    stop("The CSV batch metadata file '", path_metadata,
         "' do not exists. Provide a valid path to where the metadata file is located.")
  }
  path_metadata <- normalizePath(path_metadata)
  
  counts_path <- file.path(args[[2]])
  if (!file.exists(counts_path))
  {
    stop("The CSV count file '", counts_path,
         "' do not exists. Provide a valid path to where the count file is located.")
  }
  counts_path <- normalizePath(counts_path)
  
  if (length(args) > 2)
  {
    method <- args[[3]] 
    bio_cutoff <- as.numeric(args[[4]])
    if (is.na(bio_cutoff) || bio_cutoff < 0) {
      warning("The provided bioactivity cutoff have an invalid value or it is smaller than 0.",
              "It will be set to 0.",
              "If this behaviour is not wanted, rerun with a valid cutoff.")
    }
  }
}

output_path <- dirname(counts_path)
out_name <- sub("\\.[^\\.]*$", "", perl = T,  basename(counts_path))

metadata <- readMetadataTable(path_metadata)

# get the bioactivity columns
bioactivity <- names(metadata)[which(startsWith(names(metadata), "BIOACTIVITY_"))]

if (length(bioactivity) == 0)
{
  warning("No bioactivity score column with prefix 'BIOACTIVITY_' was present in the ", 
          basename(path_metadata), " metadata table. Correlation disabled.", call. = FALSE)
  q("no")
}

# check bioactivity columns
bioactivity <- unlist(sapply(bioactivity, function(bio_column) {
  metadata[, bio_column] <<- suppressWarnings(as.numeric(metadata[, bio_column]))
  
  if (any(is.na(metadata[, bio_column]))) 
  {
    warning("The bioactivity scores of the '", bio_column, 
            "' column contains not numeric values that will be set to zero. ",
            "If this behaviour is not wanted, refactor your metadata file and retry.",
            call. = F)
  }
  
  metadata[is.na(metadata[,bio_column]), bio_column] <<- 0 # set missing activities to zero
  
  if (any(metadata[,bio_column] < bio_cutoff)) 
  {
    warning("The bioactivity scores of the '", bio_column, 
            "' column contains values below the bioactivity cutoff (",
            bio_cutoff,
            ") and they will be set to zero. ",
            "If this behaviour is not wanted, refactor your metadata file and retry.",
            call. = F)
  }
  metadata[metadata[,bio_column] < bio_cutoff, bio_column] <<- 0 # apply bioactivity cut off
  
  if (all(metadata[,bio_column] == 0))
  {
    warning("The bioactivity scores of the '", bio_column, 
            "' column only contains zero values and will not be used. ",
            "If this behaviour is not wanted, refactor your metadata file and retry.",
            call. = F)
    NA
  } else {
    bio_column
  }
}, USE.NAMES = FALSE))
bioactivity <- bioactivity[!is.na(bioactivity)]

if (length(bioactivity) == 0)
{
  warning("No VALID bioactivity score column was present in the metadata CSV file. ",
          "To run the bioactivity correlation at least one of each one of the following ", 
          "columns must be present in the metadata file: \n",
          "  - BIOACTIVITY_<name>: (numeric values) must contain the samples bioactivity scores to be used in the correlation computation;
  - COR_<name>: (0 or 1's values) must contain at least one correlation group. \n\n", 
          "Where <name> is a unique identifier for each bioactivity score or correlation group that the user might want to add.", call. = F)
  q("no")
}

# get the list of correlations
corr_sets <- which(startsWith(names(metadata), "COR_"))

# check if exists any correlation group
if (length(corr_sets) == 0)
{
  warning("No correlation group column with prefix 'COR_' was present in the ", 
          out_name, " metadata CSV file. Correlation disabled.", call. = FALSE)
  q("no")
}

# check corr group columns, convert missing or NA values to 0
corr_sets <- unlist(sapply(corr_sets, function(corr_column) {
  # cat(corr_column,"\n")
  metadata[((metadata[,corr_column] == "") | 
              (is.na(metadata[,corr_column]))),corr_column] <<- 0
  if (!all(metadata[,corr_column] %in% c(0, 1)))
  {
    warning("The correlation group of the '", names(metadata)[corr_column], 
            "' column contains values different from 0 and 1, and thus will not be used. ",
            "If this behaviour is not wanted, refactor your metadata file and retry.",
            call. = F)
    NA
  } else {
    corr_column
  }
}, USE.NAMES = FALSE))
corr_sets <- corr_sets[!is.na(corr_sets)]

# check if exists any correlation group
if (length(corr_sets) == 0)
{
  warning("No VALID correlation group column was present in the metadata CSV file. ",
          "To run the bioactivity correlation at least one of each one of the following ", 
          "columns must be present in the metadata file: \n",
          "  - BIOACTIVITY_<name>: (numeric values) must contain the samples bioactivity scores to be used in the correlation computation;
  - COR_<name>: (0 or 1's values) must contain at least one correlation group. \n\n", 
          "Where <name> is a unique identifier for each bioactivity score or correlation group that the user might want to add.", call. = F)
  q("no")
}

count_samples <- read.csv(counts_path, stringsAsFactors = F, comment.char = "")

# add data type to the samples count column names
if (all(paste0(metadata$SAMPLE_CODE, "_area") %in% names(count_samples)))
{
  metadata$SAMPLE_CODE <- paste0(metadata$SAMPLE_CODE, "_area")
} else if (all(paste0(metadata$SAMPLE_CODE, "_spectra") %in% names(count_samples))) {
  metadata$SAMPLE_CODE <- paste0(metadata$SAMPLE_CODE, "_spectra")
} else {
  stop("Some samples of the metadata file are not present in the given count file.")
}
samplesNames <- metadata$SAMPLE_CODE

corr_bioactivity <- data.frame(metadata[, bioactivity], 
                               row.names = samplesNames)
colnames(corr_bioactivity) <- bioactivity

# compute the BCC for each corr_bioactivity columns
for (bio in bioactivity)
{
  bio_suffix <- sub("BIOACTIVITY", "", bio)
  # retrieve the columns within each correlations
  # remove correlations with less or equal than one sample selected with a valid bioactivity
  corr_samples_code <- lapply(corr_sets, function(i)
  {
    samples <- metadata[as.numeric(metadata[, i]) == 1 & 
                                !is.na(metadata[, i]), "SAMPLE_CODE"]
    
    if (length(samples) > 1) {
      # check if bioactivity is valid
      if (all(corr_bioactivity[samples, bio] == 0)) #
      {
        warning("Correlation group ", names(metadata)[[i]], 
                " only selected samples with the bioactivity score equals ", 
                "zero and thus will be skipped.", call. = F)
        NULL
      } else if (all(corr_bioactivity[samples, bio] == corr_bioactivity[samples, bio][[1]])) {
        warning("Correlation group ", names(metadata)[[i]], 
                " selected samples with the same bioactivity score ", 
                "(standard deviation equals zero) and thus will be skipped.", call. = F)
        NULL
      } else {
        samples
      }
    } else {
      warning("Correlation group ", names(metadata)[[i]], 
              " do not have more than one sample selected and thus will be skipped.", call. = F)
      NULL
    }
  })
  # remove invalid correlations
  invalid_corr <- sapply(corr_samples_code, is.null)
  corr_samples_code <- corr_samples_code[!invalid_corr]
  corr_sets <- corr_sets[!invalid_corr]
  
  # apply the biocorrelation
  count_samples <- bind_cols(count_samples, 
    lapply(seq_along(corr_samples_code), function(i)
    {
      # shaTest <- lapply(seq_len(nrow(count_samples)), function(x) {
      #   shapiro.test(as.numeric(count_samples[x, corr_samples_code[[i]]]))$p.value
      # })
      #which(shaTest > 0.05)
      if (!all(corr_samples_code[[i]] %in% names(count_samples)))
      {
        warning("The samples of the correlation group ", names(metadata)[corr_sets[[i]]], 
                " are not present in the count file '", counts_path, 
                "'. Skipping this group.", call. = F)
        return(NULL)
      }
      
      # suppressing warnings due to sd = 0 or div/0 when sum = 0
      corRes <- transmute(count_samples, 
                          !!paste0(names(metadata)[corr_sets][[i]], bio_suffix) := 
                            cor(t(select(count_samples,!!corr_samples_code[[i]])), 
                                                 corr_bioactivity[corr_samples_code[[i]], bio],
                                method = method)) 
      
      # if the bioactivity is not all equal, treat contant counts with sd = 0 substituting NAs
      if (sd(corr_bioactivity[corr_samples_code[[i]], bio]) != 0 && 
           sum(corr_bioactivity[corr_samples_code[[i]], bio]) != 0)
      {
        corNA <- is.na(corRes)
        countSumZero <- (rowSums(count_samples[, corr_samples_code[[i]]]) == 0)
        corRes[!countSumZero & corNA] <- "CTE"
      }
      
      corRes
    }))
}

# check if any correlation was performed
if (!any(startsWith(names(count_samples), "COR_")))
{
  stop("No correlation was computed for ", 
       out_name, 
       ". Probably a wrong count file was provided or a wrong selection of ", 
       "samples for the correlation groups was provided in the metadatada table.")
}

# get the max correlation score of each group
corrs_max <- sapply(which(startsWith(names(count_samples), "COR_")), 
                    function(j) {
                      max(count_samples[count_samples[,j]!= "CTE",j], na.rm = TRUE)
                    })

# samples order
samplesOrder <- match(names(count_samples)[which(names(count_samples) %in% samplesNames)],
                      samplesNames) 

# write the bioactivity and the max corrs value in the first row 
# skip the samplesNames[[1]] pos, the "BIOACTIVITY" tag and the bioactivity action (inhibition or activation) 
# to print the bioactivity value by sample
# skip the columns before the first corr to print the max corr by group
write.table(t(c(rep("", min(match(samplesNames, names(count_samples))) - 2), 
                bioactivity[1], corr_bioactivity[samplesOrder, bioactivity[1]], 
                rep("", which(startsWith(names(count_samples), "COR_"))[[1]]- max(match(samplesNames, names(count_samples))) -1),
                corrs_max)), 
            file = file.path(output_path, paste0(out_name, "_corr_", method, "_bioAct.csv")),
            sep = ",", col.names = FALSE, row.names = FALSE)
# write the bioactivities in the first rows
for (bio in bioactivity[-1])
{
  suppressWarnings(write.table(t(c(rep("", min(match(samplesNames, names(count_samples))) - 2), 
                  bio, corr_bioactivity[samplesOrder, bio])), 
              file = file.path(output_path, paste0(out_name, "_corr_", method, "_bioAct.csv")),
              sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)) # suppressing warning because of the param 'append'
}
# write the quantification table
suppressWarnings(write.table(count_samples, 
          file = file.path(output_path, paste0(out_name, "_corr_", method, "_bioAct.csv")), 
          row.names = FALSE, append = TRUE, sep = ",")) # suppressing warning because of the param 'append'
# write the quantification table without bioactivity data on the top
suppressWarnings(write.table(count_samples, 
                             file = file.path(output_path, paste0(out_name, "_corr_", method, ".csv")), 
                             row.names = FALSE, sep = ",")) # suppressing warning because of the param 'append'

# check if more than 10 warnings were emmited and show all unique ones
if (length(warnings()) > 10)
  unique(warnings())