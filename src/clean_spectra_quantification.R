## ----load-libs, message = FALSE--------------------------------------------
cat("Loading packages Rcpp, readr, dplyr, CPP functions...\n")
suppressPackageStartupMessages(library(readr))
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

Rcpp::sourceCpp(file.path(script_path(), 'triangular_matrix_R.cpp'))
Rcpp::sourceCpp(file.path(script_path(), 'read_mgf_peak_list_R.cpp'))
Rcpp::sourceCpp(file.path(script_path(), 'norm_dot_product.cpp'))
source(file.path(script_path(), "count_peak_area.R"))
source(file.path(script_path(), "read_metadata_table.R"))
source(file.path(script_path(), "writeMgfData_NP3.R"))

options(digits=10) # increase precision
options(readr.show_progress = FALSE)

# default args values
mz_tol <- 0.05
sim_tol <- 0.55 # sim to join
rt_tol <- 2
bin_size <- 0.025
blank_depth <- 3
scale_factor <- 0.5
table_limit_size <- 3000  # set max number of rows to process in a chunck
mz_rt_digits <- 4  # number of digits to round the mzConsensus and the retention times
ion_mode <- "+" # + or -
# A positive numeric value to scale the interquartile range (IQR) of the 
# blank spectra basePeakInt distribution and allow spectra with a basePeakInt
# value below this distribution median plus IQR*bflag_cutoff to be joined with 
# a blank spectrum without relying on the similarity value.
# Or FALSE to disable it.
# The IQR is the range between the 1st quartile (25th quantile) and the 3rd 
# quartile (75th quantile). The spectra with a basePeakInt value <= median + 
# IQR*bflag_cutoff of the blank spectra basePeakInt distribution and BFLAG TRUE 
# will be joined to a blank spectrum in the clean step. This cutoff will affect 
# the spectra with BFLAG TRUE that would not get joined to a blank spectra when 
# relying only on the similarity cutoff. This is a turn around to the fact that 
# blank spectra have low quality spectra and thus can not fully rely on the similarity values.
BFLAG_cutoff_factor <- 1.5 
NOISE_cutoff_factor <- FALSE

RMSE <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  sqrt(mean((x - y)^2))
}

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Two arguments must be supplied to clean and annotate the counts:\n", 
       " 1 - Path to the output data folder, inside the outs directory of the clustering result folder. ", 
       "It should contain the mgf folder, the peak area count CSV and the spectra count CSV. The data name will be extracted from here.;\n", 
       " 2 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks; if this is a join_jobs result the metadata must be set as '-1';\n",
       " 3 - Path to the pre processed data folder were the MGFs were created. Used to compute the peak areas;\n",
       " 4 - The precursor m/z tolerance in Da;\n",
       " 5 - The similarity tolerance to JOIN mass dissipation spectra;\n",
       " 6 - The retention time tolerance to enlarge the peaks before comparisions;\n",
       " 7 - The bin size in Da to consider two fragmented peaks m/z's the same;\n",
       " 8 - The scale factor to be used to aggregated the peaks of joined clusters;\n",
       " 9 - The ionization mode, one of -1 or 1 (default to 1);\n",
       " 10 - The bflag_cutoff a positive numeric or FALSE to disable it. If positive numeric, allow to join the spectra with a low basePeakInt (<= median + bflag_cutoff*IQR) with a blank spectra without relying on the similarity value;\n",
       " 10 - The NOISE_cutoff a positive numeric or FALSE to disable it. If positive numeric, allow to remove the spectra with a low basePeakInt  (<= median + NOISE_cutoff*IQR);\n",
       " 11 - The maximum number of spectra (rows) to be processed at a time.",
       call.=FALSE)
} else {
  output_path <- file.path(args[[1]])
  if (!dir.exists(output_path))
  {
    stop("The job output folder '", output_path, 
         "' do not exists. Provide a valid path to where the job final result is located.")
  }
  output_name <- basename(output_path)
  
  path_area_count <- file.path(output_path, "count_tables", paste0(output_name,"_peak_area.csv"))
  if (!file.exists(path_area_count))
  {
    stop("The count file '", path_area_count,
         "' do not exists. Provide a valid output path to where csv file with the count ",
         "of peak area is located.")
  }
  
  path_spectra_count <- file.path(output_path, "count_tables", paste0(output_name,"_spectra.csv"))
  if (!file.exists(path_spectra_count))
  {
    stop("The count file '", path_spectra_count,
         "' do not exists. Provide a valid output path to where csv file with the count ",
         "of spectra is located.")
  }
  
  path_batch_metadata <- args[[2]]
  if (path_batch_metadata == "-1") {
    # this is a joinig job process, the metadata should be extracted from the 
    # default path to the original samples metadata - located in the output path
    joining_jobs <- TRUE
    path_batch_metadata <- file.path(output_path, "../..", 
                                            "original_samples_METADATA.csv")
    if (!file.exists(path_batch_metadata))
    {
      stop("The metadata with the original samples from the joined jobs '", path_batch_metadata,
           "' do not exists. Provide a valid output path to where the csv file with  ",
           "the default samples metadata is located.")
    }
  } else {
    joining_jobs <- FALSE
    path_batch_metadata <- file.path(args[[2]])
    if (!file.exists(path_batch_metadata))
    {
      stop("The CSV batch metadata file '", path_batch_metadata, 
           "' do not exists. Provide a valid path to where the metadata is located.")
    }
  }
  
  if (!joining_jobs) {
    # ignore the pre processed data dir when joining jobs, 
    # this info will be extracted from the jobs metadata
    processed_data_path <- file.path(args[[3]])
    if (!dir.exists(processed_data_path))
    {
      stop("The processed data folder '", processed_data_path, 
           "' do not exists. Provide a valid path to where the pre processed MGFs are located.")
    }
  }
  
  if (length(args) > 11) {
    mz_tol <- as.numeric(args[[4]])
    if (is.na(mz_tol))
      stop("The m/z tolerance must be a numerica value. Wrong value informed: ",mz_tol)
    sim_tol <- as.numeric(args[[5]]) # sim to join
    if (is.na(sim_tol))
      stop("The similarity tolerance must be a numerica value. Wrong value informed: ",sim_tol)
    rt_tol <- as.numeric(args[[6]])
    if (is.na(rt_tol))
      stop("The retention time tolerance must be a numerica value. Wrong value informed: ",rt_tol)
    bin_size <- as.numeric(args[[7]])
    if (is.na(bin_size))
      stop("The bin size must be a numerica value. Wrong value informed: ",bin_size)
    
    scale_factor <- as.numeric(args[[8]])
    if (!is.na(scale_factor) || !is.null(scale_factor))
    {
      scale_factor <- as.numeric(scale_factor[[1]])
      if (scale_factor < 0)
      {
        warning("Invalid scale factor provided, it must be a numeric value greater or equal to 0: \n", 
                "  - 0 : ln natural logarithm of the intensities\n",
                "  - 1 : no scale\n",
                "  - x : x > 0 pow of the intensities to x. (e.g. x = 0.5 square root)\n",
                "Scale factor equals 0.5 (sqrt) will be selected by default.", call. = FALSE)
        scale_factor <- 0.5
      }
    } else {
      scale_factor <- 1 # no scale
    }
    
    ion_mode <- as.numeric(args[[9]])
    if (!(ion_mode %in% c(-1,1)))
      stop("The ion_mode arg must be a numeric value indicating the precursor ", 
           "ion mode. One of the following valid numeric values corresponding ", 
           "to a ion adduct type: '1' = positive, or '-1' = negative",  call. = FALSE)
    ion_mode <- ifelse(ion_mode > 0, "+", "-")
    
    BFLAG_cutoff_factor <- args[[10]]  # check if it is FALSE
    if (BFLAG_cutoff_factor != "FALSE") { 
      # check if it is numeric
      BFLAG_cutoff_factor <- as.numeric(args[[10]]) 
      if (is.na(BFLAG_cutoff_factor) || BFLAG_cutoff_factor < 0)
        stop("The BFLAG_cutoff arg must be a positive numeric value or FALSE to disable it. ",
             "If positive numeric value, indicates that the features with low basePeakInt and BFLAG TRUE should be joined with ",
             "a blank mz without relying on the similarity value.",  call. = FALSE)
    } else {
      BFLAG_cutoff_factor <- as.logical(args[[10]]) 
    }
    
    NOISE_cutoff_factor <- args[[11]]  # check if it is FALSE
    if (NOISE_cutoff_factor != "FALSE") { 
      # check if it is numeric
      NOISE_cutoff_factor <- as.numeric(args[[11]]) 
      if (is.na(NOISE_cutoff_factor) || NOISE_cutoff_factor < 0)
        stop("The NOISE_flag arg must be a positive numeric value or FALSE to disable it. ",
             "If positive numeric value, indicates that the features with low basePeakInt should be removed.",  
             call. = FALSE)
    } else {
      NOISE_cutoff_factor <- as.logical(args[[11]]) 
    }
    
    table_limit_size <- as.numeric(args[[12]])
    if (table_limit_size < 100) {
      warning("The max number of spectra (rows) to be processed at a time was too small. ",
              "Setting it to 100.", call. = FALSE)
      table_limit_size <- 100
    }
  } 
}
#print(args)

inverse_scale <- function(ints, scale_factor)
{
  if (scale_factor == 0) # log scale was applied
  {
    # inverse of 1.0 + log(1.0 + multVal * ints[i]))
    ints <- expm1(ints-1)
  } else {
    ints <- ints^(1/scale_factor)
  }
  return(round(ints,5))
}

scaleInts <- function(ints, scale_factor)
{
  if (scale_factor == 0) # log scale was applied
  {
    # following 1.0 + log(1.0 + multVal * ints[i]))
    ints <- 1.0 + log(1.0+ints)
  } else {
    ints <- ints ^ scale_factor
  }
  return(round(ints,5))
}

# aggregate the rows of the similarity table for the joined spectra ids
# use the maximum among values to aggregate the similarities
aggregate_sim_table <- function(joined_ids, scans_order, rm_rows, col_types, 
                                nscans, sim_file, table_limit_size) 
{
  # print(joined_ids)
  n_scans <- length(scans_order)
  
  # obtain number of chunks to be used to read the sim table
  n_chunks <- ceiling(n_scans/table_limit_size)
  
  last_row_joined <- sapply(joined_ids, max)
  keep_rows <- which(!rm_rows)
  x_rows <- as.numeric(names(joined_ids))
  # add header
  write_csv(as.data.frame(matrix(scans_order[keep_rows], nrow = 1)), 
            path = file.path(output_path, 
                             "molecular_networking/similarity_tables", 
                             paste0("similarity_table_", output_name, "_tmp.csv")),
            col_names = FALSE)
  
  # for each join, aggregate the similarities
  to_row <- 1
  cat(paste0("            |", paste0(rep(" ", n_chunks), collapse = ""), "|\n            |", collapse = ""))
  for (i in seq_len(n_chunks)) 
  {
    cat("=")
    from_row <- to_row
    to_row <- min(max(round(n_scans/n_chunks*(i)), to_row + round(n_scans/n_chunks)), n_scans)
    if (n_scans - to_row == 1) # prevent resting one line
      to_row <- n_scans
    
    # check if the ids that are inside this chunk aggregate rows that are also inside,
    # if not enlarge chunk until all x_rows ids are included
    repeat {
      x_rows_chunk <- x_rows >= from_row & x_rows <= to_row
      to_row <- max(last_row_joined[x_rows_chunk][last_row_joined[x_rows_chunk] >= to_row], to_row)
      # check if new x_rows were included after enlarging the chunk, repeat until no new inclusion is possible
      if (!any(xor(x_rows_chunk, (x_rows >= from_row & x_rows <= to_row))))
        break
    }
    n_joins <- sum(rm_rows[1:from_row]) # number of joins until this row is reached - row 1 is the header
    n_rows <- to_row - from_row
    sim_chunk <- as.matrix(read_csv(sim_file, col_names = FALSE, skip = from_row, 
                                    n_max = n_rows, col_types = col_types))
    sim_chunk <- sim_chunk[,-1] # rm scans number col
    
    if (!is.double(sim_chunk))
      sim_chunk <- as.double(sim_chunk)
    
    # mirror similarity table simetric for the case where there are rows between the ids that
    # will be joined, garanting completeness
    # mirror from column from_row to to_row and all rows
    if (!is.null(nrow(sim_chunk)))
      mirror_matrix_tri_upper(sim_chunk, from_row, to_row-1)
    
    # aggregate columns
    # only if the selected cols for aggregation are greater or equal than from_row
    # else, only rm the aggregated ones - their values are NA here, no need of aggregating
    select_cols_ge_from_row <- (x_rows-1 >= from_row)
    if (any(x_rows-1 >= from_row)) {
      sim_chunk[,(x_rows-1)[select_cols_ge_from_row]] <- matrix(apply(sim_chunk, 1, function(x) {
        sapply(joined_ids[select_cols_ge_from_row], function(j) max(x[j-1], na.rm = T))}),
        ncol = length((x_rows-1)[select_cols_ge_from_row]), byrow = TRUE)
    }
    sim_chunk <- sim_chunk[, !rm_rows[-1]] # rm first column of scans number
    
    # if any joined id is inside the read chunk, aggregate the rows
    # only aggregate rows with at least one not NA value, which are the columns >= from_row
    # or with at least one not NA value
    if (any(x_rows_chunk)) 
    {
      sim_chunk[x_rows[x_rows_chunk]-from_row,] <- apply(sim_chunk, 2, function(x) {
        sapply(joined_ids[x_rows_chunk], 
               function(j) {
                 if (!all(is.na(x[j-from_row]))) {
                   max(x[j-from_row], na.rm = T)
                 } else {
                   NA
                 }})
      })
      # remove already aggregated rows
      sim_chunk <- sim_chunk[!rm_rows[(from_row+1):to_row], ]
    }
    
    # make lower triangle equals zero
    zero_matrix_tri_down(sim_chunk, from_row - n_joins, from_row - n_joins + nrow(sim_chunk)-1)
    
    # zero all columns with index <= from_row - n_joins - 1, these should be NA
    # when (from_row - n_joins - 1) > 1, this is not the first chunk and there will be NA columns left
    if ((from_row - n_joins - 1) >= 1) {
      sim_chunk[,1:(from_row - n_joins - 1)] <- 0
    }
    # check if there is still any NA value in the matrix, if yes something went wrong
    if (any(is.na(sim_chunk)) || any(is.infinite(sim_chunk))) {
      stop("ERROR when aggregating the similarity table. ",
           "There is still a NA or infinite value in the table, some indexing may have",
           " gone wrong - all values should be valid and greater or equal than zero here.",
           " Please contact the dev team.")
    }
    
    # add the scans order as the first column
    sim_chunk <- cbind(scans_order[(from_row+1):to_row][!rm_rows[(from_row+1):to_row]], 
                       sim_chunk)
    # store the current chunk of the similarity table
    write_csv(as.data.frame(sim_chunk), 
              path = file.path(output_path, 
                               "molecular_networking/similarity_tables", 
                               paste0("similarity_table_", output_name, "_tmp.csv")),
              col_names = FALSE, append = TRUE)
    rm(sim_chunk)
    if (to_row == n_scans)
      break;
  }
  cat("|\n")
  # check number of resulting rows and cols
  res <- readr::count_fields(file.path(output_path,
                                       "molecular_networking/similarity_tables",
                                       paste0("similarity_table_", output_name, "_tmp.csv")), 
                             skip = 1, tokenizer = tokenizer_csv())
  if (length(res) != nscans || any(res != nscans + 1)) {
    unlink(file.path(output_path,
                     "molecular_networking/similarity_tables",
                     paste0("similarity_table_", output_name, "_tmp.csv")),
           force = TRUE)
    unlink(file.path(output_path,
                     "molecular_networking/similarity_tables",
                     paste0("similarity_table_", output_name, "_aggMax.csv")),
           force = TRUE)
    stop("Wrong similarity table dimension after aggregation, something went wrong. Process Aborted.",
         " Please contact the dev team.")
  }
  # if number of rows is correct, then save clean table and remove the tmp
  file.copy(from = file.path(output_path,
                             "molecular_networking/similarity_tables",
                             paste0("similarity_table_", output_name, "_tmp.csv")), 
            to = file.path(output_path,
                           "molecular_networking/similarity_tables",
                           paste0("similarity_table_", output_name, "_aggMax.csv")),
            overwrite = TRUE)
  unlink(file.path(output_path,
                   "molecular_networking/similarity_tables",
                   paste0("similarity_table_", output_name, "_tmp.csv")),
         force = TRUE)
}

# function to merge the count table by column
# rules to merge each column
merge_counts <- function(col_name, x)
{
  switch(col_name,
         msclusterID = min(as.numeric(x$msclusterID)),
         numJoins =  sum(x$numJoins) + nrow(x) - 1,
         mzConsensus=,rtMean=
           round(weighted.mean(x[[col_name]], x$sumInts), mz_rt_digits),
         rtMin=ifelse(any(x$rtMin == 0 & x$rtMax == 1000000), 0, round(weighted.mean(x[[col_name]], x$sumInts), mz_rt_digits)),
         rtMax=ifelse(any(x$rtMin == 0 & x$rtMax == 1000000), 1000000, round(weighted.mean(x[[col_name]], x$sumInts), mz_rt_digits)),
         numSpectra =,BLANKS_TOTAL =,BEDS_TOTAL=,CONTROLS_TOTAL=,sumInts = sum(x[[col_name]]),
         basePeakInt=max(as.numeric(x[[col_name]])),
         BEDFLAG=,BFLAG =,CFLAG = any(as.logical(x[[col_name]])),
         HFLAG=,DESREPLICATION=,scans=,joinedJobsIDs=
           ifelse(any(!is.na(x[[col_name]])), # if there is a not NA value paste it
                  paste(x[[col_name]][!is.na(x[[col_name]])], collapse = ";"), 
                  NA),
         peaksList=,peaksInt= # if joining blank mzs, only keep the most intense spectra
           ifelse(any(x$rtMin == 0), x[[col_name]][which.max(x$sumInts)],
                  ifelse(any(!is.na(x[[col_name]])), # if there is a not NA value paste it
                         paste(x[[col_name]][!is.na(x[[col_name]])], collapse = ";"), 
                         NA)),
         peakIds= # if there is a not NA value paste it, remove duplicates
           ifelse(any(!is.na(x[[col_name]])), 
                  paste(unique(unlist(strsplit(x[[col_name]][!is.na(x[[col_name]])], ";"))), collapse = ";"),
                  NA),
         joinedIDs = paste(c(x$msclusterID[is.na(x[[col_name]])],
                             x[[col_name]][!is.na(x[[col_name]])]), collapse = ";"),
         ifelse(is.numeric(x[[col_name]]) && !startsWith(col_name, 'gnps_'), # do not sum gnps results, concatenate them
                sum(x[[col_name]], na.rm = TRUE), # sum all the counts cols (_area and _spectra)
                ifelse(any(!is.na(x[[col_name]])), # if there is a not NA value paste it, delim = ;
                       paste(x[[col_name]][!is.na(x[[col_name]])], collapse = ";"), 
                       NA))) # cocatenate string fields (e.g. from gnps)
}


# create the aggMax similarity file - copy from the clustering sim table
if (!file.exists(file.path(output_path, 
                           "molecular_networking/similarity_tables", 
                           paste0("similarity_table_", output_name, ".csv"))))
{
  stop("The pairwise similarity table file '", file.path(output_path, 
                                                         "molecular_networking/similarity_tables", 
                                                         paste0("similarity_table_", output_name, ".csv")),
       "' do not exists. ",
       "Rerun the job or provide a valid path to an output directory.")
} else {
  sim_file <- file.path(output_path, 
                        "molecular_networking/similarity_tables", 
                        paste0("similarity_table_", output_name, "_aggMax.csv"))
  if (!file.copy(file.path(output_path, 
                           "molecular_networking/similarity_tables", 
                           paste0("similarity_table_", output_name, ".csv")),
                 sim_file, overwrite = TRUE))
  {
    stop("The pairwise similarity table file '", file.path(output_path, 
                                                           "molecular_networking/similarity_tables", 
                                                           paste0("similarity_table_", output_name, ".csv")),
         "' could not be copied. Give write privileges to the current user.")
  } 
}

batch_metadata <- readMetadataTable(path_batch_metadata)

# read count files
ms_spectra_count <- suppressMessages(read_csv(path_spectra_count, 
                                              guess_max = 5000,
                                              col_types = cols(.default="?", 
                                                               msclusterID="i")))
ms_spectra_count <- ms_spectra_count[, !startsWith(names(ms_spectra_count), "tremolo_")]
ms_spectra_count$joinedIDs <- NA
ms_spectra_count$numJoins <- 0
# multiple the peakInts by the sumInts to correctly scale the joined spectra peaks
ms_spectra_count$peaksInt <- sapply(seq_len(nrow(ms_spectra_count)), function(i){
  peaksInt <- as.numeric(strsplit(ms_spectra_count$peaksInt[[i]], ";")[[1]])
  # remove scale from intensities
  peaksInt <- inverse_scale(peaksInt, scale_factor)
  # weight by the sumInts
  peaksInt <- peaksInt*ms_spectra_count$sumInts[[i]]
  paste0(peaksInt, collapse=";")
})
# disable the BFLAG and NOISE cutoff if there is no blank sample
if ('BLANKS_TOTAL' %in% names(ms_spectra_count)) {
  blanks_flag <- TRUE
} else {
  blanks_flag <- FALSE
}
# if no blanks, disable BFLAG cutoff
if (!blanks_flag) {
  BFLAG_cutoff_factor <- FALSE
  BFLAG_cutoff <- -1
  # NOISE_cutoff <- -1
  cat("\nBFLAG cutoff : disabled\n")
  # cat("\nNOISE cutoff : disabled\n")
} 
if (BFLAG_cutoff_factor != FALSE || NOISE_cutoff_factor != FALSE) {
  # compute the summary of the basePeakInt distribution for blank mzs if any, or using the complete distribution
  if (blanks_flag) {
    summary_basePeakInt <- summary(ms_spectra_count$basePeakInt[ms_spectra_count$BLANKS_TOTAL > 0])
  } else {
    summary_basePeakInt <- summary(ms_spectra_count$basePeakInt)
  }
  if (BFLAG_cutoff_factor == FALSE) {
    # if there is no blank, variable BFLAG_cutoff is already defined;
    # then do not print the bflag disabled message again
    if (!exists('BFLAG_cutoff')) {
      BFLAG_cutoff <- -1
      cat("\nBFLAG cutoff : disabled\n")
    }
  } else {
    # BFLAG_cutoff is a numeric value, compute the cutoff as the basePeakInt median + bflag_cutoff_factor*(q75-q25) of the blank mzs
    BFLAG_cutoff <- summary_basePeakInt[['Median']] + BFLAG_cutoff_factor*(summary_basePeakInt[['3rd Qu.']]-summary_basePeakInt[['1st Qu.']]) 
    cat("\nBFLAG cutoff : basePeakInt median +",BFLAG_cutoff_factor,"* (q75-q25) =", BFLAG_cutoff,"\n")
  }
  if (NOISE_cutoff_factor == FALSE) {
    NOISE_cutoff <- -1
    cat("\nNOISE cutoff : disabled\n")
  } else {
    # NOISE_cutoff is a numeric value, compute the cutoff as the basePeakInt median + NOISE_cutoff_factor*(q75-q25) of the blank mzs
    NOISE_cutoff <- summary_basePeakInt[['Median']] + NOISE_cutoff_factor*(summary_basePeakInt[['3rd Qu.']]-summary_basePeakInt[['1st Qu.']]) 
    cat("\nNOISE cutoff : basePeakInt median +",NOISE_cutoff_factor,"* (q75-q25) =", NOISE_cutoff,"\n")
  }
  remove(summary_basePeakInt)
} else {
  if (!exists('BFLAG_cutoff')) {
    BFLAG_cutoff <- -1
    cat("\nBFLAG cutoff : disabled\n")
  }
  NOISE_cutoff <- -1
  cat("\nNOISE cutoff : disabled\n")
}

not_count_columns <- which(!(names(ms_spectra_count) %in% paste0(batch_metadata$SAMPLE_CODE, "_spectra"))) 
count_columns <- which(names(ms_spectra_count) %in% paste0(batch_metadata$SAMPLE_CODE, "_spectra"))

# order in the columns and lines: skip header[1] and add -1 to avoid first column
scans_order <- c(-1, unlist(read.csv(file.path(output_path, 
                                               "molecular_networking/similarity_tables", 
                                               paste0("similarity_table_", output_name, ".csv")), 
                                     stringsAsFactors = FALSE, nrows= 1,
                                     strip.white = TRUE, header = FALSE)[-1]))
nscans <- length(scans_order) - 1 # rm first column with scan numbers

order_table <- match(scans_order[-1], ms_spectra_count$msclusterID)
if (any(is.na(order_table))) {
  stop("Wrong matching between the pairwise similarity table and the provided count table. Something went wrong in the pairwise similarity table computation.")
}
ms_spectra_count <- ms_spectra_count[order_table,]

# apply the noise cutoff filter before the clean step if the number of rows is
# greater than 15k, what can cause a slow processing
# apply the Noise cutoff filter based on the basePeakInt <= noisy_cutoff
if (NOISE_cutoff >= 0 && any(ms_spectra_count$basePeakInt <= NOISE_cutoff) && 
    nrow(ms_spectra_count) > 15000) {
  cat("\n  ** Applying the noise cutoff and removing all spectra with a basePeakInt value <=",
      NOISE_cutoff,"**\n")
  spectra_to_keep <- (ms_spectra_count$basePeakInt > NOISE_cutoff)
  # ms_area_count <- ms_area_count[spectra_to_keep,]
  ms_spectra_count <- ms_spectra_count[spectra_to_keep,]
  # also remove the spectra from the similarity table
  # consider the first column which have the scans numbers
  spectra_to_keep <- c(TRUE, spectra_to_keep)
  col_types <- rep("d", length(spectra_to_keep))
  col_types[!spectra_to_keep] <- "-"
  col_types[1] <- 'c' # first column with scan numbers as char
  col_types <- paste0(col_types, collapse = "") 
  # read pairwise sim limited by the max chunk size and write it again removing
  # the columns and the rows of the removed spectra
  rows_read <- 0
  # when the number of rows read is equal the spectra_to_keep length, finish the 
  # process - the header is also counted
  while (rows_read < length(spectra_to_keep)) {
    pairwise_sim <- read_csv(file.path(output_path, 
                                       "molecular_networking/similarity_tables", 
                                       paste0("similarity_table_", output_name, "_aggMax.csv")), 
                             n_max = table_limit_size, skip = rows_read, 
                             col_names = F, col_types = col_types)
    nrow_pariwise_sim <- nrow(pairwise_sim)
    pairwise_sim <- pairwise_sim[spectra_to_keep[(rows_read+1):min(rows_read+table_limit_size,length(spectra_to_keep))],]
    rows_read <- rows_read + nrow_pariwise_sim
    write_csv(pairwise_sim, 
              path = file.path(output_path, 
                               "molecular_networking/similarity_tables", 
                               paste0("similarity_table_", output_name, "_tmp.csv")),
              col_names = FALSE, append = TRUE)
    
  }
  # copy tmp similarity table to the final aggMax with the removed spectra similarities
  file.copy(from = file.path(output_path,
                             "molecular_networking/similarity_tables",
                             paste0("similarity_table_", output_name, "_tmp.csv")), 
            to = file.path(output_path,
                           "molecular_networking/similarity_tables",
                           paste0("similarity_table_", output_name, "_aggMax.csv")),
            overwrite = TRUE)
  unlink(file.path(output_path,
                   "molecular_networking/similarity_tables",
                   paste0("similarity_table_", output_name, "_tmp.csv")),
         force = TRUE)
  # recompute col_types - remove columns with col_types equals '-'
  col_types <- gsub(pattern = '-', replacement = '', x = col_types)
  # order in the columns and lines: skip header[1] and add -1 to avoid first column
  scans_order <- scans_order[spectra_to_keep]
  nscans <- length(scans_order) - 1
  cat("\n  * Done removing", sum(!spectra_to_keep), "noise spectra  *\n\n")
  remove(spectra_to_keep, pairwise_sim, rows_read, nrow_pariwise_sim)
}

if (nscans != nrow(ms_spectra_count)) {
  stop("Wrong pairwise similarity table. It is not compatible with the provided count table.")
}
# store the scan number of the spectra after the clean step joins
assigned_scans <- scans_order 

cat("\n** Cleanning counts, removing redudancies and aggregating data **\n")
if (!dir.exists(file.path(output_path, "count_tables", "clean")))
  dir.create(file.path(output_path, "count_tables", "clean"), showWarnings = FALSE)
t0 <- Sys.time()
# if there is a blank sample, run a previous step only joining the bflags to the blank clusters
if (blanks_flag) {
  step_join <- 0
} else {
  step_join <- 1
}
num_joins_total <- 0
# stores the number of joins by step and the total number of joins by cluster
num_joins_last_step <- total_num_join_clusters <- ms_spectra_count$numJoins

repeat
{
  t1 <- Sys.time()
  cat("\n      * Step", step_join, "*\n")
  num_joins <- 0
  
  # read the pairwise table in chunks
  col_types <- strrep("d", length(scans_order))
  if (step_join <= 1)
  {
    # read pairwise sim limited by the max chunk size
    pairwise_sim <- read_csv(file.path(output_path, 
                                       "molecular_networking/similarity_tables", 
                                       paste0("similarity_table_", output_name, "_aggMax.csv")), 
                             n_max = table_limit_size, skip = 1, 
                             col_names = F, col_types = col_types)
    
    if (step_join == 1 && blanks_flag && num_joins_total > 0) {
      # started with step 0, then if there was at least one joining and
      # num_joins_total > 0 -> read count tables from last step
      # otherwise continue to use the clustering table in step 1
      ms_spectra_count <- suppressMessages(read_csv(file.path(output_path, "count_tables", "clean", 
                                                              paste0(output_name,"_spectra_clean.csv")),
                                                    guess_max = 5000,
                                                    col_types = cols(.default="?", 
                                                                     msclusterID="i")))
    }
  } else { 
    # after first step read the lines of all joined clusters
    # also use similarity of column info -> to obtain info of preceding scans
    if (length(joined_ids_step) > table_limit_size)
      joined_ids <- joined_ids_step[1:table_limit_size]
    else
      joined_ids <- joined_ids_step
    
    pairwise_sim <- read_csv_chunked(file.path(output_path, 
                                               "molecular_networking/similarity_tables", 
                                               paste0("similarity_table_", output_name, "_aggMax.csv")), 
                                     DataFrameCallback$new(function(x, pos) 
                                       subset(x, unlist(x[1]) %in% joined_ids)), 
                                     col_names = TRUE, col_types = col_types)
    # get preceding scans sim info
    pairwise_sim[,-1] <- pairwise_sim[,-1] + t(read_csv(file.path(output_path, 
                                                                  "molecular_networking/similarity_tables", 
                                                                  paste0("similarity_table_", output_name, "_aggMax.csv")), 
                                                        col_names = TRUE, paste(sapply(scans_order, 
                                                                                       function(i, x) ifelse(i %in% x, "d", "-"), 
                                                                                       joined_ids), 
                                                                                collapse = "")))
    
    # read count tables from last step
    ms_spectra_count <- suppressMessages(read_csv(file.path(output_path, "count_tables", "clean", 
                                                            paste0(output_name,"_spectra_clean.csv")),
                                                  guess_max = 5000,
                                                  col_types = cols(.default="?", 
                                                                   msclusterID="i")))
  }
  
  # print(table_limit_size)
  progress_joins <- unique(trunc(c(seq(from = 1, to = nscans, by = nscans/25), nscans)))
  n_progress <- length(progress_joins)
  cat(paste0("  |", paste0(rep(" ", n_progress), collapse = ""), "|\n  |", collapse = ""))
  
  i <- 1
  while (i <= nscans) 
  {
    #cat("i: ", i, "\n")
    if (n_progress > 0 && i == progress_joins[[1]]) {
      progress_joins <- progress_joins[-1]
      n_progress <- n_progress - 1
      cat("=")
    }
    
    # get next cluster scan number and cluster info
    cluster <- ms_spectra_count[i,]
    scan_num <- cluster[[1]]
    
    # after step 1 only check clusters that changed - that were joined
    if (step_join > 1 && cluster$numJoins == 0) 
    { 
      i <- i + 1
      next()
    } else if (step_join == 0 && cluster$BLANKS_TOTAL == 0) 
    {
      # if step_join == 0 and cluster not blank, go to next
      # only join blank clusters in step 0
      i <- i + 1
      next()
    }
    
    sim_i <- which_eq(pairwise_sim[[1]], scan_num, 1)
    if (length(sim_i) == 0) # read more lines to find scan_num
    {
      if (step_join <= 1)
      {
        sim_i <- which_eq(scans_order, scan_num, 1) # scan line
        pairwise_sim <- read_csv(file.path(output_path, 
                                           "molecular_networking/similarity_tables", 
                                           paste0("similarity_table_", output_name, "_aggMax.csv")), 
                                 n_max = table_limit_size, skip = sim_i-1, 
                                 col_names = F, col_types = col_types)
        
        sim_i <- which_eq(pairwise_sim[[1]], scan_num, 1)
        if (length(sim_i) == 0) # read more lines
          stop("Error with the new similarity table, wrong scans order.")
        
      } else { # after first step only read joined ids similarity
        joined_ids <- joined_ids_step[joined_ids_step >= scan_num]
        if (length(joined_ids) > table_limit_size)
          joined_ids <- joined_ids[1:table_limit_size]
        
        pairwise_sim <- read_csv_chunked(file.path(output_path, 
                                                   "molecular_networking/similarity_tables", 
                                                   paste0("similarity_table_", output_name, "_aggMax.csv")), 
                                         DataFrameCallback$new(function(x, pos) 
                                           subset(x, unlist(x[1]) %in% joined_ids)), 
                                         col_names = TRUE, col_types = col_types)
        # get preceding scans sim info
        pairwise_sim[,-1] <- pairwise_sim[,-1] + t(read_csv(file.path(output_path, 
                                                                      "molecular_networking/similarity_tables", 
                                                                      paste0("similarity_table_", output_name, "_aggMax.csv")), 
                                                            col_names = TRUE, paste(sapply(scans_order, 
                                                                                           function(i, x) ifelse(i %in% x, "d", "-"), 
                                                                                           joined_ids), 
                                                                                    collapse = "")))
        sim_i <- which_eq(pairwise_sim[[1]], scan_num, 1)
        if (length(sim_i) == 0) # read more lines
          stop("Error with the new similarity table, wrong scans order.")
      }
    }
    
    if (pairwise_sim[sim_i,1] != scan_num)
      stop("Wrong scans order in the pairwise similarity.")
    
    # get the clusters ids that have a similarity with the current cluster including itself
    adj_clusters <- assigned_scans[which_ge(unlist(pairwise_sim[sim_i,-1]), sim_tol, 1)] # add one to scape first column
    
    # get the spectra that are in the same peak as the current cluster:
    # mz diff <= mz_tol, 
    # the rt center of one spectrum is contained in the other spectrum rt range within the rt_tol,
    cluster_peak <- ms_spectra_count[abs(ms_spectra_count$mzConsensus - cluster$mzConsensus) <= mz_tol &
                                       ((ms_spectra_count$rtMean >= cluster$rtMin - rt_tol & 
                                           ms_spectra_count$rtMean <= cluster$rtMax + rt_tol) |
                                          (cluster$rtMean >= ms_spectra_count$rtMin - rt_tol &
                                             cluster$rtMean <= ms_spectra_count$rtMax + rt_tol)),]
    
    # if not bflag check peak center and boundaries deviation, remove peaks not aligned
    # the spectra peak center deviation is <= 4 * rt_tol or peak boundaries deviation <= 2*rt_tol
    if (!blanks_flag || !cluster$BFLAG) {
      cluster_peak <- cluster_peak[abs(cluster$rtMean-cluster_peak$rtMean) <= 4*rt_tol |
                                     apply(cluster_peak[,c("rtMin", "rtMax")], 1,
                                           function(x,y)  RMSE(x, y),
                                           y = c(cluster$rtMin, cluster$rtMax)) <= 2*rt_tol,]
    }
    
    # if there is a not similar cluster in the current cluster peak check if they share any MS1 peak id 
    # and add them to the adj list to be joined
    if (cluster$numSpectra < 5000) { # too many spectra will make the regular expression too big 
      non_adj_peak <- !(cluster_peak$msclusterID %in% adj_clusters)
      if (any(non_adj_peak)) {
        # do not consider fake peaks ids
        peakIds <- strsplit(cluster$peakIds,";")[[1]]
        # peaksIds <- paste0(peakIds,collapse = "|")
        peakIds <- peakIds[!startsWith(peakIds, "fake_")]
        if (length(peakIds) > 0 && length(peakIds) <= 500) {
          peaksIds <- paste0(peakIds, collapse = "|")
          non_adj_peak <- cluster_peak$msclusterID[non_adj_peak][
            grepl(pattern = peaksIds,  
                  cluster_peak$peakIds[non_adj_peak])]
          if (length(non_adj_peak) > 0)
            adj_clusters <- c(adj_clusters, non_adj_peak)
        }
      }
    }
    # remove bflags TRUE using the basePeakInt cutoff
    if (BFLAG_cutoff >= 0) {
      non_adj_peak <- NULL
      # if blank cluster, also join BFLAGS with base peak below the cutoff
      if (cluster$BLANKS_TOTAL > 0) {
        non_adj_peak <- ((!(cluster_peak$msclusterID %in% adj_clusters)) & 
                           (cluster_peak$basePeakInt <= BFLAG_cutoff))
      } else if (cluster$BFLAG && cluster$basePeakInt <= BFLAG_cutoff) {
        # if not blank cluster but BFLAG True and basePeakInt <= cutoff, also join to blank clusters in the peak
        non_adj_peak <- ((!(cluster_peak$msclusterID %in% adj_clusters)) & 
                           (cluster_peak$BLANKS_TOTAL > 0))
        
      }
      if (any(non_adj_peak)) {
        adj_clusters <- c(adj_clusters, cluster_peak$msclusterID[non_adj_peak])
      }
    }
    # filter only adj peaks
    cluster_peak <- cluster_peak[cluster_peak$msclusterID %in% adj_clusters,]
    
    if (nrow(cluster_peak) == 1) # just the current cluster, then go to next spectrum
    {
      i <- i + 1
      next()
    } else if (nrow(cluster_peak) == 0) 
    {
      stop("Error in the similarity table aggMax, the current cluster is not ",
           "adjacent to itself - diagonal probable different from 1.0. ",
           "Something went wrong in the similarity table aggMax construction.")
    }
    num_joins <- num_joins + nrow(cluster_peak) - 1
    
    # get the cluster_peak members pos in the count table
    i_count <- match(cluster_peak$msclusterID, ms_spectra_count$msclusterID)
    
    # merge counts based on spectra counts
    cluster <- lapply(names(cluster_peak), merge_counts, cluster_peak)
    ms_spectra_count[i,] <- cluster
    
    # remove merged row from count tables
    i_count <- i_count[i_count != i] # remove the cluster pos from the cluster_peak members
    ms_spectra_count <- ms_spectra_count[-i_count,] 
    
    num_joins_last_step <- num_joins_last_step[-i_count]
    total_num_join_clusters <- total_num_join_clusters[-i_count]
    
    # update assigned scans of joined cluster
    assigned_scans[assigned_scans %in% cluster_peak$msclusterID] <- cluster[[1]]
    
    nscans <- nscans - nrow(cluster_peak) + 1 # remove scans that were merged from the total count
    if (n_progress > 0 && progress_joins[n_progress] > nscans) # rm from progress the joined scans
    {
      cat("=")
      progress_joins <- progress_joins[-n_progress]
      n_progress <- n_progress - 1
    }
    
    # decrement number of rows preceding the current that were removed - prevent skipping rows
    i <- i + 1 - sum(i_count < i)
  }
  cat("|\n")
  rm(pairwise_sim)
  
  # reset number of joins by m/z with the joined clusters of last step 
  ms_spectra_count$numJoins <- ms_spectra_count$numJoins - num_joins_last_step
  
  # reset number of joined clusters by m/z in this step
  num_joins_last_step <- ms_spectra_count$numJoins
  # update number of total joins by cluster after this step
  total_num_join_clusters <- total_num_join_clusters + num_joins_last_step
  
  t2 <- Sys.time()
  cat("        * Joined", num_joins, "similar clusters in", 
      round(t2-t1, 2), units(t2-t1), "*\n")
  
  # aggregate similarity values of joined idxs
  if (num_joins > 0)
  {
    cat("          * Aggregating similarity table with the joined Ids *\n")
    # aggregate the sim table row and cols, get the duplicated ids (to be removed) 
    # and match the joined ids with their position
    joined_idx <- duplicated(assigned_scans)
    numJoins_idx <- (ms_spectra_count$numJoins > 0)
    joined_ids <- lapply(ms_spectra_count$joinedIDs[numJoins_idx], 
                         function(x) {
                           last_joined_scans <- match(as.numeric(strsplit(x, split = ";")[[1]]), scans_order)
                           last_joined_scans <- last_joined_scans[!is.na(last_joined_scans)] # remove NA introduced due to already joined scans from last step
                           last_joined_scans})
    names(joined_ids) <- match(ms_spectra_count$msclusterID[numJoins_idx],
                               scans_order)
    
    # free memory space before next step and the sim table aggregation
    joined_ids_step <- ms_spectra_count$msclusterID[numJoins_idx]
    write_csv(ms_spectra_count, path = file.path(output_path, "count_tables", "clean", 
                                                 paste0(output_name, "_spectra_clean.csv")))
    rm(ms_spectra_count, numJoins_idx)
    
    aggregate_sim_table(joined_ids, scans_order, joined_idx, col_types, 
                        nscans, sim_file, table_limit_size)
    
    # remove duplicated from assigned scans and from scans_order
    scans_order <- scans_order[!joined_idx]
    assigned_scans <- assigned_scans[!joined_idx]
    t1 <- Sys.time()
    cat("          * Done aggregating similarity table in", 
        round(t1-t2, 2), units(t1-t2), "*\n")
  }
  
  num_joins_total <- num_joins_total + num_joins
  
  # stop joining if no join was made in the last step and this is not the blanks step (step_join == 0)
  if (step_join >= 1 && (num_joins == 0 || step_join == 10))
    break()
  
  step_join <- step_join + 1
}

# get total number of joins by cluster
ms_spectra_count$numJoins <- total_num_join_clusters
# scale the spectra peaksInt by dividing by the sumInts
ms_spectra_count$peaksInt <- sapply(seq_len(nrow(ms_spectra_count)), function(i){
  peaksInt <- as.numeric(strsplit(ms_spectra_count$peaksInt[[i]], ";")[[1]])
  peaksInt <- peaksInt/ms_spectra_count$sumInts[[i]]
  # rescale the peaksInt of the not joined clusters
  if (is.na(ms_spectra_count$joinedIDs[[i]])) {
    peaksInt <- scaleInts(peaksInt, scale_factor)
  }
  paste0(peaksInt, collapse=";")
})

tf <- Sys.time()
cat("\n  * Done reducing", num_joins_total, "similar clusters in", 
    round(tf-t0, 2), units(tf-t0), "*\n\n")
rm(cluster, cluster_peak, total_num_join_clusters, num_joins, num_joins_total, 
   num_joins_last_step)
# warnings()
cat("\n  ** Computing Peak Area by Sample ** \n\n")
ms_area_count <- ms_spectra_count
names(ms_area_count)[count_columns] <- sub(pattern = "_spectra",
                                           replacement = "_area", fixed = TRUE,
                                           x = names(ms_area_count)[count_columns])
if (!joining_jobs) 
{
  # this is the default flow of np3, compute the peak areas using the provided
  # pre processed data path
  peak_areas_base_peak_int <- compute_peak_area(processed_data_path,
                                                ms_area_count$msclusterID,
                                                lapply(ms_area_count$scans, function(x) strsplit(x, ";")[[1]]),
                                                lapply(ms_area_count$peakIds, function(x) strsplit(x, ";")[[1]]),
                                                batch_metadata)
} else {
  # this is the join_jobs flow, compute the peak areas using all the original samples
  # pre processed data
  peak_areas_base_peak_int <- compute_peak_areas_joined_jobs(output_path, 
                                                             ms_area_count$msclusterID,
                                                             ms_area_count$scans,
                                                             ms_area_count$peakIds, 
                                                             ms_area_count$joinedJobsIDs)
}
# assign the peak areas following the count columns order 
# and also add the base peak intensity (last column)
ms_area_count[,count_columns] <- peak_areas_base_peak_int[,
                                                          match(names(ms_area_count)[count_columns],
                                                                names(peak_areas_base_peak_int)[-ncol(peak_areas_base_peak_int)])]
ms_spectra_count$basePeakInt <- ms_area_count$basePeakInt <- peak_areas_base_peak_int$basePeakInt
rm(peak_areas_base_peak_int)

# compute max area and the mean precursor intensity of the final clusters
ms_area_count$maxArea <- ms_spectra_count$maxArea <- apply(ms_area_count[,count_columns], 1, max)
ms_area_count$meanInt <- ms_spectra_count$meanInt <- ms_area_count$sumInts / ms_area_count$numSpectra

#stop("Check sim aggMax!!")

# apply the Noise cutoff filter based on the basePeakInt <= noisy_cutoff
if (NOISE_cutoff >= 0 && any(ms_spectra_count$basePeakInt <= NOISE_cutoff)) {
  order_table <- match(scans_order[-1], ms_spectra_count$msclusterID)
  if (any(is.na(order_table))) {
    stop("Wrong matching between the pairwise similarity table and the provided count table. Something went wrong in the pairwise similarity table computation.")
  }
  ms_spectra_count <- ms_spectra_count[order_table,]
  ms_area_count <- ms_area_count[order_table,]
  cat("\n  ** Applying the noise cutoff and removing all spectra with a basePeakInt value <=",
      NOISE_cutoff,"**\n")
  spectra_to_keep <- (ms_spectra_count$basePeakInt > NOISE_cutoff)
  ms_area_count <- ms_area_count[spectra_to_keep,]
  ms_spectra_count <- ms_spectra_count[spectra_to_keep,]
  # also remove the spectra from the similarity table
  # consider the first column which have the scans numbers
  spectra_to_keep <- c(TRUE, spectra_to_keep)
  col_types <- rep("d", length(spectra_to_keep))
  col_types[!spectra_to_keep] <- "-"
  col_types[1] <- 'c'
  col_types <- paste0(col_types, collapse = "") 
  # read pairwise sim limited by the max chunk size and write it again removing
  # the columns and the rows of the removed spectra
  rows_read <- 0
  # when the number of rows read is equal the spectra_to_keep length, finish the 
  # process - the header is also counted
  while (rows_read < length(spectra_to_keep)) {
    pairwise_sim <- read_csv(file.path(output_path, 
                                       "molecular_networking/similarity_tables", 
                                       paste0("similarity_table_", output_name, "_aggMax.csv")), 
                             n_max = table_limit_size, skip = rows_read, 
                             col_names = F, col_types = col_types)
    nrow_pariwise_sim <- nrow(pairwise_sim)
    pairwise_sim <- pairwise_sim[spectra_to_keep[(rows_read+1):min(rows_read+table_limit_size,length(spectra_to_keep))],]
    rows_read <- rows_read + nrow_pariwise_sim
    write_csv(pairwise_sim, 
              path = file.path(output_path, 
                               "molecular_networking/similarity_tables", 
                               paste0("similarity_table_", output_name, "_tmp.csv")),
              col_names = FALSE, append = TRUE)
    
  }
  # copy tmp similarity table to the final aggMax with the removed spectra similarities
  file.copy(from = file.path(output_path,
                             "molecular_networking/similarity_tables",
                             paste0("similarity_table_", output_name, "_tmp.csv")), 
            to = file.path(output_path,
                           "molecular_networking/similarity_tables",
                           paste0("similarity_table_", output_name, "_aggMax.csv")),
            overwrite = TRUE)
  unlink(file.path(output_path,
                   "molecular_networking/similarity_tables",
                   paste0("similarity_table_", output_name, "_tmp.csv")),
         force = TRUE)
  # recompute col_types - remove columns with col_types equals '-'
  col_types <- gsub(pattern = '-', replacement = '', x = col_types)
  # order in the columns and lines: skip header[1] and add -1 to avoid first column
  scans_order <- scans_order[spectra_to_keep]
  cat("\n  * Done removing", sum(!spectra_to_keep), "noise spectra  *\n\n")
  remove(spectra_to_keep, pairwise_sim, rows_read, nrow_pariwise_sim)
}

cat("\n  ** Checking joined ids, recomputing samples types indicators and aggregating peak list of joined ids **\n")

# check if joined ID are ok: min ID is in msclusterID column and the others IDs are not
# and aggregate peaksList and peaksInt by ordering and joining adjacent peaks
joined_idx <- which(!is.na(ms_area_count$joinedIDs))
ms_spectra_count[joined_idx, c("peaksList", "peaksInt")] <- 
  ms_area_count[joined_idx, c("peaksList", "peaksInt")] <- 
  Reduce(rbind, lapply(joined_idx, function(i) {
    # check joined scans are ok  
    scan_nums <-  as.numeric(strsplit(ms_area_count$joinedIDs[[i]], split = ";")[[1]])
    min_scan <- min(scan_nums)
    scan_nums <- scan_nums[scan_nums != min_scan]
    
    if (!(min_scan %in% ms_area_count$msclusterID))
      stop("Error in the cleanning step. Missing scan number in the resulting counts table.")
    if (any(scan_nums %in% ms_area_count$msclusterID))
      stop("Error in the cleanning step. Cleaned scan number is still present in the counts table.")
    
    # join adj peaks
    peaks <- as.numeric(strsplit(ms_area_count$peaksList[[i]], split = ";")[[1]])
    ints <- as.numeric(strsplit(ms_area_count$peaksInt[[i]], split = ";")[[1]])
    # order peaks
    peaksOrder <- order(peaks)
    peaks <- peaks[peaksOrder]
    ints <- ints[peaksOrder]
    
    sapply(joinAdjPeaksScalee(peaks, ints, bin_size, -1, scale_factor, 0), 
           function(x) {
             x <- round(x,5)
             paste(x,collapse = ";")
           })
  }), init = NULL)
# check if the inverse scaled intensities sum 1000
check_ints_sum <- sapply(ms_area_count$peaksInt[joined_idx], function(x) {
  ints <- as.numeric(strsplit(x, ";")[[1]])
  # remove scale from intensities
  ints <- inverse_scale(ints, scale_factor)
  round(sum(ints)) == 1000
})
if (!all(check_ints_sum)) {
  stop("Bad scaling of the peak lists' intensities. Some spectra do not have the inverse scales intensities summing 1000 (the normalized value).")
}
cat("\n  ** Saving the clean MGF file **\n\n")
# sort by msclsuterID, this order will be applied to the scans in the MGF file
ms_area_count <- arrange(ms_area_count, msclusterID)
ms_spectra_count <- arrange(ms_spectra_count, msclusterID)
# write the peaks list and ints to the mgf file - the intensities will be 
# inversed scaled to be saved with no scaling
writeMgfDataFile_NP3_table(ms_count_table = ms_area_count,
                           file_MGF = file.path(output_path, "mgf",
                                                paste0(output_name, "_clean.mgf")),
                           output_name = output_name,
                           charge = ion_mode,
                           scale_factor = scale_factor)

# compute Blanks total and controls total and beds total
blanks_code <- batch_metadata[batch_metadata$SAMPLE_TYPE == "blank", "SAMPLE_CODE"]
if (length(blanks_code) > 0)
{
  if (length(blanks_code) == 1) {
    ms_area_count$BLANKS_TOTAL <- ms_area_count[,  paste0(blanks_code,"_area")][[1]]
    ms_spectra_count$BLANKS_TOTAL <- ms_spectra_count[,  paste0(blanks_code,"_spectra")][[1]]
  } else {
    ms_area_count$BLANKS_TOTAL <- rowSums(ms_area_count[, paste0(blanks_code,"_area")])
    ms_spectra_count$BLANKS_TOTAL <- rowSums(ms_spectra_count[, paste0(blanks_code,"_spectra")])
  }
}
controls_code <- batch_metadata[batch_metadata$SAMPLE_TYPE == "control", "SAMPLE_CODE"]
if (length(controls_code) > 0)
{
  if (length(controls_code) == 1)
  {
    ms_area_count$CONTROLS_TOTAL <- ms_area_count[, paste0(controls_code,"_area")][[1]]
    ms_spectra_count$CONTROLS_TOTAL <- ms_spectra_count[, paste0(controls_code,"_spectra")][[1]]
  } else {
    ms_area_count$CONTROLS_TOTAL <- rowSums(ms_area_count[, paste0(controls_code,"_area")])
    ms_spectra_count$CONTROLS_TOTAL <- rowSums(ms_spectra_count[, paste0(controls_code,"_spectra")])
  }
  rm(controls_code)
}
bed_controls_code <- batch_metadata[batch_metadata$SAMPLE_TYPE == "bed", "SAMPLE_CODE"]

if (length(bed_controls_code) > 0)
{
  if (length(bed_controls_code) == 1)
  {
    ms_area_count$BEDS_TOTAL <- ms_area_count[, paste0(bed_controls_code,"_area")][[1]]
    ms_spectra_count$BEDS_TOTAL <- ms_spectra_count[, paste0(bed_controls_code,"_spectra")][[1]]
  } else {
    ms_area_count$BEDS_TOTAL <- rowSums(ms_area_count[, paste0(bed_controls_code,"_area")])
    ms_spectra_count$BEDS_TOTAL <- rowSums(ms_spectra_count[, paste0(bed_controls_code,"_spectra")])
  }
  rm(bed_controls_code)
}

# remove duplicated entries due to merge
if ("DESREPLICATION" %in% names(ms_area_count))
{
  hflag_rows <- !is.na(ms_area_count$HFLAG) | !is.na(ms_area_count$DESREPLICATION)
  ms_spectra_count[hflag_rows, c("DESREPLICATION", "HFLAG")] <- 
    ms_area_count[hflag_rows, c("DESREPLICATION", "HFLAG")] <- Reduce(rbind, 
                                                                      lapply(which(hflag_rows), 
                                                                             function(i)
                                                                             {
                                                                               if (is.na(ms_area_count$HFLAG[[i]]))
                                                                                 hflag <- NA
                                                                               else
                                                                                 hflag <- paste0(unique(strsplit(ms_area_count$HFLAG[[i]], ";")[[1]]), collapse = ";")
                                                                               if (is.na(ms_area_count$DESREPLICATION[[i]]))
                                                                                 desrep <- NA
                                                                               else
                                                                                 desrep <- paste0(unique(strsplit(ms_area_count$DESREPLICATION[[i]], ";")[[1]]), collapse = ";")
                                                                               
                                                                               return(c(desrep, hflag))
                                                                             }), init = NULL)
}

cat("\n  ** Adding blanks neighborhood information **\n")
# compute distance to a blank
if (length(blanks_code) > 0)
{
  any_blank <- TRUE
  # free memmory of one count table
  write_csv(ms_spectra_count, path = file.path(output_path, "count_tables", "clean", 
                                               paste0(output_name, "_spectra_clean.csv")))
  rm(ms_spectra_count)
  
  # explore blank distS
  blanks_neighbor <- list()
  ms_area_count$BLANK_DIST <- NA
  # get directly connected nodes to a blank
  blank_ids <- ms_area_count$msclusterID[ms_area_count$BLANKS_TOTAL > 0]
  ms_area_count$BLANK_DIST[ms_area_count$BLANKS_TOTAL > 0] <- 0 # set blanks dist = 0
  
  # dist 1
  # get blanks similarity rows
  pairwise_sim_blanks <- read_csv_chunked(file.path(output_path, 
                                                    "molecular_networking/similarity_tables", 
                                                    paste0("similarity_table_", output_name, "_aggMax.csv")), 
                                          DataFrameCallback$new(function(x, pos) subset(x, unlist(x[1]) %in% blank_ids)), 
                                          col_names = TRUE, col_types = col_types)
  # use similarity proportional to the neighbor similarity* their similarity to each neighbor
  blanks_neighbor[[1]] <- which(colSums(pairwise_sim_blanks[,-1] >= sim_tol) > 0)
  # get blanks similarity cols
  pairwise_sim_blanks <- read_csv(file.path(output_path, 
                                            "molecular_networking/similarity_tables", 
                                            paste0("similarity_table_", output_name, "_aggMax.csv")), 
                                  col_names = TRUE, col_types = paste(sapply(scans_order, 
                                                                 function(i, x) ifelse(i %in% x, "d", "-"), 
                                                                 blank_ids), 
                                                          collapse = ""))
  blanks_neighbor[[1]] <- unique(c(blanks_neighbor[[1]],
                                   which(rowSums(pairwise_sim_blanks >= sim_tol) > 0)))
  blanks_neighbor[[1]] <- blanks_neighbor[[1]][!(blanks_neighbor[[1]] %in% match(blank_ids, scans_order[-1]))]
  # set blank dist to 1 to not blank nodes
  ms_area_count$BLANK_DIST[is.na(ms_area_count$BLANK_DIST) & 
                             ms_area_count$msclusterID %in% scans_order[-1][blanks_neighbor[[1]]]] <- 1
  
  for (depth in 2:blank_depth) 
  {
    blank_idx <- blanks_neighbor[[depth-1]]
    if (length(blank_idx) == 0) # if no new neighbor -> stop
      break()
    # get blanks similarity rows
    pairwise_sim_blanks <- read_csv_chunked(file.path(output_path, 
                                                      "molecular_networking/similarity_tables", 
                                                      paste0("similarity_table_", output_name, "_aggMax.csv")), 
                                            DataFrameCallback$new(function(x, pos) 
                                              subset(x, unlist(x[1]) %in% scans_order[-1][blank_idx])), 
                                            col_names = TRUE, col_types = col_types)
    blanks_neighbor[[depth]] <- which(colSums(pairwise_sim_blanks[,-1] >= sim_tol) > 0)
    # get blanks similarity cols
    pairwise_sim_blanks <- read_csv(file.path(output_path, 
                                              "molecular_networking/similarity_tables", 
                                              paste0("similarity_table_", output_name, "_aggMax.csv")), 
                                    col_names = TRUE, col_types = paste0("-", paste(sapply(seq_along(scans_order[-1]), 
                                                                               function(i, x) ifelse(i %in% x, "d", "-"), 
                                                                               blank_idx), 
                                                                        collapse = "")))
    blanks_neighbor[[depth]] <- unique(c(blanks_neighbor[[depth]],
                                         which(rowSums(pairwise_sim_blanks >= sim_tol) > 0)))
    blanks_neighbor[[depth]] <- blanks_neighbor[[depth]][!(blanks_neighbor[[depth]] %in% unlist(blanks_neighbor[1:(depth-1)]))]
    # set blank dist to depth to not assigned dists
    ms_area_count$BLANK_DIST[is.na(ms_area_count$BLANK_DIST) &
                               ms_area_count$msclusterID %in% scans_order[-1][blanks_neighbor[[depth]]]] <- depth
  }
  
  rm(blank_ids, blank_idx, blank_depth, blanks_neighbor, blanks_code, pairwise_sim_blanks)
  ms_spectra_count <- suppressMessages(read_csv(file.path(output_path, "count_tables", "clean",
                                                          paste0(output_name,"_spectra_clean.csv")),
                                                guess_max = 5000,
                                                col_types = cols(.default="?", 
                                                                 msclusterID="i")))
  
  ms_spectra_count$BLANK_DIST <- ms_area_count$BLANK_DIST
} else {
  any_blank <- FALSE
}

# sort by msclsuterID
ms_area_count <- arrange(ms_area_count, msclusterID)
ms_spectra_count <- arrange(ms_spectra_count, msclusterID)

# round mzConsensus and rts
ms_spectra_count$mzConsensus <- ms_area_count$mzConsensus <- round(ms_area_count$mzConsensus, mz_rt_digits)
ms_spectra_count$rtMean <- ms_area_count$rtMean <- round(ms_area_count$rtMean, mz_rt_digits)
ms_spectra_count$rtMin <- ms_area_count$rtMin <- round(ms_area_count$rtMin, mz_rt_digits)
ms_spectra_count$rtMax <- ms_area_count$rtMax <- round(ms_area_count$rtMax, mz_rt_digits)

# write cleanned data without annotation
write_csv(ms_area_count, path = file.path(output_path, "count_tables", "clean", 
                                          paste0(output_name, "_peak_area_clean.csv")))
write_csv(ms_spectra_count, path = file.path(output_path, "count_tables", "clean", 
                                             paste0(output_name, "_spectra_clean.csv")))
rm(ms_spectra_count, ms_area_count)
# not removing aggMax sim table - it is used in the testing scripts
# unlink(file.path(output_path,
#                  "molecular_networking/similarity_tables",
#                  paste0("similarity_table_", output_name, "_aggMax.csv")),
#        force = TRUE)
t0 <- Sys.time()
cat("\n    * Done in", round(t0-tf, 2), units(t0-tf), "*\n\n")