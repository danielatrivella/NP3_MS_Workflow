# pre-processing step
# Merge LC-MS into MS/MS tandem mass data, use xcms CentWave algorithm to detect 
# LC-MS peaks and MsnBase to read MS/MS in tandem and merge peak info
# import functions to write to mgf the MS/MS data with LC-MS peak info
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

suppressPackageStartupMessages(library(dplyr))
source(file.path(script_path(),"../read_metadata_table.R"))
source(file.path(script_path(),"pre_process_statistics.R"))
#source(file.path(script_path(),"writeMgfData_NP3.R"))
Rcpp::sourceCpp(file.path(script_path(),
                          '../read_mgf_peak_list_R.cpp'))
# TODO add arg
parallel_cores <- 1

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  stop("Three arguments must be supplied to extract the chromatograms of raw data collections:\n", 
       " 1 - The job name.",
       " 2 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and sample class;\n",
       " 3 - Path to the raw data folder;\n",
       " 4 - Retention time tolerance in seconds - used to simulate a chromatography peak when no match is found;\n",
       " 5 - Precursor mass tolerance in Da;\n",
       " 6 - Precursor mass tolerance in ppm;\n",
       " 7 - Ion mode coresponding to a ion adduct type: '1' = [M+H]+ , '0' = [M]+-, or '-1' = [M-H].\n",
       call. = FALSE)
} else {
  # default values
  rt_tol <- 3   # tolerance in Rt in seconds
  mz_tol <- 0.05   # tolerance in mz in Da
  ppm_tol   <- 15   # ppm tol - typically  set  to  a  generous multiple  of  the  mass  accuracy  of  the  mass  spectrometer
  peak_width <- c(2,10)   #  Chromatographic  peak  width  range  wmin, wmax in seconds
  snthresh  <- 0 #0   #Signal to noise ratio threshold 
  prefilter <- c(1,750) # 1,750 # only keep mzs that have at least x peaks with intensity y : (x,y)
  min_fraction <- 0.0 #0.3    # defining the minimum fraction of samples in at least one sample group in which the peaks have to be present to be considered as a peak group (feature).
  bw_win <- 2    # defining the bandwidth (standard deviation ot the smoothing kernel) to be used. This argument is passed to the [density() method.
  bin_size <- 0.05   #  defining the size of the overlapping slices in mz dimension.
  max_features <- 100   # with the maximum number of peak groups to be identified in a single mz slice.
  noise_tol <- 500   #  allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity < noise are omitted from ROI detection).
  mzCenter_fun <- "wMeanApex3"   # Name of the function to calculate the m/z center of the chromatographic peak. Allowed are: "wMean": intensity weighted mean of the peak's m/z values, "mean": mean of the peak's m/z values, "apex": use the m/z value at the peak apex, "wMeanApex3": intensity weighted mean of the m/z value at the peak apex and the m/z values left and right of it and "meanApex3": mean of the m/z value of the peak apex and the m/z values left and right of it.
  n_samples_batch_align <- 3
  processed_data_dir <- "processed_data"
  integrate_method <- 2L   # Integration method. For integrate = 1 peak limits are found through descent on the mexican hat filtered data, for integrate = 2 the descent is done on the real data. The latter method is more accurate but prone to noise, while the former is more robust, but less exact.
  fitgauss <- FALSE   # whether or not a Gaussian should be fitted to each peak. This affects mostly the retention time position of the peak.
  ion_mode <- 1
  overwrite_preprocess <- FALSE
  
  data_name <- args[[1]]
  path_batch_metadata <- file.path(args[[2]])
  if (!file.exists(path_batch_metadata))
  {
    stop("The CSV batch metadata file '", path_batch_metadata, 
         "' do not exists. Provide a valid path to where the metadata file is located.")
  }
  path_batch_metadata <- normalizePath(path_batch_metadata)
  
  path_raw_data <- file.path(args[[3]])
  if (!dir.exists(path_raw_data))
  {
    stop("The raw data file folder '", path_raw_data, 
         "' do not exists. Provide a valid path to where the raw data is located.")
  }
  path_raw_data <- normalizePath(path_raw_data)
  
  rt_tol <- as.numeric(args[[4]]) # tolerance of Rt in seconds
  mz_tol <- as.numeric(args[[5]]) # tolerance of mz in Da
  ppm_tol   <- as.numeric(args[[6]]) # ppm tol - typically  set  to  a  generous multiple  of  the  mass  accuracy  of  the  mass  spectrometer
  ion_mode <- as.numeric(args[[7]])
  
  if (!(ion_mode %in% c(-1,1)))
    stop("The ion_mode arg must be a numeric value indicating the precursor ", 
         "ion mode. One of the following valid numeric values corresponding ", 
         "to a ion adduct type: '1' = [M+H]+ or '-1' = [M-H]-",  call. = FALSE)
  ion_mode <- paste0(abs(ion_mode), ifelse(ion_mode > 0, "+", "-"))
  
  processed_data_dir <- file.path(args[[8]])
  overwrite_preprocess <- as.logical(args[[9]])
  
  if (length(args) >= 21)
  {
    peak_width <- as.numeric(strsplit(args[[10]], ",")[[1]]) #  Chromatographic  peak  width  range  wmin, wmax in seconds
    if (length(peak_width) != 2)
    {
      stop("The peak width arg must be two numeric values separated by a comma ", 
           "with the expected approximate peak width in chromatographic space.",
           " Given as a range 'min, max' in seconds.",  call. = FALSE)
    }
    snthresh  <- as.numeric(args[[11]]) #Signal to noise ratio threshold 
    
    prefilter <- as.numeric(strsplit(args[[12]], ",")[[1]]) # only keep mzs that have at least x peaks with intensity y : (x,y)
    if (length(prefilter) != 2)
    {
      stop("The pre-filter arg must be two numeric values separated by a comma ", 
           "'k, I' specifying the prefilter step for the first analysis step ",
           "(ROI detection) of the XCMS centWave algorithm. Mass traces are ", 
           "only retained if they contain at least k peaks with intensity >= I.",  call. = FALSE)
    }
    min_fraction <- as.numeric(args[[13]])
    bw_win <- as.numeric(args[[14]]) 
    bin_size <- as.numeric(args[[15]])
    max_features <- as.numeric(args[[16]])
    noise_tol <- as.numeric(args[[17]])
    mzCenter_fun <- args[[18]]
    n_samples_batch_align <- as.numeric(args[[19]])
    
    integrate_method <- as.numeric(args[[20]])
    if (!(integrate_method %in% c(1,2)))
      stop("Wrong integration method value, retry with a value equals 1 or 2. ", 
           "For integrate = 1 peak limits are found through descent on the ", 
           "mexican hat filtered data, for integrate = 2 the descent is done ", 
           "on the real data. The latter method is more accurate but prone to ", 
           "noise, while the former is more robust, but less exact.")
    fitgauss <- as.logical(args[[21]])
  } else {
    # default values from most of the parms -  pre-process executed from run
    args[[10]] <- paste(peak_width,collapse = ",")
    args[[11]] <- snthresh
    args[[12]] <- paste(prefilter, collapse = ",")
    args[[13]] <- min_fraction
    args[[14]] <- bw_win
    args[[15]] <- bin_size
    args[[16]] <- max_features
    args[[17]] <- noise_tol
    args[[18]] <- mzCenter_fun
    args[[19]] <- n_samples_batch_align
    args[[20]] <- integrate_method
    args[[21]] <- fitgauss
  }
  names(args) <- c("data_name", "metadata", "raw_data_path", "rt_tolerance", 
                   "mz_tolerance", "ppm_tolerance", "ion_mode", 
                   "processed_data_name", "processed_data_overwrite",
                   "peak_width", 
                   "snthresh", "pre_filter", "min_fraction", "bw", "bin_size", 
                   "max_features", "noise", "mz_center_fun", 
                   "max_samples_batch_align", "integrate_method", "fit_gauss")
  prettyNum(args)
}

metadata <- readMetadataTable(path_batch_metadata, path_raw_data) 

# table header to store not fragmented mzs
dt_ms1_no_ms2_header <-  c("peakIds", "mz", "rtMean", "rtMin", 
                    "rtMax", "ms2_matches", paste0(metadata$SAMPLE_CODE, "_area"))

# create dir to store processed data
if (!dir.exists(file.path(path_raw_data, processed_data_dir))) {
  dir.create(file.path(path_raw_data, processed_data_dir))
  # write table header to store not fragmented MS1 mzs
  write.table(t(no_ms2_htable), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir,
                                "MS1_list_no_MS2.csv"),
              row.names = FALSE, col.names = FALSE)
  # and to write the fragmented mzs 
  write.table(t(no_ms2_htable), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir,
                                "MS1_list_with_MS2.csv"),
              row.names = FALSE, col.names = FALSE)
} else if (!overwrite_preprocess) {
  # dir exists and overwrite = F, then only process missing files 
  # do not overwrite already processed samples
  overwrite_preprocess <- !file.exists(file.path(path_raw_data, processed_data_dir,
                                                 paste0(metadata$SAMPLE_CODE, 
                                                        "_peak_info.mgf")))
  # if everything was already preprocessed, just exit
  if (sum(overwrite_preprocess) == 0) { 
    cat("All raw LC-MS/MS files were already pre processed!\n\n")
    if (file.exists(file.path(path_raw_data, processed_data_dir,"logPreProcessStatisticsWarning")))
    {
      if (.Platform$OS.type == "unix") {
        file.show(file.path(path_raw_data, processed_data_dir,"logPreProcessStatisticsWarning"))
      } else {
        file.show(file.path(path_raw_data, processed_data_dir,"logPreProcessStatisticsWarning"), pager="console")
      }
    } else {
      if (.Platform$OS.type == "unix") {
        file.show(file.path(path_raw_data, processed_data_dir,"logPreProcessStatistics"))
      } else {
        file.show(file.path(path_raw_data, processed_data_dir,"logPreProcessStatistics"), pager="console")
      }
      
    }
    q("no")
  }
  # if the pre processing data of at least one sample is present (it was already preprocessed),
  # skip alignment
  if (sum(overwrite_preprocess) < nrow(metadata)) {
    cat("\n* Using already pre processed files, only pre processing the missing ones. ",
        "Alignment disabled (partial pre processing). *\n\n")
    n_samples_batch_align <- 0
    missing_samples_preprocess <- paste0(metadata$SAMPLE_CODE[overwrite_preprocess], "_area")
    # do not overwrite table header to store not fragmented mzs
    # extract current table header and append the new samples that will be processed to the no_ms2_htable
    previous_no_ms2_htable <- read.csv(file.path(path_raw_data, processed_data_dir,
                                   "MS1_list_no_MS2.csv"), stringsAsFactors = FALSE, 
                                   comment.char = "", strip.white = TRUE)
    no_ms2_htable <-  c(names(previous_no_ms2_htable), 
                        missing_samples_preprocess)
    previous_no_ms2_htable[,missing_samples_preprocess] <- 0
    write.table(previous_no_ms2_htable, sep = ",", 
                file = file.path(path_raw_data, processed_data_dir,
                                                "MS1_list_no_MS2.csv"), 
                row.names = FALSE)
    rm(previous_no_ms2_htable)
    # add extra headers to the table "MS1_list_with_MS2.csv" as well, rm bflag if present
    previous_ms1_w_ms2_htable <- read.csv(file.path(path_raw_data, processed_data_dir,
                                                    "MS1_list_with_MS2.csv"), 
                                          stringsAsFactors = FALSE, 
                                          comment.char = "", strip.white = TRUE)
    previous_ms1_w_ms2_htable[,missing_samples_preprocess] <- 0
    previous_ms1_w_ms2_htable[,"BFLAG"] <- NULL
    write.table(previous_ms1_w_ms2_htable, sep = ",", 
                file = file.path(path_raw_data, processed_data_dir,
                                 "MS1_list_with_MS2.csv"), 
                row.names = FALSE)
  } else {
    # complete overwrite would be required, but it is disabled - error
    # another pre processed dir must be chosen
    stop("There is another pre process result in the provided processed_data_name ('", 
         processed_data_dir,"') inside the raw_data_path, ",
         "with a different set of samples from the current metadata. ",
         "Please chose a different pre processed directory name or remove the existing one and retry.",
         " Aborting...\n\n")
  }
  
  metadata <- metadata[overwrite_preprocess,]
} else {
  # dir exists and overwrite = TRUE
  unlink(file.path(path_raw_data, processed_data_dir), recursive = TRUE)
  dir.create(file.path(path_raw_data, processed_data_dir))
  # write table header to store not fragmented MS1 mzs
  write.table(t(no_ms2_htable), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir,
                                "MS1_list_no_MS2.csv"),
              row.names = FALSE, col.names = FALSE)
  # and to write the fragmented mzs 
  write.table(t(no_ms2_htable), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir,
                                "MS1_list_with_MS2.csv"),
              row.names = FALSE, col.names = FALSE)
}

cat("Loading packages XCMS, Spectra, MsExperiment, MsBackendMgf...\n")
suppressPackageStartupMessages(library(xcms))
suppressPackageStartupMessages(library(Spectra))
suppressPackageStartupMessages(library(MsExperiment))
library(MsBackendMgf)


# TODO function start here

# check if number of cores is valid
if (parallel_cores < 1) {
  cat("Invalid parallel_cores value (", parallel_cores, 
      "). The number of cores to be used for parallel processing must be at least 2 ",
      "or 1 for disabling parallelization. The number of cores was set to 1 by default.\n")
  warning("Invalid parallel_cores parameter value.")
  parallel_cores <- 1    
} else if (parallel_cores > 1 && !require(parallel)) {
  cat("Parallel library is not available, disabling parallelization. ",
      "The number of cores was set to 1 by default.\n")
  warning("Invalid parallel_cores parameter value.")
  parallel_cores <- 1 
} else if (parallel_cores > 1 && !is.na(detectCores(logical = FALSE)) && parallel_cores >= detectCores(logical = FALSE)) {
  cat("The number of paralell_cores =",parallel_cores,
      " is greater or equal than the number of physical cores =",detectCores(logical = FALSE),
      " in your machine. This is not recommended,",
      "the number of parallel clusters will be limited by the number of physical cores minus two =",
      detectCores(logical = FALSE)-2, ".\n")
  warning("Invalid parallel_cores parameter value.")
  parallel_cores <- detectCores(logical = FALSE)-2
}
if (parallel_cores == 1) {
  # Disable parallel processing
  register(SerialParam())
} else {
  if (.Platform$OS.type == "unix") {
      register(bpstart(MulticoreParam(parallel_cores)))
  } else {
      ## Use socket based parallel processing on Windows systems
      register(bpstart(SnowParam(parallel_cores)))
      # or test
      # For Windows/Mac/Linux safely
      #register(MulticoreParam(workers = 2))
  }
}

# check if all files exists
if (!all(file.exists(metadata$FILENAME)))
{
  stop("Could not pre process the '", data_name, " data'. The following raw data files do not exists:\n", 
       paste(metadata$FILENAME[!file.exists(metadata$FILENAME)], collapse = "\n"))
}

# check if all batches have at least one non  blank sample
if (n_samples_batch_align > 0 && length(unique(metadata[metadata$SAMPLE_TYPE != "blank", 
                                                              "DATA_COLLECTION_BATCH"])) != 
    length(unique(metadata$DATA_COLLECTION_BATCH)))
{
  warning("To compute the misalignment between the samples all data collections ",
          "must have at least one not blank sample. ",
          "Data collection batch ", paste(unique(metadata$DATA_COLLECTION_BATCH)[
            which(!(unique(metadata$DATA_COLLECTION_BATCH) %in% 
                      unique(metadata[metadata$SAMPLE_TYPE != "blank", 
                                            "DATA_COLLECTION_BATCH"])))], 
            collapse = ", "), 
          " do not have a not blank sample. Misalignment computation disabled.", 
          immediate. = TRUE, call. = TRUE)
  n_samples_batch_align <- 0
}

# write parms values used
write.table(t(data.frame(data_name = data_name, path_batch_metadata = path_batch_metadata, 
                         path_raw_data = path_raw_data, rt_tol = rt_tol, mz_tol = mz_tol,
                         ppm_tol = ppm_tol, ion_mode = ion_mode, 
                         peak_width = paste0(peak_width, collapse = ","), 
                         snthresh = snthresh, prefilter =  paste0(prefilter, collapse = ","), 
                         min_fraction = min_fraction, bw_win = bw_win, bin_size = bin_size, 
                         max_features = max_features, noise_tol = noise_tol, 
                         mzCenter_fun = mzCenter_fun, n_samples_batch_align = n_samples_batch_align,
                         processed_data_dir = processed_data_dir, integrate_method = integrate_method, 
                         fitgauss = fitgauss, overwrite_preprocess = overwrite_preprocess)), 
            sep = ",", col.names  = FALSE, 
            file =  file.path(path_raw_data, processed_data_dir, "parameters.csv"))

# reorder the metadata to put the blank samples first
metadata <- metadata[c(which(metadata$SAMPLE_TYPE == "blank"), 
                       which(metadata$SAMPLE_TYPE != "blank")),]

#### using just xcms
## Create a CentWaveParam object.
cwp <- CentWaveParam(ppm = ppm_tol, noise = noise_tol, peakwidth = peak_width, 
                     snthresh = snthresh, prefilter = prefilter, mzdiff = mz_tol, 
                     fitgauss = fitgauss, mzCenterFun = mzCenter_fun, 
                     firstBaselineCheck = FALSE, integrate = integrate_method)
# create group param for alignment
pgp <- PeakGroupsParam(minFraction = min_fraction, span=0.5, extraPeaks=5)
rm(overwrite_preprocess, integrate_method, path_batch_metadata, snthresh, prefilter,
   noise_tol, mzCenter_fun, fitgauss)

cat("*** Starting to extract peak info for", data_name, "***\n")


# write the misalignment file header
# this file will store the max diff in the batches alignment
if (n_samples_batch_align > 0 & 
    !file.exists(file.path(path_raw_data, processed_data_dir, "samples_misalignment.csv"))) {
  write.table(data.frame("batch", "samples_aligned", "rtTolMin", "rtTol1stQu","rtTolMedian", 
                         "rtTolMean", "rtTol3rdQu","rtTolMax"), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir, "samples_misalignment.csv"),
              row.names = FALSE, col.names = FALSE)
}


# write the header of the table that will store the MS2 mzs without a MS1 peak correspondence
# do not overwrite here if using previous result, check if file exists
if (!file.exists(file.path(path_raw_data, processed_data_dir,
                           "MS2_no_MS1_peak_match.csv"))) {
  write.table(data.frame("sample_code", "no_match", "mz", "rt", "int", "blank"), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir,
                                "MS2_no_MS1_peak_match.csv"),
              row.names = FALSE, col.names = FALSE)
}
total_spectra <- list()
blank_mzs <- c() # store the MS2 mzs that appear in blank samples
sample_count <- 0
total_ions_processed <- 0
ti <- Sys.time()
error_reading <- c()

# store the number of mzs that had x MS1 peaks extracted for them, where
# x from 0 to 300 peaks or more (aggregated in the 301 position)
# this indicates possible number of isomers (in each index) by m/z (count)
mzs_number_peaks <- numeric(length = 301) 

##### Match MS1 peak info with tandem MS2 data and 
##### rewrite raw data to mgf with peak annotation

# this is the strat of the pre process for each sample
t1 <- Sys.time()
cat("\n*** Starting processing samples ***\n\n")

# get the not blank samples indexes to keep for alignment - misalignment computation
# max to n_samples_batch_align by batch

if (n_samples_batch_align == 0) # do not align
{
  sample_idxs_keep <- c()
} else {
  sample_idxs_keep <- unlist(sapply(unique(metadata$DATA_COLLECTION_BATCH), function(batch) {
    idxs_keep <- which(metadata$DATA_COLLECTION_BATCH == batch & 
            metadata$SAMPLE_TYPE != "blank")
    # if more not blank samples in the batch then the max to keep, 
    # takes a sample from it with the max allowed size containing the first not blank sample
    if (length(idxs_keep) > n_samples_batch_align)
    {
      # select the samples in sequence
      idxs_keep <- idxs_keep[1:n_samples_batch_align] 
      #sample(sample_idxs_keep[-1], n_samples_batch_align-1)) # uniformely select the samples
      # sample_idxs_keep <- sample_idxs_keep[1:n_samples_batch_align]
    }
    idxs_keep
  }))
}


# load all data into a msExperiment
ms_data <- readMsExperiment(spectraFiles = metadata$FILENAME, sampleData = metadata)

# extract peak info of all samples
raw_datas <- lapply(seq_len(nrow(metadata)), function(sample_idx)
{
  #gc()
  cat("\n* Extracting peak info of sample", 
          metadata$SAMPLE_CODE[[sample_idx]], "*\n")
  #sample_count <<- sample_count + 1
  ti <- Sys.time()
  
  #ms_data <- readMsExperiment(spectraFiles = metadata$FILENAME[[sample_idx]], sampleData = pd[sample_idx,])
  sample_code <- metadata$SAMPLE_CODE[[sample_idx]]
  # extract ms1 peaks using CentWave algorithm
  ms_peaks <- findChromPeaks(ms_data[sample_idx], param = cwp)
  
  # return all ms2 associated with a found peak, 
  # multiple MS2 ions are assigned to all MS1 peaks that encompasses it in RT -
  # this may create duplicated MS2 to different chrom_peak_id
  # column into is the Intensity (Original / Total) and intb is intensity (Baseline-corrected)
  # intb is better for quantification, and into for checking raw instrument response
  # maxo is the peak Height - The absolute intensity at the very top of the peak.
  # into is the peak Area - The raw total area under the peak curve (includes baseline noise).
  # intb is the peak Area - The true area under the peak curve (with baseline noise subtracted).
  # sn is the noise Ratio - How many times taller the peak apex is compared to the background noise.
  ms2_peaks <- chromPeakSpectra(ms_peaks, msLevel = 2L, expandRt = rt_tol,
                                skipFilled = TRUE,
                                expandMz = mz_tol, return.type = "Spectra",
                                chromPeakColumns=c("rt", "mz", "intb",
                                                   "rt", "rtmin", "rtmax"), 
                                method = "all")
  # remove duplicated MS2 scans - duplicated acquisitionNum (scans from raw data)
  # compute the error in rt within the MS1 match and keep the duplicated correspondence with smaller error
  # we need all unique MS2 here
  # compute rt rmse from the MS2 rt and the MS1 rt boundaries without rt_tol
  # and compute rt diff from the MS2 rt and the MS1 rtmean
  ms2_peaks$rt_match_rmse <- sqrt(((ms2_peaks$chrom_peak_rtmin + rt_tol - ms2_peaks$rtime)^2 + (ms2_peaks$chrom_peak_rtmax - rt_tol - ms2_peaks$rtime)^2) / 2)
  # if the original peak bondaries encompasses the MS2 ion, set its rmse to 0
  ms2_peaks$rt_match_rmse[(ms2_peaks$rtime >= ms2_peaks$chrom_peak_rtmin + rt_tol) & 
                            (ms2_peaks$rtime <= ms2_peaks$chrom_peak_rtmax - rt_tol)] <- 0
  ms2_peaks$rt_match_diff <- abs(ms2_peaks$chrom_peak_rt - ms2_peaks$rtime)
  ms2_peaks$rt_match_error <- ms2_peaks$rt_match_rmse + ms2_peaks$rt_match_diff
  # remove duplicates keeping the one with smaller rt_match_error
  order_rt_match_error <- order(ms2_peaks$rt_match_error)
  unique_ms2_peak_by_error <- which(!duplicated(ms2_peaks$acquisitionNum[order_rt_match_error]))
  ms2_peaks <- ms2_peaks[order_rt_match_error[unique_ms2_peak_by_error]]
  # restore original MS1 peak boundaries, only enlarge it if the MS2 match used the rt_tol
  ms2_peaks$chrom_peak_rtmin <- pmin(ms2_peaks$chrom_peak_rtmin + rt_tol, ms2_peaks$rtime)
  ms2_peaks$chrom_peak_rtmax <- pmax(ms2_peaks$chrom_peak_rtmax - rt_tol, ms2_peaks$rtime)
  # get number of ms2 with ms1 peak match
  n_with_ms1 <- length(ms2_peaks)
  
  # for completude, rescue the MS2 without a MS1 peak match
  # Grab EVERY MS2 spectrum from your raw data
  ms2_all <- spectra(filterMsLevel(ms_data[sample_idx], msLevel = 2))
  n <- length(ms2_all)
  cat("Total number of MS2 ions = ",n, "\n") # LET = 4312
  cat("Number of MS2 with a MS1 peak correspondence = ",n_with_ms1, " (", 
      round(n_with_ms1/n*100,2), "%)\n") # 2793  ( 64.77 %)
  
  # filter only the MS2 not associated with a MS1 peak - fake peaks
  ms2_no_peaks <- ms2_all[!(ms2_all$acquisitionNum %in% ms2_peaks$acquisitionNum)]
  # create a fake MS1 peak for the MS2 without a MS1 match
  ms2_no_peaks$chrom_peak_rtmin <- ms2_no_peaks$rtime - rt_tol
  ms2_no_peaks$chrom_peak_rtmax <- ms2_no_peaks$rtime + rt_tol
  ms2_no_peaks$chrom_peak_intb <- ms2_no_peaks$precursorIntensity
  n_no_ms1 <- length(ms2_no_peaks) 
  cat("Number of MS2 without a MS1 correspondence = ",n_no_ms1, " (", 
      round(n_no_ms1/n*100,2), "%)\n") # 1519  ( 35.23 %)
  
  # trunc peak areas
  ms2_peaks$chrom_peak_intb <- trunc(ms2_peaks$chrom_peak_intb)
  # set peak ID as the chromatographic id + _ + sample_code + _ + peak_area
  ms2_peaks$peakID <- paste0(sapply(ms2_peaks$chrom_peak_id, 
                                        function(x) as.numeric(substring(x,3)),
                                        USE.NAMES = FALSE),
                                 "_",sample_code, "_", ms2_peaks$chrom_peak_intb)
  # for ms2 with no ms1 set as fake_ + sample_code + _ + peak_area
  ms2_no_peaks$peakID <- paste0("fake_",sample_code, "_", ms2_no_peaks$chrom_peak_intb)
  # concatenate all ms2 ions
  ms2_peaks_all <- c(ms2_peaks,ms2_no_peaks)
  # sort by acquisitionNum
  ms2_peaks_all <- ms2_peaks_all[order(ms2_peaks_all$acquisitionNum)]
  n_processed <- length(ms2_all)
  # check if length is kept here, for completude all ms2 ions must be here
  if (n_processed != n) {
    stop("Error: A different number of MS2 ions was obtained after removing ",
         "duplicated peaks from MS1 match and creating fake peaks for MS2 without MS1.",
         "Total number of raw MS2 = ", n, ". Total number of pre processed MS2 = ", n_processed)
  }
  
  # make sure scanIndex matches acquisitionNum
  ms2_peaks_all$scanIndex <- ms2_peaks_all$acquisitionNum
  # set title
  ms2_peaks_all$title <- paste0(sample_code,".",ms2_peaks_all$scanIndex)
  # set charge
  ms2_peaks_all$charge <- ion_mode
  
  # replace any zero precursor intensity spectrum with the sum of ints
  precursorIntensity_zero <- which(ms2_peaks_all$precursorIntensity <= 0)
  if (length(precursorIntensity_zero) > 0) {
    ms2_peaks_all$precursorIntensity[precursorIntensity_zero] <- ionCount(ms2_peaks_all[precursorIntensity_zero])
  }
  
  # retrieve the list of MS1 peaks
  ms1_peaks <- chromPeaks(ms_peaks)
  m <- nrow(ms1_peaks)
  cat("Total number of MS1 peaks = ",m, "\n") # 2821
  # ms1_peaks[ms2_peaks$chrom_peak_id[1],]
  # retrieve list of ms1 peaks without a MS2 correspondence
  ms1_peak_id_with_ms2 <- unique(ms2_peaks$chrom_peak_id)
  ms1_peak_id_no_ms2 <- rownames(ms1_peaks)[!(rownames(ms1_peaks) %in% ms1_peak_id_with_ms2)]
  # remove joined ms2 lists
  #rm(ms2_peaks,ms2_no_peaks)
  
  ms1_peaks_no_ms2 <- ms1_peaks[ms1_peak_id_no_ms2,]
  m_no_ms2 <- nrow(ms1_peaks_no_ms2)
  cat("Number of MS1 without a MS2 correspondence = ",m_no_ms2, " (", 
      round(m_no_ms2/m*100,2), "%)\n") # 2187  ( 77.53 %)
  
  ms1_peaks_with_ms2 <- ms1_peaks[ms1_peak_id_with_ms2,]
  m_with_ms2 <- nrow(ms1_peaks_with_ms2)
  cat("Number of MS1 with a MS2 correspondence = ",m_with_ms2, " (", 
      round(m_with_ms2/m*100,2), "%)\n") # 634  ( 22.47 %)
  
  # TODO filter ms1_peaks_with_ms2 using sn? remove no matched peaks with small signal to noise
  #nrow(ms1_peaks_no_ms2[ms1_peaks_no_ms2[,"sn"] >= max(snthresh,5),])
  
  
  # Create a custom field map
  custom_map <- c(scanIndex="SCANS",
                  title="TITLE",
                  chrom_peak_rtmin="RTMIN",
                  chrom_peak_rtmax="RTMAX",
                  peakID="PEAK_ID",
                  chrom_peak_intb="PEAK_AREA",
                  charge="CHARGE",
                  peaksCount="NUM_PEAKS",
                  rtime="RTINSECONDS",
                  precursorIntensity="PRECURSOR_INTENSITY",
                  precursorMz="PEPMASS")  #spectraVariableMapping(MsBackendMgf()))
  
  # Export with custom mappings
  file_MGF <- file.path(path_raw_data, processed_data_dir, 
                        paste0(metadata$SAMPLE_CODE[[sample_idx]], 
                               "_peak_info.mgf"))
  export(selectSpectraVariables(ms2_peaks_all, names(custom_map)), 
         backend = MsBackendMgf(), 
         file = file_MGF, 
         mapping = custom_map)
  tend <- Sys.time()
  cat("    * Done in", round(tend-t1, 2), units(tend-t1), "*\n")
  
  # TODO just check if verbose >= 10 - debug
  check_mgf(file_MGF, tags=custom_map, scans_list=ms2_peaks_all$scanIndex)
  
  
  ## TODO continue here check whats missing in the new code
  # TODO add blank baseline
  # TODO check if any precint is zero ~ raw_int[raw_int == 0] <- 1 # prevent zero count when spectra is present
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # read the raw data MS1
  raw_data <- tryCatch(readMSData(metadata$FILENAME[[sample_idx]], 
                                  mode = "onDisk", msLevel. = 1L, 
                                  pdata = new("NAnnotatedDataFrame", pd[sample_idx,])), 
                       error = function(e) message(e, " Error readMSData(msLevel = 1) from ",
                                                   metadata$FILENAME[[sample_idx]]))
  
  ## Perform the peak detection using the settings defined above
  # and get the detected MS1 peaks of each sample in a MSnExp object
  peakInfo <- tryCatch(findChromPeaks(raw_data, param = cwp), 
                       error = function(e) message(e, " Error findChromPeaks()"))
  # write.csv(chromPeaks(peakInfo, ppm = ppm_tol), 
  #           file =  file.path(path_raw_data, processed_data_dir, "ms1PeakInfo.csv"),
  #           row.names = FALSE)
  rm(raw_data)
  # read the tandem data MS2
  raw_data <- tryCatch(readMSData(file.path(metadata$FILENAME[[sample_idx]]), 
                                  msLevel. = 2L, verbose = 1), 
                       error = function(e) {
                         message(e, " Error reading MS2 data from ", 
                                                   file.path(metadata$FILENAME[[sample_idx]]))
                         NULL
                         })
  if (is.null(raw_data) || is.null(peakInfo)) {
    error_reading <- c(error_reading, metadata$FILENAME[[sample_idx]])
    return(NULL)
  }
  raw_prec <- precursorMz(raw_data)
  raw_rt   <- rtime(raw_data)
  raw_int <- ionCount(raw_data) # precursor intensity
  raw_int[raw_int == 0] <- 1 # prevent zero count when spectra is present
  
  # set min MS1 mz with the min ms2 mz minus a gap of 25 Daltons
  ms1_min_mz <- min(raw_prec) - 25
  # set the max MS1 mz with the max ms2 mz times 3 (triple charge variants) plus a gap
  ms1_max_mz <- max(raw_prec) * 3 + 25
  
  # increment number of processed ions
  total_ions_processed <<- total_ions_processed + length(raw_prec)
  
  # get peaks width info for the MS2 mz range 
  # sample_chromatogram <- tryCatch(as.data.frame(chromPeaks(peakInfo, ppm = ppm_tol)) %>%
  #                                   filter(mzmin >= min(raw_prec) - mz_tol, 
  #                                          mzmax <= max(raw_prec) + mz_tol) %>%
  #                                   arrange(mz), 
  #                                 error = function(e) message(e, " Error chromPeaks()"))
  sample_chromatogram <- tryCatch(as.data.frame(chromPeaks(peakInfo, ppm = ppm_tol)) %>%
                                    filter(mzmin >= ms1_min_mz, 
                                           mzmax <= ms1_max_mz) %>%
                                    arrange(mz), 
                                  error = function(e) message(e, " Error chromPeaks()"))
  sample_chromatogram$ms2_matches <- 0
  sample_chromatogram$intb[sample_chromatogram$intb == 0] <- 1
  sample_chromatogram$id <- paste0(row.names(sample_chromatogram), "_", metadata$SAMPLE_CODE[[sample_idx]])
  #group peaks jagged
  # gsample <- groupChromPeaks(peakInfo, PeakDensityParam(sampleGroups = c(1),
  #                                                       minFraction = min_fraction, binSize = bin_size,
  #                                                       bw = 2, maxFeatures = max_features))
  # sample_chrom_grouped <- tryCatch(as.data.frame(featureDefinitions(gsample)) %>%
  #                                    filter(mzmin >= min(raw_prec) - mz_tol,
  #                                           mzmax <= max(raw_prec) + mz_tol,
  #                                           npeaks > 1) %>%
  #                                    arrange(mzmed),
  #                                  error = function(e) message(e, " Error chromPeaks()"))
  fake_id <- paste0("fake_", metadata$SAMPLE_CODE[[sample_idx]])
  no_mz <- c()  # store the index of the ms2 spectra without a ms1 mz match
  no_rt <- c()  # store the index of the ms2 spectra without a ms1 rt match
  tandem_peak_info <- bind_rows(lapply(seq_along(raw_prec), function(i)
  {
    # try to find the mass and select the less distant in rt
    matches_mzs <- sample_chromatogram[raw_prec[[i]] >= sample_chromatogram$mzmin - mz_tol & 
                                         raw_prec[[i]] <= sample_chromatogram$mzmax + mz_tol, 
                                       c("id", "mz", "rt", "rtmin", "rtmax", "intb")]
    
    # return match
    if (nrow(matches_mzs) == 0) { # no valid mz range found, return fake MS2 peak info
      # message("No mz ", raw_prec[[i]], " rt ", raw_rt[[i]])
      if (is.null(raw_int[[i]]) || raw_int[[i]] > 1)
        no_mz <<- c(no_mz, i)
      
      # no_mzs <<- c(no_mzs,raw_prec[[i]])
      if (metadata$SAMPLE_TYPE[[sample_idx]] != "blank") {
        list(id=fake_id, mz=raw_prec[[i]],rt=raw_rt[[i]],rtmin=max(raw_rt[[i]]-mean(peak_width)/2, 0),
             rtmax=raw_rt[[i]]+mean(peak_width)/2, intb = raw_int[[i]])
      } else {
        list(id=fake_id, mz=raw_prec[[i]],rt=raw_rt[[i]],rtmin=0,
             rtmax=1000000, intb = raw_int[[i]])
      }
    } else {
      #matches_mzss <- matches_mzs
      matches_mzs <- matches_mzs[raw_rt[[i]] >= matches_mzs$rtmin - rt_tol & 
                                   raw_rt[[i]] <= matches_mzs$rtmax + rt_tol, ]
      
      if (nrow(matches_mzs) == 0) { # no valid rt range found, return fake MS2 peak info +- rt_tol
        if (is.null(raw_int[[i]]) || raw_int[[i]] > 1)
          no_rt <<- c(no_rt, i)
        # message("mz ", raw_prec[[i]], " rt ", raw_rt[[i]])
        # print(matches_mzss)
        if (metadata$SAMPLE_TYPE[[sample_idx]] != "blank") {
          list(id=fake_id, mz=raw_prec[[i]],rt=raw_rt[[i]],rtmin=max(raw_rt[[i]]-mean(peak_width)/2, 0),
               rtmax=raw_rt[[i]]+mean(peak_width)/2, intb = raw_int[[i]])
        } else {
          list(id=fake_id, mz=raw_prec[[i]],rt=raw_rt[[i]],rtmin=0,
               rtmax=1000000, intb = raw_int[[i]])
        }
      } else if (nrow(matches_mzs) > 1) { # if tied in rt range, get the closest in rt 
        # select best match
        # matches_mzs$diff_mz <- abs(matches_mzs$mz -  raw_prec[[i]])
        # matches_mzs <- top_n(matches_mzs, desc(1), diff_mz)
        # matches_mzs$diff_mz <- NULL
        #print(i)
        matches_mzs <- matches_mzs[order(abs(matches_mzs$rt -  raw_rt[[i]]), 
                                         decreasing = FALSE, na.last = TRUE), ]
        j <- which(matches_mzs[1, "id"] == sample_chromatogram$id)
        sample_chromatogram$ms2_matches[j] <<- sample_chromatogram$ms2_matches[j] + 1
        if (metadata$SAMPLE_TYPE[[sample_idx]] != "blank") {
          matches_mzs[1,]
        } else {
          matches_mzs[1,4] <- 0
          matches_mzs[1,5] <- 1000000
          matches_mzs[1,]
        }
      } else { # only one match in mz and rt ranges
        #print(i)
        j <- which(matches_mzs$id == sample_chromatogram$id)
        sample_chromatogram$ms2_matches[j] <<- sample_chromatogram$ms2_matches[j] + 1
        if (metadata$SAMPLE_TYPE[[sample_idx]] != "blank") {
          matches_mzs
        } else {
          matches_mzs[1,4] <- 0
          matches_mzs[1,5] <- 1000000
          matches_mzs
        }
      }
    }
  }))
  #save not ms2_matches peaks
  sample_chromatogram <- filter(sample_chromatogram, 
                                mzmin >= ms1_min_mz - mz_tol, 
                                mzmax <= ms1_max_mz + mz_tol)
  sample_chromatogram$intb <- trunc(as.numeric(sample_chromatogram$intb))
  sample_chromatogram <- sample_chromatogram[, c("id", "mz", "rt", "rtmin", 
                                                 "rtmax", "ms2_matches", "intb")]
  # set sample peak area column name
  names(sample_chromatogram) <- c("peakIds", "mz", "rtMean", "rtMin", 
                                  "rtMax", "ms2_matches",
                                  paste0(metadata$SAMPLE_CODE[[sample_idx]], "_area"))
  # add count columns of all samples
  sample_chromatogram[,no_ms2_htable[(!(no_ms2_htable %in% names(sample_chromatogram)))]] <- 0
  # order columns
  sample_chromatogram <- sample_chromatogram[,no_ms2_htable]
  # round numeric columns
  sample_chromatogram$mz <- round(as.numeric(sample_chromatogram$mz) , 4)
  sample_chromatogram$rtMean <- round(as.numeric(sample_chromatogram$rtMean) , 4)
  sample_chromatogram$rtMin <- round(as.numeric(sample_chromatogram$rtMin) , 4)
  sample_chromatogram$rtMax <- round(as.numeric(sample_chromatogram$rtMax) , 4)
  
  if (metadata$SAMPLE_TYPE[[sample_idx]] == "blank") {
    sample_chromatogram$rtMin <- 0
    sample_chromatogram$rtMax <- 1000000
    if (!all(tandem_peak_info$rtmin == 0) || !all(tandem_peak_info$rtmax == 1000000))
      stop("Error BASELINE RT RANGE for blank sample ", metadata$SAMPLE_CODE[[sample_idx]])
  } else {
    # compute number of isomers (MS1 peaks by m/z)
    mzs_number_peaks <<- mzs_number_peaks + 
      number_mzs_by_number_peaks(sample_chromatogram[,c('peakIds', 'mz')],
                                 blank_mzs, mz_tol)
  }
  # print(number_mzs_by_number_peaks(sample_chromatogram[,c('peakIds', 'mz')], 
  #                                  blank_mzs, mz_tol))
  # save the matched and not matches MS1 mzs
  write.table(filter(sample_chromatogram, ms2_matches > 0), 
              sep = ",",
              file.path(path_raw_data, processed_data_dir,
                        "MS1_list_with_MS2.csv"),
              row.names = FALSE, col.names = FALSE,
              append = TRUE)
  sample_chromatogram <- filter(sample_chromatogram, ms2_matches == 0)
  sample_chromatogram$peakIds <- paste0(sample_chromatogram$peakIds, "_NOPEAKS")
  write.table(sample_chromatogram, 
              sep = ",",
              file.path(path_raw_data, processed_data_dir,
                        "MS1_list_no_MS2.csv"),
              row.names = FALSE, col.names = FALSE,
              append = TRUE)
  rm(sample_chromatogram)

  # write MS2 data with peak info and a unique groupID equals the scan number
  # using intb -  Integrated baseline-corrected intensity - which is the $into - original integrated peak intensity -
  # minus local noise background - corrected peak area
  tryCatch(writeMgfDataFile_NP3(raw_data, 
                                file_MGF = file.path(path_raw_data, processed_data_dir, 
                                                paste0(metadata$SAMPLE_CODE[[sample_idx]], 
                                                       "_peak_info.mgf")), 
                                RTMIN = tandem_peak_info$rtmin,
                                RTMAX = tandem_peak_info$rtmax,
                                INTO = raw_int,
                                TITLE = metadata$SAMPLE_CODE[[sample_idx]],
                                PEAK_ID = tandem_peak_info$id,
                                PEAK_AREA = trunc(tandem_peak_info$intb),
                                CHARGE = ion_mode,
                                tags = c("BEGIN IONS", "SCANS", "RTINSECONDS",
                                         "RTMIN", "RTMAX", "PEAK_AREA", "PEAK_ID",
                                         "CHARGE", "PRECURSOR_INTENSITY", "PEPMASS",
                                         "END IONS")), 
           error = function(e) stop(e, "Error writting MGF."))

  # set number of MS2 spectra in the current sample
  total_spectra[metadata$SAMPLE_CODE[[sample_idx]]] <<- length(tandem_peak_info$id)
  
  rm(raw_data, tandem_peak_info)
  tf <- Sys.time()
  
  cat("Number of MS2 spectra without a MS1 peak m/z range correspondence:", 
      length(no_mz), "\n")
  cat("Number of MS2 spectra without a MS1 peak retention time width correspondence:", 
      length(no_rt), "\n")
  cat("DONE for sample", metadata$SAMPLE_CODE[[sample_idx]],
          "in", round(tf-ti, 2), units(tf-ti), "!\n")
  cat("(", sample_count, "/", nrow(metadata), ")\n", sep = "")
  
  # store the MS2 mzs that did not have a MS1 peak match
  write.table(data.frame(sample = rep(metadata$SAMPLE_CODE[[sample_idx]],
                                      length(no_rt)+length(no_mz)),
                         no_match = c(rep("no_mz", length(no_mz)),
                                      rep("no_rt", length(no_rt))),
                         mz = raw_prec[c(no_mz, no_rt)],
                         rt = raw_rt[c(no_mz, no_rt)], 
                         int = raw_int[c(no_mz, no_rt)], 
                         blank = rep(ifelse(metadata$SAMPLE_TYPE[[sample_idx]] == "blank",
                                            1, 0),
                                     length(no_rt)+length(no_mz))), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir,
                                "log_MS2_no_MS1peak_match.csv"),
              row.names = FALSE, append = TRUE, col.names = FALSE)
  
  if (metadata$SAMPLE_TYPE[[sample_idx]] == "blank") {
    # blank_mzs <<- unique(c(blank_mzs, round(raw_prec, decimalplaces(mz_tol)+2)))
    blank_mzs <<- almost.unique(c(blank_mzs, 
                                  almost.unique(raw_prec, mz_tol)), 
                                mz_tol)
    names(blank_mzs) <- NULL
  }
  
  # write.table(paste("finished sample idx", metadata$SAMPLE_CODE[[sample_idx]]), "memory_use", append = T)
  # write.table(gc(), "memory_use", append = T)
  gc()
  
  if (sample_idx %in% sample_idxs_keep)
    return(peakInfo)
  else
    return(NULL)
})
  
t2 <- Sys.time()
cat("\n** Done extracting peak info of samples in", 
        round(t2-t1, 2), units(t2-t1), "**\n")

# TODO update how to compute misalignment here
# get indexes of not blank samples in the batch
sample_idxs_keep <- match(sample_idxs_keep, 
                          which(metadata$DATA_COLLECTION_BATCH == batch))

# if just one not blank sample, do not need alignment
if (length(sample_idxs_keep) == 0) # do not align
{
  return(NULL)
} else if (length(sample_idxs_keep) == 1) # do not need alignment of only one sample
{
  write.table(data.frame(batch, NA,NA,NA,NA,NA,NA,NA), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir, 
                                "samples_alignment.csv"),
              row.names = FALSE, col.names = FALSE, append = TRUE)
  # batches_rttol <<- c(batches_rttol, 0)
  return(raw_datas[[sample_idxs_keep]])
}

# write.table(paste("start sligning batch ", batch), "memory_use", append = T)
# write.table(gc(), "memory_use", append = T)
gc()

cat("\n*** Starting samples alignment of data collection batch", batch, "***\n\n")
# if more than one not blank sample by batch, make alignment
pdp <- PeakDensityParam(sampleGroups = rep(batch, length(sample_idxs_keep)),
                        minFraction = min_fraction, binSize = bin_size, 
                        bw = bw_win, maxFeatures = max_features)
xdata <- tryCatch(groupChromPeaks(do.call(concatenate_XCMSnExp, 
                                          raw_datas[sample_idxs_keep]), 
                                  param = pdp), 
                  error = function(e) {
                    message(e, ". Alignment disabled for batch ", batch, 
                            ". Error in groupChromPeaks()")
                    write.table(data.frame(batch, "Error"), sep = ",",
                                file =  file.path(path_raw_data, processed_data_dir, 
                                                  "samples_alignment.csv"),
                                row.names = FALSE, col.names = FALSE, append = TRUE)
                    return(NULL)
                    })
rm(raw_datas)

if (!is.null(xdata))  # alignment enabled
{
  # RT correction - align
  res_adjustRt <- tryCatch(adjustRtime(xdata, param = pgp), 
                    error = function(e) {
                      message(e, ". Alignment disabled for batch ", batch, 
                              ". Error in adjustRtime())")
                      write.table(data.frame(batch, "Error"), sep = ",",
                                  file =  file.path(path_raw_data, processed_data_dir, 
                                                    "samples_alignment.csv"),
                                  row.names = FALSE, col.names = FALSE, append = TRUE)
                      return(NULL)
                      })
  if (!is.null(res_adjustRt)) # alignment disable
  {
    write.table(t(data.frame(c(batch, 
                               paste(metadata$SAMPLE_CODE[sample_idxs_keep], collapse = ','),
                               summary(abs(rtime(res_adjustRt, adjusted = FALSE) - 
                                                rtime(res_adjustRt, adjusted = TRUE)))))), sep = ",",
                file =  file.path(path_raw_data, processed_data_dir, "samples_alignment.csv"),
                row.names = FALSE, col.names = FALSE, append = TRUE)
    # plot alignment to compare diffs
    # plotAdjustedRtime(res_adjustRt, col = rainbow(length(sample_idxs_keep)),
    #                   peakGroupsCol = "grey", peakGroupsPch = 1)
    # legend(x="bottomright", legend=res_adjustRt$sample_name, 
    #        fill = rainbow(length(sample_idxs_keep)), cex = 0.7)
    #plot(chromatogram(res_adjustRt))
  }
}

t1 <- Sys.time()
cat("DONE! Time elapsed", round(t1-t2, 2), units(t1-t2), "\n")
# return the data of one sample representant by batch, the first not blank sample
xdata1 <- tryCatch(filterFile(xdata, 1), 
                   error = function(e) message(e, " Error filterFile(xdata, 1))")) 
rm(xdata)
gc()
return(xdata1)

tf <- Sys.time()


if (total_ions_processed == 0)
{
  stop("\nERROR pre processing the samples. No MS2 spectra was found!\n",
       "Please check if your raw data contains both MS1 and MS2 spectra, ",
       "which are mandatory for the NP3 MS workflow pipeline.\n")
}

# remove NULL
batches_data <- batches_data[!sapply(batches_data, is.null)]
# check if should align
if (is.null(unlist(batches_data))) # no alignment
{
  cat("\nDONE pre processing all samples!\nTotal number of MS2 spectra processed:",
          total_ions_processed, 
          "\nTime elapsed", 
          round(tf-ti, 2), units(tf-ti), "\n")
  
  check_preprocess_correspondence(path_raw_data, processed_data_dir, 
                                  total_spectra, blank_mzs, mzs_number_peaks,
                                  mz_tol)
  
  q("no")
} else {
  cat("\nDONE pre processing and aligning all samples by data collection batch!\n",
          "Total number of MS2 spectra processed:",total_ions_processed, 
          "\nTime elapsed", 
          round(tf-ti, 2), units(tf-ti), "\n")
}
# names(batches_rttol) <- unique(metadata$DATA_COLLECTION_BATCH)

# align all data for the integration step
if (length(batches_data) > 1)
{
  cat("\n*** Starting alignment of all", data_name, "samples ***\n\n")
  # align all not blank data
  pdp <- PeakDensityParam(sampleGroups = unique(metadata[metadata$SAMPLE_TYPE!='blank',
                                                         'DATA_COLLECTION_BATCH']),
                          minFraction = min_fraction, binSize = bin_size, bw = bw_win, 
                          maxFeatures = max_features*length(batches_data)/2)
  xdata <- tryCatch(groupChromPeaks(do.call(concatenate_XCMSnExp, 
                                            batches_data), 
                                    param = pdp), 
                    error = function(e) {
                      message(e, ". Alignment disabled for all. Error in groupChromPeaks()")
                      return(NULL)
                      })
  if (is.null(xdata)) # alignment disable
  {  
    cat("Alignment of all samples disabled. \n")
  } else {
    # RT correction
    xdata <- tryCatch(adjustRtime(xdata, param = pgp), 
                      error = function(e) {
                        message(e,  ". Alignment disabled for all. Error in adjustRtime()")
                        return(NULL)
                        })
    if (is.null(xdata)) # alignment disable
    {  
      cat("Alignment of all samples disabled.\n")
    } else {
    # plot alignment to compare diffs
    # plotAdjustedRtime(xdata, col = rainbow(length(xdata$sample_name)),
    #                   peakGroupsCol = "grey", peakGroupsPch = 1)
    # legend(x="bottomright", legend=xdata$sample_name, 
    #        fill = rainbow(length(xdata$sample_name)), cex = 0.7)
    # write the alignment of all samples
    write.table(t(data.frame(c("all", 
                               paste(sapply(unique(metadata[metadata$SAMPLE_TYPE!='blank',
                                                            'DATA_COLLECTION_BATCH']), 
                                            function (x) metadata[metadata$DATA_COLLECTION_BATCH == x,'SAMPLE_CODE'][[1]]),
                                     collapse = ','),
                               summary(abs(rtime(xdata, adjusted = FALSE) - 
                                            rtime(xdata, adjusted = TRUE)))))), sep = ",",
                file =  file.path(path_raw_data, processed_data_dir, "samples_alignment.csv"),
                row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
  # 
  # batches_rttol <- c(batches_rttol, 
  #                    all = max(abs(rtime(xdata, adjusted = FALSE) - 
  #                                    rtime(xdata, adjusted = TRUE))))
} else {
  # not enought sample to align all
  write.table(data.frame("all", NA, NA, NA, NA, NA, NA, NA), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir, 
                                "samples_alignment.csv"),
              row.names = FALSE, col.names = FALSE, append = TRUE)
  # batches_rttol <- c(batches_rttol, all = 0)
}

t3 <- Sys.time()
cat("DONE! Time elapsed", round(t3-tf, 2), units(t3-tf), "\n")

# write.csv(data.frame(batch = names(batches_rttol), 
#                      rtTol = unlist(batches_rttol, use.names = F)),
#           file =  file.path(path_raw_data, processed_data_dir, "samples_alignment.csv"),
#           row.names = FALSE)
cat("\nExtracted peak info and aligned samples of", data_name, "in", 
        round(t3-ti, 2), units(t3-ti), "\n")

check_preprocess_correspondence(path_raw_data, processed_data_dir, 
                                total_spectra, blank_mzs, mzs_number_peaks,
                                mz_tol)

if (length(error_reading) > 0)
{
  message("\nError reading the following mzXML files. The respectiveful MGFs were not created:\n",
          "- ", paste0(error_reading, collapse = ";"))
}
