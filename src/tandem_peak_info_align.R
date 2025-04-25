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
source(file.path(script_path(),"read_metadata_table.R"))
source(file.path(script_path(),"writeMgfData_NP3.R"))

decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

# ms2_id_list is expected to have id (peakID) and mz columns -> tandem_peak_info[,c('id', 'mz')]
number_mzs_by_number_peaks <- function(ms2_id_list, blank_mzs, mz_tol) {
  max_num_peaks <- 300 # maximum number of ms1 peaks to be counted
  # compute unique m/zs
  mzs_unique <- sort(almost.unique(ms2_id_list$mz, mz_tol))
  # remove mzs from blank sample
  for (mz_blank in blank_mzs) {
    mzs_unique <- mzs_unique[!(mzs_unique >= mz_blank - mz_tol &
                               mzs_unique <= mz_blank + mz_tol)]
  }
  # count number of peaks by m/z
  num_ms1_peaks <- sapply(mzs_unique, function (x) {
    length(unique(ms2_id_list[ms2_id_list$mz >= x - mz_tol & 
                                ms2_id_list$mz <= x + mz_tol, 'peakIds']))
    })
  # tabule number of m/z with a given number of peaks; aggregate into one 
  # the m/zs with a number of MS1 peaks above the maximum expected
  table_num_mz_by_num_peaks <- c(tabulate(num_ms1_peaks, nbins = max_num_peaks),
                                 sum(num_ms1_peaks > max_num_peaks))
  return(table_num_mz_by_num_peaks)
}


check_preprocess_correspondence <- function(path_raw_data, processed_data_dir, 
                                            total_spectra, blank_mzs,
                                            num_mzs_number_peak, mz_tol)
{
  # count number of samples using the total_spectra vector, more reliable
  n_samples <- length(total_spectra)
  #print(blank_mzs)
  # check if the pre process was successful in all correspondences, 
  # and if not send a warning
  mzs_nomatch <- suppressMessages(readr::read_csv(file.path(path_raw_data, processed_data_dir,
                                           "log_MS2_no_MS1peak_match.csv"),guess_max = 5000))
  cat("\n\n* Pre-processing - Statistics of the percentage of MS2 spectra without a MS1 peak correspondence *\n\n")
  if (nrow(mzs_nomatch) == 0) {
    cat("All MS2 spectra had a MS1 peak correspondence!!")
    return()
  }
  
  # add blank information to the list of ms1 peaks with a ms2 correspondence
  # and to the log of ms2 without ms1
  ms1_w_ms2_list <- suppressMessages(readr::read_csv(file.path(path_raw_data, processed_data_dir,
                                                               "MS1_list_with_MS2.csv"), guess_max = 5000))
  ms1_w_ms2_list$BFLAG <- FALSE
  if (length(blank_mzs) > 0) 
  {
    # remove mzs from blank sample
    for (mz_blank in blank_mzs) {
      mzs_nomatch[mzs_nomatch$blank == 0 & 
                    mzs_nomatch$mz >= mz_blank - mz_tol & 
                    mzs_nomatch$mz <= mz_blank + mz_tol, "blank"] <- 2
      ms1_w_ms2_list[!ms1_w_ms2_list$BFLAG &
                       ms1_w_ms2_list$mz >= mz_blank - mz_tol & 
                       ms1_w_ms2_list$mz <= mz_blank + mz_tol, "BFLAG"] <- TRUE
    }
  
    # update the log table with the mzs close to a blank mz, blank == 2
    readr::write_csv(mzs_nomatch, file.path(path_raw_data, processed_data_dir,
                                            "log_MS2_no_MS1peak_match.csv"))
  }
  # update the list of ms1 peaks with a ms2 correspondence with blank info
  readr::write_csv(ms1_w_ms2_list, file.path(path_raw_data, processed_data_dir,
                                          "MS1_list_with_MS2.csv"))
  remove(ms1_w_ms2_list)
  # filter only mzs from the blank samples and mzs different from the blank mzs in the other samples
  # in other words, remove no matchs from blank m/zs that appear in samples that are not blank
  summary_MS2_no_MS1 <- as.data.frame.matrix(table(mzs_nomatch[mzs_nomatch$blank != 2,
                                                               c("sample_code", "no_match")]))
  summary_MS2_no_MS1$total_MS2 <- sapply(row.names(summary_MS2_no_MS1), function(x) total_spectra[[x]])
  # if there is no mzs without a mz or a rt match, set the proper column to zero
  if (is.null(summary_MS2_no_MS1$no_mz)) {
    summary_MS2_no_MS1$no_mz <- 0
  }
  if (is.null(summary_MS2_no_MS1$no_rt)) {
    summary_MS2_no_MS1$no_rt <- 0
  }
  summary_MS2_no_MS1['rate_MS2_no_MS1_mz'] <- round(summary_MS2_no_MS1$no_mz/summary_MS2_no_MS1$total_MS2, 3) * 100
  summary_MS2_no_MS1['rate_MS2_no_MS1_rt'] <- round(summary_MS2_no_MS1$no_rt/summary_MS2_no_MS1$total_MS2, 3) * 100
  names(summary_MS2_no_MS1)[1:2] <- c("MS2_no_MS1_mz", "MS2_no_MS1_rt")
  summary_MS2_no_MS1 <- tibble::rownames_to_column(summary_MS2_no_MS1, var='SAMPLE_CODE')
  
  # check if the summary includes all sample codes, 
  # if not include the missing ones without no correspondece (good pre-processed samples)
  if (!all(names(total_spectra) %in% summary_MS2_no_MS1$SAMPLE_CODE)) {
    summary_MS2_no_MS1 <- bind_rows(summary_MS2_no_MS1,
              lapply(names(total_spectra)[!(names(total_spectra) %in% summary_MS2_no_MS1$SAMPLE_CODE)], 
           function(x) list(SAMPLE_CODE = x, MS2_no_MS1_mz = 0, MS2_no_MS1_rt= 0,
                            total_MS2 = 0, rate_MS2_no_MS1_mz=0.0, rate_MS2_no_MS1_rt=0.0)))
  }
  
  # print statistics of blank and not blank samples with the rates of no correspondence by sample
  summary_MS2_no_MS1 <- arrange(summary_MS2_no_MS1, rate_MS2_no_MS1_rt, rate_MS2_no_MS1_mz)
  print(summary_MS2_no_MS1[,c('SAMPLE_CODE', 'rate_MS2_no_MS1_mz', 'rate_MS2_no_MS1_rt')], width=200)
  # print legend
  cat('\n==========\n',
      '\n- Variables Description\n',
      '  - rate_MS2_no_MS1_mz \t: Percentage of MS2 spectra without a MS1 m/z range correspondence (%)\n',
      '  - rate_MS2_no_MS1_rt \t: Percentage of MS2 spectra without a MS1 retention time range correspondence (%)\n',
      '\n==========\n\n')
  
  # Only proceed the checking if there is at least one not blank sample
  if (!any(mzs_nomatch$blank == 0)) {
    return()
  }
  # filter only not blank mzs
  mzs_nomatch <- mzs_nomatch[mzs_nomatch$blank == 0,1:5]
  mzs_nomatch <- arrange(mzs_nomatch, desc(int))
  summary_MS2_no_MS1 <- as.data.frame.matrix(table(mzs_nomatch[,c("sample_code", "no_match")]))
  summary_MS2_no_MS1$total_MS2 <- sapply(row.names(summary_MS2_no_MS1), function(x) total_spectra[[x]])
  # if there is no mzs without a mz or a rt match, set the proper column to zero
  if (is.null(summary_MS2_no_MS1$no_mz)) {
    summary_MS2_no_MS1$no_mz <- 0
  }
  if (is.null(summary_MS2_no_MS1$no_rt)) {
    summary_MS2_no_MS1$no_rt <- 0
  }
  summary_MS2_no_MS1['rate_MS2_no_MS1_mz'] <- round(summary_MS2_no_MS1$no_mz/summary_MS2_no_MS1$total_MS2, 3) * 100
  summary_MS2_no_MS1['rate_MS2_no_MS1_rt'] <- round(summary_MS2_no_MS1$no_rt/summary_MS2_no_MS1$total_MS2, 3) * 100
  names(summary_MS2_no_MS1)[1:2] <- c("MS2_no_MS1_mz", "MS2_no_MS1_rt")
  summary_MS2_no_MS1 <- tibble::rownames_to_column(summary_MS2_no_MS1, var='SAMPLE_CODE')
  
  # round the mzs using the number of decimal places of the tolerance + 1 to group close by mzs below the tolerance
  mzs_nomatch$mz <- round(mzs_nomatch$mz, decimalplaces(mz_tol)+1)
  # if any MS2_noMatchingMzMS1 > 0.05 warning mz_tol and/or ppm_tol is too small to allow the correspondence
  if (any(summary_MS2_no_MS1['rate_MS2_no_MS1_mz'] > 5))
  {
    cat("\n***********\n* WARNING *  There are", 
        sum(summary_MS2_no_MS1['rate_MS2_no_MS1_mz'] > 5),
        "not blank samples with more than 5% of the MS2 spectra without correspondence in the MS1 m/z range\n***********\n")
    samples_nomatch <- summary_MS2_no_MS1$SAMPLE_CODE[summary_MS2_no_MS1['rate_MS2_no_MS1_mz'] > 5]
    mzs_nomatch_by_sample <- data.frame(table(mzs_nomatch[
      mzs_nomatch$sample_code %in% samples_nomatch & mzs_nomatch$no_match == 'no_mz', 
      c("sample_code", "mz")]), stringsAsFactors = FALSE)
    mzs_nomatch_by_sample <- mzs_nomatch_by_sample[mzs_nomatch_by_sample$Freq > 0,]
    
    mzs_nomatch_by_sample$mz <- as.numeric(as.character(mzs_nomatch_by_sample$mz))
    mzs_nomatch_by_sample$sample_code <- as.character(mzs_nomatch_by_sample$sample_code)
    mzs_nomatch_by_sample <- left_join(mzs_nomatch_by_sample, 
                                       mzs_nomatch[!duplicated(mzs_nomatch[,c('sample_code','mz')]),
                                                   c("sample_code", "mz","rt","int")],
                                       by=c("sample_code", "mz"))
    
    mzs_nomatch_by_sample <-  as.data.frame(mzs_nomatch_by_sample %>% 
                                              arrange(desc(Freq), desc(int), desc(mz)) %>% 
                                              group_by(sample_code) %>% slice(1:5))
    # List of MS2 spectra m/zs without a MS1 m/z range match - top 5 of the samples with more than 5% of missing correspondence
    cat("\n--> Top 5 list of MS2 spectra m/zs without a MS1 m/z range match - by sample with more than 5% of missing correspondence\n\n")
    # SAMPLE X
    # mz | number of spectra without a MS1 m/z range correspondence
    names(mzs_nomatch_by_sample) <- c("SAMPLE_CODE", "MS2 precursor m/z", "num_MS2_no_MS1_mz", "rt", "prec_int")
    print(mzs_nomatch_by_sample)
    cat('\n==========\n',
        '\n- Variable Description\n',
        '  - num_MS2_no_MS1_mz \t: Number of MS2 spectra with the respective m/z without a MS1 m/z range correspondence\n',
        '  - rt \t\t\t\t: The retention time of the most intense MS2 spectra with the respective m/z\n',
        '  - prec_int \t\t\t: The the precursor intensity of the most intense MS2 spectra with the respective m/z\n',
        '\n==========\n\n')
    
    cat("\nWARNING: the pre_process command mz_tol OR ppm_tol parameters can be too small to allow the correspondence. These parameters can not be representative of the true values found in the data. \n")
    cat("-> We strongly recommend the user to visualize the chromatograms of some, if possible all) of the above m/zs to define better values for the ppm_tol parameter or better evaluate the problem (if it's really an issue). The instrument accuracy and the data quality could also be causing these high rate of missing correspondence.\n")
    cat("-> The NP3 MS workflow command 'chr' can be used to visualize the chromatograms.\n")
  }
  
  # if any MS2_noMatchingRtMS1 > 0.05 warning the rt_tol parameter value could be too small or the peak_width parameter defining the chromatographic  peak  width  range  (wmin, wmax in seconds) is probably not representative of the true values found in the data. We strongly recommend the user to visualize the chromatograms of some, if possible all) of the following m/zs to define better values for the peak_width parameter.
  if (any(summary_MS2_no_MS1['rate_MS2_no_MS1_rt'] > 5))
  {
    cat("\n***********\n* WARNING *  There are", 
        sum(summary_MS2_no_MS1['rate_MS2_no_MS1_rt'] > 5),
        "not blank samples with more than 5% of the MS2 spectra without correspondence in the MS1 RT range\n***********\n")
    samples_nomatch <- summary_MS2_no_MS1$SAMPLE_CODE[summary_MS2_no_MS1['rate_MS2_no_MS1_rt'] > 5]
    mzs_nomatch_by_sample <- data.frame(table(mzs_nomatch[
      mzs_nomatch$sample_code %in% samples_nomatch & mzs_nomatch$no_match == 'no_rt', 
      c("sample_code", "mz")]))
    mzs_nomatch_by_sample <- mzs_nomatch_by_sample[mzs_nomatch_by_sample$Freq > 0,]
    
    mzs_nomatch_by_sample$mz <- as.numeric(as.character(mzs_nomatch_by_sample$mz))
    mzs_nomatch_by_sample$sample_code <- as.character(mzs_nomatch_by_sample$sample_code)
    mzs_nomatch_by_sample <- left_join(mzs_nomatch_by_sample, 
                                       mzs_nomatch[!duplicated(mzs_nomatch[,c('sample_code','mz')]),
                                                   c("sample_code", "mz","rt","int")],
                                       by=c("sample_code", "mz"))
    
    mzs_nomatch_by_sample <-  as.data.frame(mzs_nomatch_by_sample %>% 
                                              arrange(desc(Freq), desc(int), desc(mz)) %>% 
      group_by(sample_code) %>% slice(1:5))
    # List of MS2 spectra m/zs without a MS1 retention time range match - top 5 of the samples with more than 5% of missing correspondence
    cat("\n--> Top 5 list of MS2 spectra m/zs without a MS1 retention time range match - by sample with more than 5% of missing correspondence\n\n")
    # SAMPLE X
    # mz | number of spectra without a MS1 peak width correspondence
    names(mzs_nomatch_by_sample) <- c("SAMPLE_CODE", "MS2 precursor m/z", "num_MS2_no_MS1_rt", "rt", "prec_int")
    print(mzs_nomatch_by_sample)
    cat('\n==========\n',
        '\n- Variable Description\n',
        '  - num_MS2_no_MS1_rt \t: Number of MS2 spectra with the respective m/z without a MS1 retention time range correspondence\n',
        '  - rt \t\t\t\t: The retention time of the most intense MS2 spectra with the respective m/z\n',
        '  - prec_int \t\t\t: The the precursor intensity of the most intense MS2 spectra with the respective m/z\n',
        '\n==========\n\n')
    
    
    cat("\nWARNING: the pre_process command rt_tol parameter value can be too small to allow finding the correspondences OR the peak_width parameter defining the chromatographic peak width range (wmin, wmax in seconds) can probably not be representative of the true values found in the data - this could cause the XCMS::CentWave algorithm to not work properly.\n")
    cat("-> We strongly recommend the user to visualize the chromatograms of some (if possible all) of the above MS2 spectra m/zs to define better values for the peak_width parameter or better evaluate the problem (if it's really an issue). The instrument accuracy and the data quality could also be causing these high rate of missing correspondence.\n")
    cat("-> The NP3 MS workflow command 'chr' can be used to visualize the chromatograms.\n")
  }
  
  cat('\n==========\n',
      'Summary of the rate_MS2_no_MS1_mz value for all not blank samples\n\n')
  print(summary(summary_MS2_no_MS1$rate_MS2_no_MS1_mz))
  cat('\n Number of not blank samples with the rate_MS2_no_MS1_mz value above 5%:',
      sum(summary_MS2_no_MS1['rate_MS2_no_MS1_mz'] > 5), '/',n_samples)
  cat('\n==========\n')
  cat('\n==========\n',
      'Summary of the rate_MS2_no_MS1_rt value for all not blank samples\n\n')
  print(summary(summary_MS2_no_MS1$rate_MS2_no_MS1_rt))
  cat('\n Number of not blank samples with the rate_MS2_no_MS1_rt value above 5%:',
      sum(summary_MS2_no_MS1['rate_MS2_no_MS1_rt'] > 5), '/',n_samples)
  cat('\n==========\n\n')
  
  # plot the number of mzs by the number of MS1 peaks (isomers)
  
  png(filename = file.path(path_raw_data, processed_data_dir,
                           "log_number_mzs_by_number_peaks.png"))
  par(mar = c(5.1,5.1,4.1,1.1))
  xlim_pos <- max(which(!(num_mzs_number_peak == 0)))+3
  barplot(Num_mzs ~ Num_peaks, data = cbind(Num_mzs = c(num_mzs_number_peak+1)[1:xlim_pos], 
                                            Num_peaks = 1:xlim_pos),
          xlab = "Number of MS1 Peaks", ylim=c(1,max(num_mzs_number_peak)),
          ylab = "Number of m/zs\nwith the given Number of MS1 Peaks", 
          ann = TRUE, axes=TRUE, axisnames = TRUE, log="y", 
          main="Number of m/zs in all Samples\nwith the given Number of MS1 peaks (Putative Isomers)")
  dev.off()
  num_mzs_geq30 <- sum(num_mzs_number_peak[30:length(num_mzs_number_peak)])
  if (num_mzs_geq30/sum(num_mzs_number_peak)*100 > 5)
  {
    cat("\n***********\n* WARNING *  There are", 
        num_mzs_geq30/sum(num_mzs_number_peak)*100,
        "% of m/zs in all samples with more than 30 MS1 peaks (possible isomers).",
        "\nThis can indicate a BAD integration that split real MS1 peaks into multiple smaller peaks.\n",
        "To prevent splitting real peaks, reset the 'peak_width' parameter according to what is found in your data.\n",
        "The full distribution of the number of m/zs by the number of MS1 peaks is present in the 'log_number_mzs_by_number_peaks.png' file.\n",
        "Check the matched MS1 peaks in the 'MS1_list_with_MS2.csv' file of the pre-processing result to see if this is really an issue.",
        "If you expect to have a big number of isomers in your samples, you can ignore this warning.",
        "\n***********\n\n")
  }
}

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
  snthresh  <- 0   #Signal to noise ratio threshold 
  prefilter <- c(1,750) # only keep mzs that have at least x peaks with intensity y : (x,y)
  min_fraction <- 0.3    # defining the minimum fraction of samples in at least one sample group in which the peaks have to be present to be considered as a peak group (feature).
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
no_ms2_htable <-  c("peakIds", "mz", "rtMean", "rtMin", 
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
  # if at least one sample was already preprocessed skip alignment
  if (sum(overwrite_preprocess) < nrow(metadata)) {
    cat("Using already pre processed files, only pre processing the missing ones. ",
            "Alignment disabled.\n")
    
    n_samples_batch_align <- 0
    # do not overwrite table header to store not fragmented mzs
  } else {
    # write table header to store not fragmented mzs
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

cat("Loading packages XCMS, MSnbase...\n")
suppressPackageStartupMessages(library(xcms))
suppressPackageStartupMessages(library(MSnbase))

## Use socket based parallel processing on Windows systems
# if (.Platform$OS.type == "unix") {
#     register(bpstart(MulticoreParam(2)))
# } else {
#     register(bpstart(SnowParam(2)))
# }
# Disable parallel processing
register(SerialParam())

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
  warning("To align the samples all data collections must have at least one not blank sample. ",
          "Data collection batch ", paste(unique(metadata$DATA_COLLECTION_BATCH)[
            which(!(unique(metadata$DATA_COLLECTION_BATCH) %in% 
                      unique(metadata[metadata$SAMPLE_TYPE != "blank", 
                                            "DATA_COLLECTION_BATCH"])))], 
            collapse = ", "), " do not have a not blank sample. Alignment disabled.", 
          immediate. = TRUE, call. = TRUE)
  n_samples_batch_align <- 0
}

# reorder the metada to put the blank samples first
metadata <- metadata[c(which(metadata$SAMPLE_TYPE == "blank"), 
                       which(metadata$SAMPLE_TYPE != "blank")),]


pd <- data.frame(sample_name = metadata$SAMPLE_CODE,
                 sample_group = metadata$DATA_COLLECTION_BATCH,
                 class = ifelse(metadata$SAMPLE_TYPE == "blank", "blank", "sample"),
                 stringsAsFactors = FALSE) 

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

# write.table("START", "memory_use")
# write.table(gc(), "memory_use", append = T)
# write the alignment file header, this file will store the max diff in the batches alignment
if (n_samples_batch_align > 0) {
  write.table(data.frame("batch", "samples_aligned", "rtTolMin", "rtTol1stQu","rtTolMedian", 
                         "rtTolMean", "rtTol3rdQu","rtTolMax"), sep = ",",
              file =  file.path(path_raw_data, processed_data_dir, "samples_alignment.csv"),
              row.names = FALSE, col.names = FALSE)
}
# write the header of the table that will store the MS2 mzs without a MS1 peak correspondence
write.table(data.frame("sample_code", "no_match", "mz", "rt", "int", "blank"), sep = ",",
            file =  file.path(path_raw_data, processed_data_dir,
                              "log_MS2_no_MS1peak_match.csv"),
            row.names = FALSE, col.names = FALSE)
total_spectra <- list()
blank_mzs <- c() # store the MS2 mzs that appear in blank samples
sample_count <- 0
total_ions_processed <- 0
ti <- Sys.time()
error_reading <- c()
# store the number of mzs that had 0 to 300 MS1 peaks assigned to them, or more (aggregated in the 301 position)
mzs_number_peaks <- numeric(length = 301) 
##### Match MS1 peak info with tandem MS2 data and rewrite raw data to mgf with peak annotation
# batches_rttol <- c() # max diff in the batches alignment
batches_data <- lapply(unique(metadata$DATA_COLLECTION_BATCH), function(batch)
{
  # write.table(paste("start Batch ", batch), "memory_use", append = T)
  # write.table(gc(), "memory_use", append = T)
  t1 <- Sys.time()
  gc()
  cat("\n*** Starting processing samples of data collection batch", batch, "***\n\n")
  
  # get the not blank samples indexes to keep for alignment, max to n_samples_batch_align by batch
  if (n_samples_batch_align == 0) # do not align
  {
    sample_idxs_keep <- c()
  } else {
    sample_idxs_keep <- which(metadata$DATA_COLLECTION_BATCH == batch & 
                                metadata$SAMPLE_TYPE != "blank")
  }
  # if more not blank samples in the batch then the max to keep, 
  # takes a sample from it with the max allowed size containing the first not blank sample
  if (length(sample_idxs_keep) > n_samples_batch_align)
  {
    # select the samples in sequence
    sample_idxs_keep <- sample_idxs_keep[1:n_samples_batch_align] 
                          #sample(sample_idxs_keep[-1], n_samples_batch_align-1)) # uniformely select the samples
    # sample_idxs_keep <- sample_idxs_keep[1:n_samples_batch_align]
  }
  
  # extract peak info of all samples
  raw_datas <- lapply(which(metadata$DATA_COLLECTION_BATCH == batch), function(sample_idx)
  {
    # write.table(paste("start sample idx", metadata$SAMPLE_CODE[[sample_idx]]), "memory_use", append = T)
    # write.table(gc(), "memory_use", append = T)
    gc()
    cat("\n* Extracting peak info of sample", 
            metadata$SAMPLE_CODE[[sample_idx]], "*\n")
    sample_count <<- sample_count + 1
    ti <- Sys.time()
    
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
    tryCatch(writeMgfDataFile_NP3(raw_data, 
                                  file_MGF = file.path(path_raw_data, processed_data_dir, 
                                                  paste0(metadata$SAMPLE_CODE[[sample_idx]], 
                                                         "_peak_info.mgf")), 
                                  RTMIN = tandem_peak_info$rtmin,
                                  RTMAX = tandem_peak_info$rtmax,
                                  INTO = raw_int,
                                  TITLE = metadata$SAMPLE_CODE[[sample_idx]],
                                  PEAK_ID = tandem_peak_info$id,
                                  PEAK_AREA = trunc(tandem_peak_info$intb), # TODO compare intb and $into
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
  cat("\n** Done extracting peak info of data collection batch", batch, "in", 
          round(t2-t1, 2), units(t2-t1), "**\n")
  
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
})
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
