decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

# check if a MGF file have all scans in the provided list with all provided tags in the header
check_mgf <- function(file_MGF, tags, scans_list) {
  cat("Checking created MGF ", file_MGF, " .")
  # check if the written mgf was created correctly
  wmgf <- tryCatch(readLines(file_MGF),
                   warning = function(w) {
                     message("Warning opening file ", file_MGF, ". \n", w)
                   })
  # check if the count of tags is correct with he count of written spectrum
  check_tags <- sapply(tags, function(tag) 
  {
    cat(".")
    (sum(startsWith(wmgf, tag)) == n)
  })
  if (!all(check_tags))
  {
    stop("ERROR. Problem with the written MGF, the following tags do not match ",
         "with the count of written spectrum: ", 
         paste(names(check_tags)[!check_tags], collapse = ", "), 
         ". Error probably due to a connection problem, process aborted. ", 
         "Rerun the process to overwritte the corrupted file.")
  }
  rm(wmgf)
  # check if the written scans from the mgf match the msclusterIDs list from the scans list
  # remove the scans that were not written because of few peaks (sp_w == 0)
  mgf_header <- readMgfHeader(file_MGF)
  if (anyNA(match(scans_list, mgf_header$scans))) 
  {
    stop("ERROR. Problem with the written MGF, the following msclusterIDs are ",
         "missing from the list of written scans: ",
         paste(scans_list[which(is.na(match(scans_list,
                                            mgf_header$scans)))], collapse=","), 
         ". Please report this inconsistency to the development team, ",
         "a bug probable ocurried when writting the mgf file.")
  }
  # check if there is any duplicated scans, these must be unique values
  if (any(duplicated(mgf_header$scans))) {
    stop("ERROR. Problem with the written MGF, duplicated scans were written.",
         "The following scans have duplicated values: ", 
         paste(mgf_header$scans[duplicated(mgf_header$scans)],collapse=","), 
         ". Please report this inconsistency to the development team, ",
         "a bug probable ocurried when writting the mgf file.")
  }
  cat(" OK!\n")
  
  tcheck <- Sys.time()
  cat("    * Done in", round(tcheck-tend, 2), units(tcheck-tend), "*\n")
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
  # filter only the samples present in the current job
  # use only the samples preprocessed in this execution, which are present in the total_spectra list
  mzs_nomatch <- mzs_nomatch[mzs_nomatch$sample_code %in% names(total_spectra),]
  
  cat("\n\n* Pre-processing - Statistics of the percentage of MS2 spectra without a MS1 peak correspondence *\n\n")
  if (nrow(mzs_nomatch) == 0) {
    cat("Great, all MS2 spectra had a MS1 peak correspondence!!")
    return()
  }
  
  # add blank information to the list of ms1 peaks with a ms2 correspondence
  # and to the log of ms2 without ms1
  ms1_w_ms2_list <- suppressMessages(readr::read_csv(file.path(path_raw_data, processed_data_dir,
                                                               "MS1_list_with_MS2.csv"), guess_max = 5000))
  ms1_w_ms2_list$BFLAG <- FALSE
  ms1_w_ms2_list[ms1_w_ms2_list$rtMin == 0 & ms1_w_ms2_list$rtMax == 1000000, 
                 "BFLAG"] <- TRUE
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
  # if not include the missing ones without no correspondence (good pre-processed samples)
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
