Rcpp::sourceCpp(file.path(script_path(),
                          'read_mgf_peak_list_R.cpp'))

# MSnBase
# Based on the code contributed by Guangchuang Yu <guangchuangyu@gmail.com>
# Modified by Sebastian Gibb <mail@sebastiangibb.de>
# Modified by Cristina F. Bazzano <crisfbazz@gmail.com>
# tags: check if the provided list of tags were correctly written in all spectrum header of the created MGF
writeMgfDataFile_NP3 <- function(splist, file_MGF, COM = NULL, TITLE = NULL, RT = NULL,
                                 RTMIN = NULL, RTMAX = NULL, PEAK_ID = NULL, PEAK_AREA = NULL, SIZE = NULL,
                                 INTO = NULL, PEPMASS = NULL, CHARGE = NULL, 
                                 tags = c("BEGIN IONS", "SCANS", "RTINSECONDS", 
                                          "CHARGE", "PRECURSOR_INTENSITY", "PEPMASS",
                                          "END IONS"),
                                 verbose = isMSnbaseVerbose()) {
  if (class(file_MGF) == "character" && file.exists(file_MGF)) 
  {
    cat("Overwriting", file_MGF, "!\n")
    unlink(file_MGF, force = TRUE)
  } else {
    cat("Creating MGF", file_MGF, "!\n")
  }
  
  con <- file(description = file_MGF, open = "wb", encoding = "UTF-8")
  on.exit(close(con))
  
  if (is.null(COM)) {
    COM <- paste0("COM=", ifelse(length(splist) <= 1, "Spectrum", "Experiment"),
                  " exported by MSnbase with NP3 MS modifications on ", date())
  }
  cat(COM, file = con, sep = "")
  
  verbose <- verbose & length(splist) > 1
  
  if (verbose)
    pb <- txtProgressBar(min = 0, max = length(splist), style = 3)
  
  # writte all spectra
  sp_w <- sapply(seq_along(splist), function(i)
  {
    if (verbose)
      setTxtProgressBar(pb, i)
    
    writeMgfContent_NP3(splist[[i]], TITLE = TITLE, con = con, RT = RT[[i]],
                        RTMIN = RTMIN[[i]], RTMAX = RTMAX[[i]], PEAK_ID = PEAK_ID[[i]],
                        PEAK_AREA = PEAK_AREA[[i]], PEPMASS = PEPMASS[[i]],
                        INTO = INTO[[i]], SIZE = SIZE[[i]], CHARGE = CHARGE)
  })
  n_sp_w <- sum(sp_w) # get the number of written spectrum
  
  if (verbose)
    close(pb)
  
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
    (sum(startsWith(wmgf, tag)) == n_sp_w)
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
  scans_list <- unlist(sapply(which(sp_w == 1), function(i) {scanIndex(splist[[i]])}))
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
}

writeMgfContent_NP3 <- function(sp, TITLE = NULL, con, RT = NULL, RTMIN = NULL, 
                                RTMAX = NULL, PEAK_ID = NULL, PEAK_AREA = NULL, INTO = NULL,
                                PEPMASS = NULL, SIZE = NULL, CHARGE = NULL) {
  # do not write zero peak spectrum, or spectrum with all peaks intensity equals zero
  sp <- clean(sp, all=TRUE) # remove zero intensity peaks
  if (peaksCount(sp) == 0)
    return(0)
  
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  .cat("\nBEGIN IONS\n",
       "SCANS=", scanIndex(sp))
  
  if (is.null(TITLE)) {
    .cat("\nTITLE=", scanIndex(sp))
  } else {
    if (nchar(TITLE) > 200) 
      TITLE <- substr(TITLE, 1, 200)
    .cat("\nTITLE=", TITLE,".",scanIndex(sp))
  }
  
  if (!is.null(RTMIN)) {
    # if outside the peak profile, enlarge peak
    if (rtime(sp) < RTMIN)
      RTMIN <- rtime(sp) - 0.5
    
    .cat("\nRTMIN=", RTMIN)
  }
  
  if (!is.null(RTMAX)) {
    # if outside the peak profile, enlarge peak
    if (rtime(sp) > RTMAX)
      RTMAX <- rtime(sp) + 0.5
    
    .cat("\nRTMAX=", RTMAX)
  }
  
  if (!is.null(PEAK_ID)) {
    .cat("\nPEAK_ID=", PEAK_ID)
  }
  
  if (!is.null(PEAK_AREA)) {
    .cat("\nPEAK_AREA=", PEAK_AREA)
  }
  
  if (!is.null(SIZE)) {
    .cat("\nCLUSTER_SIZE=", SIZE)
  }
  
  #.cat("\nGROUPID=", acquisitionNum(sp))
  
  if (msLevel(sp) > 1) {
    if (length(precursorCharge(sp)) && !is.na(precursorCharge(sp))) {
      .cat("\nCHARGE=", ifelse(is.null(CHARGE), precursorCharge(sp), CHARGE))
    }
    .cat("\nNUM_PEAKS=", peaksCount(sp), 
         "\nRTINSECONDS=", ifelse(is.null(RT), rtime(sp), RT), 
         "\nPRECURSOR_INTENSITY=", ifelse(is.null(INTO), 
                                          ifelse(precursorIntensity(sp)>0, 
                                                 precursorIntensity(sp), 
                                                 ionCount(sp)), 
                                          INTO),
         "\nPEPMASS=", ifelse(is.null(PEPMASS), precursorMz(sp), PEPMASS))
  } else .cat("\nRTINSECONDS=", rtime(sp))
  
  .cat("\n", paste(mz(sp), intensity(sp), collapse = "\n"))
  .cat("\nEND IONS\n")
  
  return(1)
}

##### CONCATENATION FUNS

concatenate_OnDiskMSnExp <- function(...) {
  x <- list(...)
  if (length(x) == 0)
    return(NULL)
  if (length(x) == 1)
    return(x[[1]])
  ## Check that all are XCMSnExp objects.
  if (!all(unlist(lapply(x, function(z) is(z, "OnDiskMSnExp")))))
    stop("All passed objects should be 'OnDiskMSnExp' objects")
  ## Check processingQueue
  procQ <- lapply(x, function(z) z@spectraProcessingQueue)
  new_procQ <- procQ[[1]]
  is_ok <- unlist(lapply(procQ, function(z)
    !is.character(all.equal(new_procQ, z))
  ))
  if (any(!is_ok)) {
    warning("Processing queues from the submitted objects differ! ",
            "Dropping the processing queue.")
    new_procQ <- list()
  }
  ## processingData
  fls <- lapply(x, function(z) z@processingData@files)
  startidx <- cumsum(lengths(fls))
  ## featureData
  featd <- lapply(x, fData)
  ## Have to update the file index and the spectrum names.
  for (i in 2:length(featd)) {
    featd[[i]]$fileIdx <- featd[[i]]$fileIdx + startidx[i - 1]
    rownames(featd[[i]]) <- MSnbase:::formatFileSpectrumNames(
      fileIds = featd[[i]]$fileIdx,
      spectrumIds = featd[[i]]$spIdx,
      nSpectra = nrow(featd[[i]]),
      nFiles = length(unlist(fls))
    )
  }
  featd <- do.call(rbind, featd)
  featd$spectrum <- 1:nrow(featd)
  ## experimentData
  expdata <- lapply(x, function(z) {
    ed <- z@experimentData
    data.frame(instrumentManufacturer = ed@instrumentManufacturer,
               instrumentModel = ed@instrumentModel,
               ionSource = ed@ionSource,
               analyser = ed@analyser,
               detectorType = ed@detectorType,
               stringsAsFactors = FALSE)
  })
  expdata <- do.call(rbind, expdata)
  expdata <- new("MIAPE",
                 instrumentManufacturer = expdata$instrumentManufacturer,
                 instrumentModel = expdata$instrumentModel,
                 ionSource = expdata$ionSource,
                 analyser = expdata$analyser,
                 detectorType = expdata$detectorType)
  
  ## protocolData
  protodata <- lapply(x, function(z) z@protocolData)
  if (any(unlist(lapply(protodata, nrow)) > 0))
    warning("Found non-empty protocol data, but merging protocol data is",
            " currently not supported. Skipped.")
  ## phenoData
  pdata <- do.call(rbind, lapply(x, pData))
  res <- new(
    "OnDiskMSnExp",
    phenoData = new("NAnnotatedDataFrame", data = pdata),
    featureData = new("AnnotatedDataFrame", featd),
    processingData = new("MSnProcess",
                         processing = paste0("Concatenated [", date(), "]"),
                         files = unlist(fls), smoothed = NA),
    experimentData = expdata,
    spectraProcessingQueue = new_procQ)
  if (is(res, "OnDiskMSnExp"))
    res
}

concatenate_XCMSnExp <- function(...) {
  x <- list(...)
  if (length(x) == 0)
    return(NULL)
  if (length(x) == 1)
    return(x[[1]])
  ## Check that all are XCMSnExp objects.
  if (!all(unlist(lapply(x, function(z) is(z, "XCMSnExp")))))
    stop("All passed objects should be 'XCMSnExp' objects")
  new_x <- as(concatenate_OnDiskMSnExp(...), "XCMSnExp")
  ## If any of the XCMSnExp has alignment results or detected features drop
  ## them!
  x <- lapply(x, function(z) {
    if (hasAdjustedRtime(z)) {
      z <- dropAdjustedRtime(z)
      warning("Adjusted retention times found, had to drop them.")
    }
    if (hasFeatures(z)) {
      z <- dropFeatureDefinitions(z)
      warning("Feature definitions found, had to drop them.")
    }
    z
  })
  ## Combine peaks
  fls <- lapply(x, fileNames)
  startidx <- cumsum(lengths(fls))
  pks <- lapply(x, chromPeaks)
  procH <- lapply(x, processHistory)
  for (i in 2:length(fls)) {
    pks[[i]][, "sample"] <- pks[[i]][, "sample"] + startidx[i - 1]
    procH[[i]] <- lapply(procH[[i]], function(z) {
      z@fileIndex <- as.integer(z@fileIndex + startidx[i - 1])
      z
    })
  }
  pks <- do.call(rbind, pks)
  new_x@.processHistory <- unlist(procH)
  chromPeaks(new_x) <- pks
  
  # message("DONE CONCAT XCMSExp")
  if (is(new_x, "XCMSnExp"))
    new_x
  else 
    NULL
}

str_list_to_numeric <- function(str_list, sep = ';') {
  as.numeric(strsplit(str_list, sep)[[1]])
}

inverse_scale_numeric <- function(ints, scale_factor)
{
  ints <- str_list_to_numeric(ints)
  if (scale_factor == 0) # log scale was applied
  {
    # inverse of 1.0 + log(1.0 + multVal * ints[i]))
    ints <- expm1(ints-1)
  } else {
    ints <- ints^(1/scale_factor)
  }
  return(round(ints,5))
}

writeMgfDataFile_NP3_table <- function(ms_count_table, file_MGF, output_name, 
                                       scale_factor=0.5, charge='+') {
  tags <- c("BEGIN IONS", "SCANS", "RTINSECONDS",
           "RTMIN", "RTMAX", "CLUSTER_SIZE",
           "CHARGE", "PRECURSOR_INTENSITY", "PEPMASS",
           "END IONS")
  if (class(file_MGF) == "character" && file.exists(file_MGF)) 
  {
    cat("Overwriting", file_MGF, "!\n")
    unlink(file_MGF, force = TRUE)
  } else {
    cat("Creating MGF", file_MGF, "!\n")
  }
  
  con <- file(description = file_MGF, open = "wb", encoding = "UTF-8")
  on.exit(close(con))
  
  COM <- paste0("COM=", ifelse(nrow(ms_count_table) <= 1, "Spectrum", "Experiment"),
                  " exported by NP3 MS workflow on ", date())
  cat(COM, file = con, sep = "")
  
  verbose <- nrow(ms_count_table) > 1
  
  if (verbose)
    pb <- txtProgressBar(min = 0, max = nrow(ms_count_table), style = 3)
  
  # writte all spectra
  sp_w <- sapply(seq_len(nrow(ms_count_table)), function(i)
  {
    if (verbose)
      setTxtProgressBar(pb, i)
    
    writeMgfContent_NP3_table(peaksList = str_list_to_numeric(ms_count_table$peaksList[[i]]), 
                        peaksInt = inverse_scale_numeric(ms_count_table$peaksInt[[i]], scale_factor), 
                        SCANS = as.integer(ms_count_table$msclusterID[[i]]), 
                        TITLE = output_name,
                        con = con, 
                        RT = ms_count_table$rtMean[[i]],
                        RTMIN = ms_count_table$rtMin[[i]], 
                        RTMAX = ms_count_table$rtMax[[i]], 
                        PEPMASS = ms_count_table$mzConsensus[[i]],
                        INTO = ms_count_table$sumInts[[i]], 
                        SIZE = ms_count_table$numSpectra[[i]], 
                        CHARGE = charge)
  })
  n_sp_w <- sum(sp_w) # get the number of written spectrum
  
  if (verbose)
    close(pb)
  
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
    (sum(startsWith(wmgf, tag)) == n_sp_w)
  })
  if (!all(check_tags))
  {
    stop("ERROR. Problem with the written MGF, the following tags do not match ",
         "with the count of written spectrum: ", 
         paste(names(check_tags)[!check_tags], collapse = ", "), 
         ". Error probably due to a connection problem, process aborted. ", 
         "Rerun the process to overwritte the corrupted file.")
  }
  # check if the written scans number match the msclusterIDs list,
  # all scans must have been written here
  mgf_header <- readMgfHeader(file_MGF)
  if (anyNA(match(ms_count_table$msclusterID, mgf_header$scans))) 
  {
    stop("ERROR. Problem with the written MGF, the following msclusterIDs are ",
         "missing from the list of written scans: ",
         paste(ms_count_table$msclusterID[which(
           is.na(match(ms_count_table$msclusterID, 
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
}

writeMgfContent_NP3_table <- function(peaksList, peaksInt, SCANS = NULL, 
                                      TITLE = NULL, con, 
                                      RT = NULL, RTMIN = NULL, RTMAX = NULL, 
                                      INTO = NULL, PEPMASS = NULL, SIZE = NULL, 
                                      CHARGE = NULL) {
  if (length(peaksList) == 0) # do not write zero peak spectrum
    return(0)
  
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  .cat("\nBEGIN IONS\n",
       "SCANS=", SCANS)
  
  if (!is.null(TITLE)) {
    # max title size is 200
    if (nchar(TITLE) > 200) 
      TITLE <- substr(TITLE, 1, 200)
    .cat("\nTITLE=", TITLE,".",SCANS)
  }
  
  if (!is.null(RTMIN)) {
    .cat("\nRTMIN=", RTMIN)
  }
  
  if (!is.null(RTMAX)) {
    .cat("\nRTMAX=", RTMAX)
  }
  
  if (!is.null(SIZE)) {
    .cat("\nCLUSTER_SIZE=", SIZE)
  }

  if (!is.null(CHARGE)) {
    .cat("\nCHARGE=", ifelse(CHARGE=="+","1+", "1-"))
  }
  .cat("\nNUM_PEAKS=", length(peaksList), 
       "\nRTINSECONDS=", RT, 
       "\nPRECURSOR_INTENSITY=", INTO,
       "\nPEPMASS=", PEPMASS)
  
  
  .cat("\n", paste(peaksList, peaksInt, collapse = "\n"))
  .cat("\nEND IONS\n")
  
  return(1)
}
