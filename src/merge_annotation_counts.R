## ----load-libs, message = FALSE--------------------------------------------
cat("Loading packages readr and auxiliary functions...\n")
library(readr)
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
source(file.path(script_path(), "count_peak_area.R"))

# function to merge the count table by column for all annotations 
merge_counts <- function(col_name, x)
{
  # print(col_name)
  # x <- unlist(x,recursive = F)
  switch(col_name,
         msclusterID=,mzConsensus=,rtMean=,rtMin=,rtMax=,
         peaksList=,peaksInt=,adducts=,isotopes=,
         dimers=,multiCharges=,fragments=,analogs=,
         multicharge_ion=,isotope_ion=
           x[[col_name]][[1]],
         BLANK_DIST=ifelse(any(!is.na(x[[col_name]])), 
                           min(x[[col_name]], na.rm = TRUE), 
                           NA),
         BEDFLAG=,BFLAG=,CFLAG=any(as.logical(x[[col_name]])),
         HFLAG=,DESREPLICATION=,peakIds=
           ifelse(any(!is.na(x[[col_name]])), 
                  paste(unique(unlist(strsplit(x[[col_name]][
                    !is.na(x[[col_name]])], ";"))), collapse = ";"), 
                  NA),
         precursorMz=,scans=,joinedJobsIDs= paste(x[[col_name]], collapse = ";"),
         sumInts =,BLANKS_TOTAL =,BEDS_TOTAL=,CONTROLS_TOTAL=,numSpectra=,numJoins=
           sum(x[[col_name]]),
         basePeakInt=max(as.numeric(x[[col_name]])),
         joinedIDs=
           ifelse(any(!is.na(x[[col_name]])), # if there is a not NA value paste it, delim = ;
                  paste(x[[col_name]][!is.na(x[[col_name]])], collapse = ";"), 
                  NA),
         gnps_Smiles=,tremolo_SMILES=ifelse(any(!is.na(x[[col_name]])), # if there is a not NA value paste it, delim = ,
                               paste(x[[col_name]][!is.na(x[[col_name]])], collapse = ","), 
                               NA),
         mergedIDs_all= unlist(x[[col_name]]),
         mergedIDs_adducts=,precursorMz_adducts=,mergedIDs_isotopes=,precursorMz_isotopes=,
         mergedIDs_dimers=,precursorMz_dimers=,mergedIDs_multiCharges=,
         precursorMz_multiCharges=,mergedIDs_fragments=,precursorMz_fragments=
           ifelse(any(!is.na(x[[col_name]])),
                  paste(x[[col_name]][!is.na(x[[col_name]])], collapse = ";"),
                  NA),
         mergedIDs = paste(c(x$msclusterID[is.na(x[[col_name]])],
                             x[[col_name]][!is.na(x[[col_name]])]), 
                           collapse = ";"),
         ifelse(is.numeric(x[[col_name]]) && !startsWith(col_name, 'gnps_'), # do not sum gnps results, concatenate them
                sum(x[[col_name]], na.rm = TRUE), # sum all the counts cols (_area and _spectra)
                ifelse(any(!is.na(x[[col_name]])), # if there is a not NA value paste it, delim = ;
                       paste(x[[col_name]][!is.na(x[[col_name]])], collapse = ";"), 
                       NA))) # cocatenate string fields (e.g. from tremolo or gnps)
}

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Two arguments must be supplied to merge the clean count tables based on an annotation\n", 
       " 1 - Path to the output data folder, inside the outs directory of the clustering result folder. ", 
       "It should contain the count_tables folder with the clean and annotated counts. The data name will be extracted from here;\n", 
       " 2 - The annotations types to be used to merge the count tables separated by a comma: adducts,isotopes,dimers,multiCharges,fragments;\n",
       " 3 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks;\n",
       " 4 - The pre processed data path were the MGFs were created. Used to compute the peak areas;\n",
       " 5 - A boolean TRUE or FALSE indicating if only the protonated representative spectra should be merged. If FALSE merge all msclusterIDs.\n",
       call.=FALSE)
} else {
  output_path <- file.path(args[[1]])
  if (!dir.exists(output_path))
  {
    stop("The job output folder '", output_path, 
         "' do not exists. Provide a valid path to where the the job final result is located.")
  }
  output_name <- basename(output_path)
  
  path_spectra_count <- file.path(output_path, "count_tables", "clean", 
                                  paste0(output_name,"_spectra_clean_ann.csv"))
  if (!file.exists(path_spectra_count))
  {
    stop("The spectra count file '", path_spectra_count,
         "' do not exists. Provide a valid output path to where csv file with the count ",
         "of spectra is located.")
  }
  
  annotation_merge <- tolower(strsplit(args[[2]], ",")[[1]])
  if (!any(annotation_merge %in% c("adducts", "isotopes", "dimers",	
                                   "multicharges", "fragments")))
  {
    stop("The annotation types to be used to merge the count tables are invalid: '", 
         paste(annotation_merge, collapse=", "),
         "'. It should be one of: adducts, isotopes, dimers,	multicharges or fragments.")
  } else if (!all(annotation_merge %in% c("isotopes", "adducts", "dimers", 
                                          "multicharges", "fragments"))) {
    warning("Some of the provided annotation types are not valid. ", 
            "Only the following will be merged: ", paste(annotation_merge, collapse=", "),
            call. = FALSE)
  } 
  # order anns types
  annotation_merge <- c("isotopes", "adducts", "dimers", 
                        "multicharges", "fragments")[
                          c("isotopes", "adducts", "dimers", "multicharges", "fragments") %in% 
                            annotation_merge]
  
  # convert multicharges names to the expected format: multiCharges
  annotation_merge[annotation_merge=="multicharges"] <- "multiCharges"
  
  path_batch_metadata <- file.path(args[[3]])
  if (!file.exists(path_batch_metadata))
  {
    stop("The CSV batch metadata file '", path_batch_metadata, 
         "' do not exists. Provide a valid path to where the metadata is located.")
  }
  
  processed_data_path <- file.path(args[[4]])
  if (!dir.exists(processed_data_path))
  {
    stop("The processed data folder '", processed_data_path, 
         "' do not exists. Provide a valid path to where the pre processed MGFs are located.")
  }
  
  merge_protonated_representative <- as.logical(toupper(args[[5]]))
  if (is.na(merge_protonated_representative))
  {
    warning("The argument indicating if only the protonated representative ",
            "spectra should be merged had a wrong values and was converted to ",
            "NA. Setting it to TRUE by default.")
    merge_protonated_representative <- TRUE
  }
}

# read the spectra count table
ms_spectra_count <- suppressMessages(read_csv(path_spectra_count, guess_max = 5000))
if (!all(annotation_merge %in% names(ms_spectra_count))) {
  stop("The spectra count file '", path_spectra_count,
       "' do not have the annotation column '", annotation_merge, 
       "' to be used in the merge. Provide a valid path to where the csv ",
       "file with the annotated count of spectra is located.")
}

#  round up to the nearest power of 10
anns_start_idx_i <- anns_start_idx <- 10^ceiling(log10(max(ms_spectra_count$msclusterID)))

# # filter only the protonated reresentative msclusterID, if merge_protonated_representative is TRUE
if (merge_protonated_representative) {
  merge_protonated_representative <- 1
} else {
  merge_protonated_representative <- 0
}

# create all merged IDs column as a list
ms_spectra_count$mergedIDs_all <- lapply(ms_spectra_count$msclusterID, function(x) x)

# get the blank IDs
if ("BLANKS_TOTAL" %in% names(ms_spectra_count)) {
  blank_scans <- ms_spectra_count$msclusterID[ms_spectra_count$BLANKS_TOTAL > 0]
} else {
  blank_scans <- NULL
}


cat("\n** Merging counts based on the",  paste(annotation_merge, collapse=", "), "annotations **\n\n")
t0 <- Sys.time()
cat(paste0("    |", paste0(rep(" ", length(annotation_merge)), collapse = ""), "|\n    |", collapse = ""))
n_symbolic_clusters <- 0
for (ann_type in annotation_merge) {
  cat("=")
  # cat(ann_type)
  # create idx column 
  ms_spectra_count$idx <- 1:nrow(ms_spectra_count)
  ms_spectra_count$mergedIDs <- NA
  ms_spectra_count$precursorMz <- ms_spectra_count$mzConsensus
  
  # only merge annotations of not blank mzs, when there is a blank sample
  if (!is.null(blank_scans)) {
    ann_idxs <- ms_spectra_count$idx[!is.na(ms_spectra_count[[ann_type]]) & 
                                       ms_spectra_count$BLANKS_TOTAL == 0 &
                                       (ms_spectra_count$protonated_representative >= merge_protonated_representative)]
  } else {
    ann_idxs <- ms_spectra_count$idx[!is.na(ms_spectra_count[[ann_type]]) &
                                       (ms_spectra_count$protonated_representative >= merge_protonated_representative)]
  }
  
  if (length(ann_idxs) == 0) # no annotations of this type
    next()
    
  # do not add ann target if it is already present in the list of merged ids of the source
  # do not add blanks ann targets
  annotation_info <- gtools::smartbind(list=lapply(ann_idxs, function(x)
  {
    # cat(x)
    scan_source <- ms_spectra_count$msclusterID[[x]]
    anns <- strsplit(ms_spectra_count[[ann_type]][[x]], ";")[[1]]
    scan_targets <- unique(as.numeric(regmatches(anns,
                                     regexpr("(?<=\\[)[0-9]+(?=\\])",
                                             anns, perl = T))))
    # get all merged ids without symbolic
    mergedIds_real <- ms_spectra_count$mergedIDs_all[[x]][ms_spectra_count$mergedIDs_all[[x]] < anns_start_idx_i]
    
    # rm repeated targets
    # valid_targets <- which(!(scan_targets %in% c(blank_scans, mergedIds_real)))
    valid_targets <- which(!(scan_targets %in% c(mergedIds_real)))
    
    if (length(valid_targets) == 0) {
      ann_info <- data.frame(idx_source = x,
                             idx_target = paste0(match(scan_targets, ms_spectra_count$msclusterID), collapse=";"),
                             stringsAsFactors = FALSE)
    } else {
      # filter only the valid targets
      scan_targets <- scan_targets[valid_targets]
      anns <- anns[valid_targets]
      
      ann_info <- data.frame(
           idx_source = x,
           idx_target = paste0(match(scan_targets, ms_spectra_count$msclusterID), collapse=";"),
           degree = length(anns),
           stringsAsFactors = FALSE)
    }
    ann_info
  }))
  # remove not valid NAs
  annotation_info <- na.omit(annotation_info)
  annotation_info$idx_target <- lapply(annotation_info$idx_target, 
                                       function(x) as.numeric(strsplit(x, split=";")[[1]]))
  # no annotation of this type, go to next
  if (nrow(annotation_info) == 0)
    next()
  
  # suppresswarning of implicit type convertion to chars 
  # merge all annotations present in the info table
  merged_counts <- do.call(rbind,
    lapply(seq_len(nrow(annotation_info)), function(i)
    {
      merged_targets <- list(lapply(names(ms_spectra_count), merge_counts,
                                    ms_spectra_count[c(annotation_info$idx_source[[i]],
                                                    annotation_info$idx_target[[i]]),]))
    
      merged_targets <- do.call(rbind, merged_targets)
      colnames(merged_targets) <- names(ms_spectra_count)
      
      merged_targets
    }))
  
  merged_counts[,1] <- anns_start_idx:(anns_start_idx+nrow(merged_counts)-1)
  anns_start_idx <- anns_start_idx + anns_start_idx + nrow(merged_counts) - anns_start_idx
  # add symbolic cluster id to the mergedIDs
  merged_counts[,"mergedIDs_all"] <- lapply(seq_len(nrow(merged_counts)), 
                                            function(i) 
                                              unique(c(merged_counts[[i,"mergedIDs_all"]], merged_counts[[i,1]])))
  
  # at least one merge by row of the annotaiton info table
  if (nrow(merged_counts) < nrow(annotation_info))
    stop("Something went wrong in the merging step, missing new symbolic candidates.")
  
  n_symbolic_clusters <- n_symbolic_clusters + nrow(merged_counts)
  ms_spectra_count <- rbind(ms_spectra_count, merged_counts, deparse.level = 0)
  ms_spectra_count[,1:ncol(ms_spectra_count)] <- lapply(1:ncol(ms_spectra_count), 
                                                  function(j, mergedIDs_all_idx) 
                                                    {
                                                    if (j != mergedIDs_all_idx)
                                                      unlist(ms_spectra_count[,j], use.names = FALSE)
                                                    else
                                                      ms_spectra_count$mergedIDs_all
                                                    },
                                                  which(names(ms_spectra_count) == "mergedIDs_all")) 
  ms_spectra_count$precursorMz[is.na(ms_spectra_count$mergedIDs)] <- NA
  
  names(ms_spectra_count)[match(c("mergedIDs","precursorMz"),names(ms_spectra_count))] <- 
    paste0(c("mergedIDs_","precursorMz_"), ann_type)
}
cat("|\n\n")
tf <- Sys.time()
cat("  * Done creating", n_symbolic_clusters, "symbolic clusters in", 
        round(tf-t0, 2), units(tf-t0), "*\n")

if (n_symbolic_clusters == 0) {
  q('no')
}

#check merge
check_merged_IDs <- sapply(seq_len(nrow(ms_spectra_count)), function(i, mergedIDs_idx) {
  merged_IDs_i <- unlist(ms_spectra_count[i, mergedIDs_idx])
  merged_IDs_i <- unique(c(as.numeric(unlist(strsplit(merged_IDs_i[!is.na(merged_IDs_i)], ";"))),
                           ms_spectra_count[[i, 1]]))
  if (length(merged_IDs_i) == 0)
    return(TRUE)
  all(merged_IDs_i == ms_spectra_count$mergedIDs_all[[i]])
}, which(startsWith(names(ms_spectra_count), "mergedIDs") & 
           names(ms_spectra_count) != "mergedIDs_all"))
if (!all(check_merged_IDs))
  stop("Missing mergedIDs for some clusters. Something went wrong in the merge step.")

# ms_area_count$mergedIDs_all <- NULL
ms_spectra_count$mergedIDs_all <- sapply(ms_spectra_count$mergedIDs_all, paste, collapse = ";")
ms_spectra_count$idx <- NULL

# create dir and save spectra count merged
if (!dir.exists(file.path(output_path, "count_tables", "merge")))
  dir.create(file.path(output_path, "count_tables", "merge"), showWarnings = FALSE)

# compute peak area
ms_area_count <- ms_spectra_count
count_columns <- which(endsWith(names(ms_area_count), "_spectra"))
names(ms_area_count)[count_columns] <- sub(pattern = "_spectra",
                                           replacement = "_area", fixed = TRUE,
                                           x = names(ms_area_count)[count_columns])
# check if this is a joined jobs, if yes call compute peak area from joined jobs
if (!("joinedJobsIDs" %in% names(ms_spectra_count))) {
  # this is the default flow of np3, compute the peak areas using the provided
  # pre processed data path
  batch_metadata <- readMetadataTable(path_batch_metadata)
  cat("\n  ** Computing Peak Area by Sample ** \n\n")
  peak_areas_base_peak_int <- compute_peak_area(processed_data_path,
                                                ms_area_count$msclusterID,
                                                lapply(ms_area_count$scans, function(x) strsplit(x, ";")[[1]]),
                                                lapply(ms_area_count$peakIds, function(x) strsplit(x, ";")[[1]]),
                                                batch_metadata)
} else {
  # this is the join_jobs flow, compute the peak areas using all the original samples
  # pre processed data
  cat("\n  ** Computing Peak Area by Original Sample of the Joined jobs ** \n\n")
  peak_areas_base_peak_int <- compute_peak_areas_joined_jobs(output_path, 
                                                             ms_area_count$msclusterID,
                                                             ms_area_count$scans,
                                                             ms_area_count$peakIds, 
                                                             ms_area_count$joinedJobsIDs)
}
# order the peak area columns from the metadata with the order present in the count table
ms_area_count[,count_columns] <- peak_areas_base_peak_int[,
                                                          match(names(ms_area_count)[count_columns],
                                                                names(peak_areas_base_peak_int)[-ncol(peak_areas_base_peak_int)])]
ms_spectra_count$basePeakInt <- ms_area_count$basePeakInt <- peak_areas_base_peak_int$basePeakInt

write_csv(ms_spectra_count, path = file.path(output_path, "count_tables", "merge", 
                                             paste0(output_name, "_spectra_merged_ann.csv")))
write_csv(ms_area_count, path = file.path(output_path, "count_tables", "merge",
                                          paste0(output_name, "_peak_area_merged_ann.csv")))

t0 <- Sys.time()
cat("  * Done in", 
        round(t0-tf, 2), units(t0-tf), "*\n")
