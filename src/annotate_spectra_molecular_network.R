##
# Step 7 - performs the pairwise annotation of concurrent consensus spectra and 
# creates the IVAMN
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

cat("Loading packages Rcpp, readr, dplyr...\n")
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
source(file.path(script_path(),"read_metadata_table.R"))
Rcpp::sourceCpp(file.path(script_path(),
                          'norm_dot_product.cpp'))
Rcpp::sourceCpp(file.path(script_path(),
                          'triangular_matrix_R.cpp'))


# annotations values
iso_mass <- 1.0033  # mass (13C) - mass (12C)
h2o_mass <- 18.01056
nh3_mass <- 17.03052
h_e_mass <- 1.00783 # mass of proton + mass of electron
pattern_ditrimer_adduct_isotope_multiIsotopic <- "\\[[0-9]M(\\+|-).*|\\[M(\\+|-)(?!.*\\][0-9](\\+|-)|H-H2O|H-NH3).*\\](\\+|-)|.*isotopic$"

RMSE <- function(x, y){
  x <- as.numeric(x)
  y <- as.numeric(y)
  sqrt(mean((x - y)^2))
}

annFormat <- function(ann, sim, mzError, rtError, variant_ID, num_common_samples) {
  mzError <- round(as.numeric(mzError), 5)
  variant_ID <- as.integer(trimws(variant_ID))
  paste0(ann, " (sim ", sim, " - mzE ", mzError, " - rtE ", rtError,")", 
         "[", variant_ID, "]", 
         "{", num_common_samples, "}")
}

save_annotation_net <- function(i, scans_annotations, scans_order, output_path,
                                output_name) {
  #print(i)
  scan_num <- as.integer(scans_annotations[[1]][[i]])
  annotations <- scans_annotations[[2]][[i]]
  annotations_i <- list(adducts = NA,
                        isotopes = NA,
                        dimers = NA,
                        multiCharges = NA,
                        fragments = NA,
                        analogs = NA)
  
  if (!is.na(annotations))
  {
    annotations <- strsplit(annotations, ";")[[1]]
    ann_idx <- match(as.numeric(regmatches(annotations, 
                                           regexpr("(?<=\\[)[0-9]+(?=\\])", 
                                                   annotations, perl = T))), 
                     scans_order)
    ann <- regmatches(annotations, 
                      regexpr("^.*(?= \\(.*\\).*)", 
                              annotations, perl = T))
    
    for (j in seq_along(ann_idx)) {
      # fill annotations by type of variant
      if (grepl(pattern = "\\[(2M|3M)", ann[[j]])) {
        annotations_i$dimers <- ifelse(is.na(annotations_i$dimers), 
                                       annotations[[j]], 
                                       paste(annotations_i$dimers, annotations[[j]], sep = ";"))
      } else if (grepl(pattern = "\\](2\\+|3\\+)", ann[[j]])) {
        annotations_i$multiCharges <- ifelse(is.na(annotations_i$multiCharges), 
                                             annotations[[j]], 
                                             paste(annotations_i$multiCharges, annotations[[j]], sep = ";")) 
      } else if (grepl(pattern = "fragment", ann[[j]])) {
        annotations_i$fragments <- ifelse(is.na(annotations_i$fragments), 
                                          annotations[[j]], 
                                          paste(annotations_i$fragments, annotations[[j]], sep = ";"))
      } else if (grepl(pattern = "analog", ann[[j]])) {
        annotations_i$analogs <- ifelse(is.na(annotations_i$analogs), 
                                        annotations[[j]], 
                                        paste(annotations_i$analogs, annotations[[j]], sep = ";"))
      } else if (grepl(pattern = "\\[M\\+[0-9]+\\]", ann[[j]])) {
        annotations_i$isotopes <- ifelse(is.na(annotations_i$isotopes), 
                                         annotations[[j]], 
                                         paste(annotations_i$isotopes, annotations[[j]], sep = ";"))
      # TODO add neutral loss, now it is added to the adducts
      #  } else if (grepl(pattern = "\\[M\\+H-(H2O|NH3|NH3-H2O)\\]", ann[[j]])) {
      #   annotations_i$neutral_loss <- ifelse(is.na(annotations_i$isotopes), 
      #                                    annotations[[j]], 
      #                                    paste(annotations_i$isotopes, annotations[[j]], sep = ";"))
      } else { # adduct and neutral loss
        annotations_i$adducts <- ifelse(is.na(annotations_i$adducts), 
                                        annotations[[j]], 
                                        paste(annotations_i$adducts, annotations[[j]], sep = ";")) 
      }
    }
    rm_analog <- !grepl(pattern = "analog", ann)
    if (any(rm_analog)) {
      ann <- ann[rm_analog]
      annotations <- annotations[rm_analog]
      # write annotation network rows
      ann_protonate_id_sources <- as.integer(regmatches(annotations, 
                                                        regexpr("(?<=\\[)[0-9]+(?=\\])", 
                                                                annotations, perl = T)))
      ann_sim <- as.numeric(regmatches(annotations, 
                                       regexpr("(?<=\\(sim )[0-9]+(\\.?[0-9]*|e\\-[0-9]+)(?= \\- mzE)", 
                                               annotations, perl = T)))
      ann_mzError <- as.numeric(regmatches(annotations, 
                                           regexpr("(?<= mzE )[0-9]+(\\.?[0-9]*|e\\-[0-9]+)(?= \\- rtE)", 
                                                   annotations, perl = T)))
      ann_rtError <- as.numeric(regmatches(annotations, 
                                           regexpr("(?<= rtE )[0-9]+(\\.?[0-9]*|e\\-[0-9]+)(?=\\))", 
                                                   annotations, perl = T)))
      ann_numCommonSamples <- as.numeric(regmatches(annotations, 
                                           regexpr("(?<=\\{)[0-9]+(?=\\})", 
                                                   annotations, perl = T)))
      # the current scan is the target
      write.table(data.frame(ann_protonate_id_sources, scan_num, ann_sim, ann,
                             ann_mzError, ann_rtError, ann_numCommonSamples), 
                  file.path(output_path, paste0("molecular_networking/", 
                                                output_name,
                                                "_ivamn.selfloop")),
                  row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",", 
                  fileEncoding = "UTF-8")
    }
  }
  annotations_i
}

write_annotation_net_singletons_join_duplicated <- function(output_path, 
                                                            output_name, 
                                                            ms_area_clean) {
  ann_net <- suppressMessages(read_csv(file.path(output_path, 
                                                 paste0("molecular_networking/", 
                                                        output_name,
                                                        "_ivamn.selfloop")),
                                       guess_max = 5000))
  # join duplicated annotation
  duplicated_ann <- which(duplicated(ann_net[,c(1,2)]))
  if (length(duplicated_ann) > 0) {
    for (duplicated_row in duplicated_ann) {
      keeping_row <- which(ann_net$msclusterID_source == ann_net$msclusterID_source[[duplicated_row]] & 
        ann_net$msclusterID_target == ann_net$msclusterID_target[[duplicated_row]])[1]
      ann_net[keeping_row,"annotation"] <- paste(ann_net[keeping_row,"annotation"],
                                            ann_net[duplicated_row,"annotation"], sep=";")
      ann_net[keeping_row,"mzError"] <- mean(ann_net[[keeping_row,"mzError"]],
                                             ann_net[[duplicated_row,"mzError"]])
    }
    ann_net <- ann_net[-duplicated_ann,]
    ann_net[,"mzError"] <- round(ann_net[,"mzError"], 5)
  }
  
  # check if the direction of the annotations is correct
  adduct_dimer_isotope_anns <- grepl(
    pattern = pattern_ditrimer_adduct_isotope_multiIsotopic, 
    ann_net$annotation, 
    perl = TRUE)
  # adducts|dimers|trimers|isotope : M -> m
  # if (!all(ann_net$msclusterID_source[adduct_dimer_isotope_anns] > 
  #          ann_net$msclusterID_target[adduct_dimer_isotope_anns])) {
  #   print(ann_net[adduct_dimer_isotope_anns,][
  #     ann_net$msclusterID_source[adduct_dimer_isotope_anns] <
  #     ann_net$msclusterID_target[adduct_dimer_isotope_anns], c(1,2,4)])
  if (!all(ms_area_clean$mzConsensus[
    match(ann_net$msclusterID_source[adduct_dimer_isotope_anns], 
          ms_area_clean$msclusterID)] > 
    ms_area_clean$mzConsensus[
      match(ann_net$msclusterID_target[adduct_dimer_isotope_anns], 
            ms_area_clean$msclusterID)])) 
  {
    print(ann_net[adduct_dimer_isotope_anns,][
      ms_area_clean$mzConsensus[
        match(ann_net$msclusterID_source[adduct_dimer_isotope_anns], 
              ms_area_clean$msclusterID)] <
        ms_area_clean$mzConsensus[
          match(ann_net$msclusterID_target[adduct_dimer_isotope_anns], 
                ms_area_clean$msclusterID)], c(1,2,4)])
    stop("Error. A inconsistency was found in the annotations direction for adducts, ",
         "dimers/trimers and isotopes. The above edges of the annotation ",
         "network are in the wrong direction.")
  }
  # fragments|multi charge|H2O and NH3 loss : m -> M
  # added equality for the case of a fragment annotation between the same ion, 
  # which were not reduced in the clean step
  fragments_multicharge_neutralloss <- (!adduct_dimer_isotope_anns)
  if (!all(ms_area_clean$mzConsensus[
    match(ann_net$msclusterID_source[fragments_multicharge_neutralloss], 
          ms_area_clean$msclusterID)] <= 
    ms_area_clean$mzConsensus[
      match(ann_net$msclusterID_target[fragments_multicharge_neutralloss], 
            ms_area_clean$msclusterID)])) 
  {
    print(ann_net[fragments_multicharge_neutralloss,][
      ms_area_clean$mzConsensus[
        match(ann_net$msclusterID_source[fragments_multicharge_neutralloss], 
              ms_area_clean$msclusterID)] >=
        ms_area_clean$mzConsensus[
          match(ann_net$msclusterID_target[fragments_multicharge_neutralloss], 
                ms_area_clean$msclusterID)], c(1,2,4)])
    stop("Error. A inconsistency was found in the annotations direction for fragments, ",
         "multi charges and neutral losses. The above edges of the annotation ",
         "network are in the wrong direction.")
  }
  
  # find missing nodes and add them as selfloops
  singletons_IDs <- ms_area_clean$msclusterID[!(ms_area_clean$msclusterID %in% ann_net$msclusterID_source | 
                                     ms_area_clean$msclusterID %in% ann_net$msclusterID_target)]
  
  ann_net = add_row(ann_net, 
          msclusterID_source=singletons_IDs, 
          msclusterID_target=singletons_IDs, 
          cosine=1.0, 
          annotation="",
          mzError = 0,
          rtError = 0,
          numCommonSamples=NA)
  
  # write the annotation net with selfloops and duplicated links joined
  write.table(ann_net, 
              file.path(output_path, paste0("molecular_networking/", 
                                            output_name,
                                            "_ivamn.selfloop")),
              row.names = FALSE, sep = ",", 
              fileEncoding = "UTF-8")
}


# annoatate the spectra present in the clean table using numerical equivalences
# and chemical rules.
# add the annotations to the clean table and create a molecular network of 
# annotations with the results
annotate_spectra_table_network <- function(output_path,  # path to the last output path folder inside the outs folder
                                           path_batch_metadata, 
                                           modification_file="rules/np3_modifications.csv", 
                                           mz_tol=0.025, 
                                           fragment_tol=0.025, rt_tol=2,
                                           ms2_peaks_intensity_cutoff=15,  # absolute intensity cutoff for fragmented MS2 peaks scaled from 0 to 1000
                                           scale_factor=0.5, 
                                           ion_mode="+",  # + or -
                                           table_limit_size=3000  # set max number of rows to process in a chunck
                                           ) 
{
  print(tibble(variables=c('output_path', 'path_batch_metadata',
                               'modification_file', 'mz_tol', 'fragment_tol',
                               'rt_tol', 'ms2_peaks_intensity_cutoff', 'scale_factor',
                               'ion_mode', 'table_limit_size'), 
                   values=c(output_path, path_batch_metadata, modification_file, 
                            mz_tol, fragment_tol, rt_tol, ms2_peaks_intensity_cutoff,
                            scale_factor, ion_mode, table_limit_size)))
  
  # only add the annotation to the m/z acting as [M+H]+
  # add m/z diff error and rt limits diff error for the current annotation
  annotateSpectra <- function(ann, sim, prec, variant, num_common_samples)
  {
    # compute rt diff error, for blanks use the rtMean, for not blank use the peak boundaries
    if (any((c(as.numeric(prec['rtMin']), as.numeric(variant['rtMin'])) == 0) & 
        (c(as.numeric(prec['rtMax']), as.numeric(variant['rtMax'])) == 1000000))) {
      # there is a blank cluster, use only the rtMean
      variant[['rtError']] <- round(RMSE(variant['rtMean'], prec['rtMean']), 
                                    2)
    } else {
      # no blank cluster, use the peaks boundaries
      variant[['rtError']] <- round(RMSE(unlist(prec[c('rtMin', 'rtMax')]), 
                                       unlist(variant[c('rtMin', 'rtMax')])), 2)
    }
    variant_idx <- as.integer(variant[['idx']])
    # if analog annotation, add to both the prec and the variant ion
    if (grepl(pattern = "analog", ann)) {
      ms_area_clean[variant_idx, "annotation"] <<-
        ifelse(is.na(ms_area_clean[variant_idx, "annotation"]),
               # no annotation, add this one
               annFormat(ann, sim, variant[['mzError']], variant[['rtError']], 
                         prec[['msclusterID']], num_common_samples)
               ,
               # has an annotation, concat them
               paste0(ms_area_clean[variant_idx, "annotation"], 
                      ";",
                      annFormat(ann, sim, variant[['mzError']], variant[['rtError']], 
                                prec[['msclusterID']], 
                                num_common_samples)))[[1]]
      
      return(annFormat(ann, sim, variant[['mzError']], variant[['rtError']], 
                       variant[['msclusterID']], num_common_samples))
    }
    # check if prec is [M+H]+
    # if not assign annotation to variant and retur NA
    # if yes return annotation
    adduct_dimer_isotope_ann <- grepl(
      pattern = pattern_ditrimer_adduct_isotope_multiIsotopic, 
      ann, 
      perl = TRUE)
    # testing regex - dimer/trimer | adduct/isotope | muti charge isotopic pattern 
    # sapply(c(rules_mod$ion, 'fragment', '[M+H+K]2+ isotopic', '[M+3H]3+ isotopic',
    #          '[M+K-H2O]+', '[2M+K-NH3]+', '[M+Na-NH3]+', '[M+H-H2O]+', '[M+H-NH3]+',
    #          '[M+H-NH3-H2O]+', '[M-H-H2O]-', '[M+1]+', '[M+1]-'),
    #        function(x) grepl(
    #   pattern = pattern_ditrimer_adduct_isotope_multiIsotopic,
    #   x, perl = TRUE))
    # adducts|dimers|trimers|isotope : M -> m
    if (adduct_dimer_isotope_ann) {
      if (as.numeric(prec[['mzConsensus']]) <= as.numeric(variant[['mzConsensus']])) {
        # annotate prec
        annFormat(ann, sim, variant[['mzError']], variant[['rtError']], 
                  variant[['msclusterID']], num_common_samples)[[1]]
      } else {
        # annotate variant and return NA
        ms_area_clean[variant_idx, "annotation"] <<-
          ifelse(is.na(ms_area_clean[variant_idx, "annotation"]),
                 # no annotation, add this one
                 annFormat(ann, sim, variant[['mzError']], variant[['rtError']], 
                           prec[['msclusterID']], num_common_samples),
                 # has an annotation, concat them
                 paste0(ms_area_clean[variant_idx, "annotation"], 
                        ";",
                        annFormat(ann, sim, variant[['mzError']], variant[['rtError']], 
                                  prec[['msclusterID']], num_common_samples)))[[1]]
        return(NA)
      }
    } else {
      # fragments|multi charge|H2O and NH3 loss : m -> M
      if (as.numeric(prec[['mzConsensus']]) >= as.numeric(variant[['mzConsensus']])) {
        # annotate prec
        annFormat(ann, sim, variant[['mzError']], variant[['rtError']], 
                  variant[['msclusterID']], num_common_samples)[[1]]
      } else {
        # annotate variant and return NA
        ms_area_clean[variant_idx, "annotation"] <<-
          ifelse(is.na(ms_area_clean[variant_idx, "annotation"]),
                 # no annotation, add this one
                 annFormat(ann, sim, variant[['mzError']], variant[['rtError']], 
                           prec[['msclusterID']], num_common_samples),
                 # has an annotation, concat them
                 paste0(ms_area_clean[variant_idx, "annotation"], 
                        ";",
                        annFormat(ann, sim, variant[['mzError']], variant[['rtError']], 
                                  prec[['msclusterID']], num_common_samples)))[[1]]
        return(NA)
      }
    }
  }
  
  # create the annotation depending on the targets mass
  # only return the annotation where prec is the [M+H]+
  createAnnotation <- function(ann = "[2M+H]+", common_samples, 
                               prec, cluster_variants)
  {
    cluster_annotations <- apply(cluster_variants, 1, function(variant)
    {
      # print(variant)
      # get the number of common samples 
      num_common_samples <- sum(common_samples %in% ms_peaks[as.numeric(variant[['idx']]),4][[1]])
      
      # only annotate as isotopic multicharge or isotope variant if they has at least one sample in common
      # if (num_common_samples == 0 && grepl("(isotopic|\\[M\\+[0-9]+\\])", ann))
      #   return(NA)
      # else
      annotateSpectra(ann, variant[['sim']], prec, variant, num_common_samples)
    })
    # cluster_annotations
    # remove annotation of no common samples in isotopi annotations
    cluster_annotations[!is.na(cluster_annotations)] 
  }
  
  # parsing arguments
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
  batch_metadata <- readMetadataTable(path_batch_metadata)
  output_name <- basename(output_path)
  path_clean_area_count <- file.path(output_path, "count_tables", "clean",
                                     paste0(output_name,"_peak_area_clean.csv"))
  path_clean_spectra_count <- file.path(output_path, "count_tables", "clean",
                                     paste0(output_name,"_spectra_clean.csv"))
  if (!file.exists(path_clean_area_count) || 
      !file.exists(path_clean_spectra_count))
  {
    path_clean_area_count <- file.path(output_path, "count_tables", "clean",
                                       paste0(output_name,"_peak_area_clean_ann.csv"))
    path_clean_spectra_count <- file.path(output_path, "count_tables", "clean",
                                          paste0(output_name,"_spectra_clean_ann.csv"))
    if (!file.exists(path_clean_area_count) || 
        !file.exists(path_clean_spectra_count)) {
      # print in the error message the two possible names of the clean tables - 
      # without and with annotations' info - concatenate both names
      path_clean_area_count <- paste(file.path(output_path, "count_tables", "clean",
                                               paste0(output_name,"_peak_area_clean.csv")), 
                                     path_clean_area_count, sep=" or ")
      path_clean_spectra_count <- paste(file.path(output_path, "count_tables", "clean",
                                                  paste0(output_name,"_spectra_clean.csv")),
                                        path_clean_spectra_count, sep=" or ")
      stop("The clean counts files '", path_clean_area_count,
           "' and '", path_clean_spectra_count,"' do not exists. ",
           "Provide a valid output path to where the csv files with the counts ",
           "of peak area and spectra are located.")
    }
  }
  if (ms2_peaks_intensity_cutoff <= 0 || ms2_peaks_intensity_cutoff >= 1000) 
  {
    stop("The fragmented peaks absolute intensity cutoff value must be between 0 and 1000.",
         "' The provided value equals ", ms2_peaks_intensity_cutoff," is not valid.")
  }
  if (scale_factor == 0) {  # ln scale
    ms2_peaks_intensity_cutoff = 1.0 + log(1.0 + ms2_peaks_intensity_cutoff)
  } else {  # power scale
    ms2_peaks_intensity_cutoff = ms2_peaks_intensity_cutoff ^ scale_factor
  }
  
  # read input
  # fix integer column types that may trigger errors when matching ids or 
  # performing comparisons
  ms_area_clean <- suppressMessages(read_csv(path_clean_area_count, 
                                             guess_max=5000,
                                             col_types = cols(.default="?", 
                                                              msclusterID="i", 
                                                              numSpectra="i")))
  count_columns <- which(names(ms_area_clean) %in% 
                         paste0(batch_metadata$SAMPLE_CODE, "_area"))

  cat("\n** Getting list of fragmented peaks and data collection batches by cluster **\n\n")
  t0 <- Sys.time()
  # get the data collection batch of each sample
  data_collection_batch <- batch_metadata$DATA_COLLECTION_BATCH[
    match(names(ms_area_clean)[count_columns], 
          paste0(batch_metadata$SAMPLE_CODE, "_area"))]
  # get list of fragmented peaks
  # remove baseline peaks bellow an absolute cutoff of 100
  # check if the spectra has at least one fragmented peak mass greater than 
  #   its precursor mass + 5*iso_mass
  ms_peaks <-  Reduce(rbind, lapply(seq_along(ms_area_clean$peaksList), function(i) {
    mzs <- as.numeric(strsplit(ms_area_clean$peaksList[[i]], split = ";")[[1]])
    ints <- as.numeric(strsplit(ms_area_clean$peaksInt[[i]], split = ";")[[1]])
    
    # filter out peaks m/zs around the precursor m/z +- 2 C isotope mass
    # assign an empty list if not peak is left
    precursor_mz_filter <- (mzs < ms_area_clean$mzConsensus[[i]] - iso_mass*2 | 
                              mzs > ms_area_clean$mzConsensus[[i]] + iso_mass*2)
    mzs <- mzs[precursor_mz_filter]
    ints <- ints[precursor_mz_filter]
    
    if (length(mzs) == 0) {
      mzs <- numeric(0)
      ints <- numeric(0)
    } else {
      # ms2 absolute baseline
      int_filter <- (ints >= ms2_peaks_intensity_cutoff)
      if (any(int_filter)) {
        mzs <- mzs[int_filter]
        ints <- ints[int_filter]
      } else {
        # select the top 3 peaks if no peaks passed the filter
        int_filter <- order(ints, decreasing = TRUE)
        ints <- ints[int_filter][1:min(3,length(ints))]
        mzs <- mzs[int_filter][1:length(ints)]
        int_filter <- order(mzs)
        ints <- ints[int_filter]
        mzs <- mzs[int_filter]
      }
    }
    
    collection_batches <- unique(data_collection_batch[
      ms_area_clean[i , count_columns] > 0])
    
    samples <- count_columns[which(ms_area_clean[i , count_columns] > 0)]
    
    list(mzs = mzs, ints = ints, 
         bigger_peaks = ifelse(length(mzs) > 0, 
                               (mzs[length(mzs)] > ms_area_clean$mzConsensus[[i]] + iso_mass*3 + h2o_mass),
                               FALSE), # prevent setting bigger peaks as NULL
         batches = collection_batches,
         samples = samples)
  }), init = NULL)
  ms_area_clean$bigger_peaks <- unlist(ms_peaks[, "bigger_peaks"], recursive=FALSE) # can remove NULL values
  ms_peaks <- ms_peaks[,c(1,2,4,5)] # rm bigger_peaks column
  rm(data_collection_batch)
  
  tf <- Sys.time()
  cat("  * Done in", round(tf-t0, 2), units(tf-t0), "*\n\n")
  
  cat("\n** Annotating variant spectra that are in concurrent MS1 peaks **\n")
  cat("  * Searching for adducts, dimers/trimers, multi charge, neutral loss, fragments, isotopes ions *\n")
  # read modification rules
  rules_mod <- readr::read_csv(modification_file, col_types = "cdiiiid")
  # filter rules by ionization
  if (ion_mode == "+") {
    rules_mod <- rules_mod[rules_mod$charge > 0,]
    ion_mode <- h_e_mass
  } else { # negative
    rules_mod <- rules_mod[rules_mod$charge < 0,]
    rules_mod$charge <- rules_mod$charge * -1 # make charges positive for formula equivalence
    ion_mode <- -h_e_mass
  }
  adduct_max <- max(abs(rules_mod[grepl("^\\[M(?!\\+[0-9]\\])",
                                        x = rules_mod$ion, perl = TRUE) & 
                                  rules_mod$charge == 1 &
                                    rules_mod$neutral_loss == 0,"mzdiff"]))
  rm(modification_file)
  
  ms_area_clean$annotation <- NA
  ms_area_clean$idx <- 1:nrow(ms_area_clean)
  ms_area_clean$sumAreas <- rowSums(ms_area_clean[,count_columns])
  ms_area_clean$multicharge_ion <- 0
  ms_area_clean$isotope_ion <- 0
  
  any_blank <- any(batch_metadata$SAMPLE_TYPE == "blank")
  # remove count columns of ms_area_clean, only keep the used ones
  if (any_blank) {
    ms_area_clean <- ms_area_clean[,c("msclusterID","mzConsensus", "numJoins",
                                      "rtMean","rtMin","rtMax", "bigger_peaks", 
                                      "annotation", "idx", "sumAreas", "peakIds", 
                                      "BLANKS_TOTAL", "multicharge_ion","isotope_ion")]
  } else {
    ms_area_clean <- ms_area_clean[,c("msclusterID","mzConsensus", "numJoins",
                                      "rtMean","rtMin","rtMax", "bigger_peaks", 
                                      "annotation", "idx", "sumAreas", "peakIds", 
                                      "multicharge_ion","isotope_ion")]
  }
  
  # order in the columns and lines: skip header[1] and add -1 to avoid first column
  scans_order <- c(-1, unlist(read.csv(file.path(output_path, 
                                                 "molecular_networking/similarity_tables", 
                                                 paste0("similarity_table_", 
                                                        output_name, "_clean.csv")), 
                                       stringsAsFactors = FALSE, nrows= 1,
                                       strip.white = TRUE, header = FALSE)[-1]))
  nscans <- length(scans_order) - 1 
  col_types <- strrep("d", length(scans_order))
  # read similarity table
  pairwise_sim <- read_csv(file.path(output_path, 
                                     "molecular_networking/similarity_tables", 
                                     paste0("similarity_table_", output_name, "_clean.csv")), 
                           n_max = table_limit_size, skip = 1, 
                           col_names = F, col_types = col_types)
  
  # read ms1 count
  # use this count to add isotopic charge variants > 1:
  ## look in not fragmented mzs if exists a mz in the same peak with mzdiff = 1/charge
  path_ms1_count <- file.path(output_path, "count_tables", 
                              paste0(output_name,"_peak_area_MS1.csv"))
  if (!file.exists(path_ms1_count))
  {
    warning("The MS1 peaks not fragmented m/z's count file '", path_ms1_count,
            "' do not exists. Some annotation heuristics will be disabled.")
    ms_no_spectra_count <- NULL
  } else {
    ms_no_spectra_count <- suppressMessages(read_csv(path_ms1_count, 
                                                     guess_max = 5000,
                                                     col_types = cols(.default="?", 
                                                                      msclusterID="i")))
    
    # add the batches in witch each mz appears
    count_columns_ms1 <- match(paste0(batch_metadata$SAMPLE_CODE, "_area"),
                               names(ms_no_spectra_count))
    data_collection_batch <- batch_metadata$DATA_COLLECTION_BATCH[
      match(names(ms_no_spectra_count)[count_columns_ms1], 
            paste0(batch_metadata$SAMPLE_CODE, "_area"))]
    ms_no_spectra_count$batches <- lapply(seq_len(nrow(ms_no_spectra_count)), function(i) {
      unique(data_collection_batch[
        ms_no_spectra_count[i , count_columns_ms1] > 0])
    })
    # filter columns
    if (any_blank) {
      ms_no_spectra_count <- ms_no_spectra_count[,c("msclusterID", "peakIds", "mzConsensus", "rtMean", 
                                                    "rtMin", "rtMax", "sumAreas","BLANKS_TOTAL", "batches")]
    } else {
      ms_no_spectra_count <- ms_no_spectra_count[,c("msclusterID", "peakIds", "mzConsensus", "rtMean", 
                                                    "rtMin", "rtMax", "sumAreas","batches")]
    }
    
    rm(data_collection_batch, count_columns_ms1)
  }
  rm(path_ms1_count)
  
  progress_anns <- unique(trunc(c(seq(from = 1, to = nscans, by = nscans/25), nscans)))
  cat(paste0("  |", paste0(rep(" ", length(progress_anns)), collapse = ""), "|\n  |", collapse = ""))
  scans_order <- scans_order[-1] # rm heading col
  # Check if the scans from the pairwise table match the scans from the clean table
  # if there is any inconsistence, stop the processing here
  if (anyNA(match(ms_area_clean$msclusterID,scans_order)) || 
      anyNA(match(scans_order,ms_area_clean$msclusterID))) {
    stop("ERROR. An inconsistency was found between the scans list from the pairwise ",
         "similarity table and the msclusterIDs from the clean count table.",
         " There is a total of ", sum(is.na(match(ms_area_clean$msclusterID,scans_order))),
         " mismatched msclusterID in the clean count table, equal to: ", 
         paste(ms_area_clean$msclusterID[which(is.na(match(ms_area_clean$msclusterID,scans_order)))], sep=','),
         " And a total of ", sum(is.na(match(scans_order, ms_area_clean$msclusterID))),
         " mismatched scans in the clean pairwise similarity table, equal to: ", 
         paste(scans_order[which(is.na(match(scans_order,ms_area_clean$msclusterID)))], sep=','), 
         ". Report this error to the development team if there is no logical way to fix it, ", 
         "probable a bug occuried in the clean MGF creation and in the pairwise similarity table computation.")
  }
  # start annotating ionization variant spectra
  for (i in seq_len(nscans))
  {
    #cat("\nStep i: ", i)  # debug in which step the error is occurring
    if (nrow(ms_area_clean) > length(scans_order)) {
      stop("Something went wrong in the annotation step, found an inconsistence in the counts table.",
           " The counts table have more rows than it should. Error row index = ",
           i-1)
    }
    if (i %in% progress_anns) {
      cat("=")
    }
    # check if the cluster was annotated in the count table, if yes set its network row
    if (is.na(ms_area_clean$annotation[[i]])) {
      annotation <- c()
    } else {
      annotation <- strsplit(ms_area_clean$annotation[[i]], ";")[[1]]
    }
    
    # get the current cluster info, scan number and similarity row
    cluster <- ms_area_clean[i,]
    scan_num <- cluster$msclusterID
    
    # get the list of mzs in MS1 and MS2 that are in the same peak of the current cluster
    # remove already annotated clusters -> mscluster ID <= scan_num
    # only annotate blank clusters using other blanks
    if (!any_blank) { # no blanks
      cluster_peak <- ms_area_clean[ms_area_clean$msclusterID > scan_num &
                                      ((ms_area_clean$rtMean >= cluster$rtMin - rt_tol & 
                                          ms_area_clean$rtMean <= cluster$rtMax + rt_tol) |
                                         (cluster$rtMean >= ms_area_clean$rtMin - rt_tol &
                                            cluster$rtMean <= ms_area_clean$rtMax + rt_tol)), ]
      # if not bflag, check peak center and boundaries deviation, remove peaks not aligned
      # the spectra peak center deviation is <= 4 * rt_tol or peak boundaries deviation <= 2*rt_tol
      cluster_peak <- cluster_peak[abs(cluster$rtMean-cluster_peak$rtMean) <= 4*rt_tol |
                                     apply(cluster_peak[,c("rtMin", "rtMax")], 1,
                                           function(x,y)  RMSE(x, y),
                                           y = c(cluster$rtMin, cluster$rtMax)) <= 2*rt_tol,]
      if (!is.null(ms_no_spectra_count)) {
        cluster_peak_ms1 <- ms_no_spectra_count[((ms_no_spectra_count$rtMean >= cluster$rtMin - rt_tol & 
                                                    ms_no_spectra_count$rtMean <= cluster$rtMax + rt_tol) |
                                                   (cluster$rtMean >= ms_no_spectra_count$rtMin - rt_tol &
                                                      cluster$rtMean <= ms_no_spectra_count$rtMax + rt_tol)), ]
        cluster_peak_ms1 <- cluster_peak_ms1[abs(cluster$rtMean-cluster_peak_ms1$rtMean) <= 4*rt_tol |
                                       apply(cluster_peak_ms1[,c("rtMin", "rtMax")], 1,
                                             function(x,y)  RMSE(x, y),
                                             y = c(cluster$rtMin, cluster$rtMax)) <= 2*rt_tol,]
      }
    } else if (cluster$BLANKS_TOTAL > 0) { 
      # blank cluster, filter peak using only the rtMean - to contour baseline peaks
      cluster_peak <- ms_area_clean[ms_area_clean$msclusterID > scan_num & 
                                      ((ms_area_clean$rtMean >= cluster$rtMean - rt_tol*4 & 
                                          ms_area_clean$rtMean <= cluster$rtMean + rt_tol*4) |
                                         ((ms_area_clean$BLANKS_TOTAL == 0 &
                                             (cluster$rtMean >= ms_area_clean$rtMin - rt_tol &
                                                cluster$rtMean <= ms_area_clean$rtMax + rt_tol)) | 
                                            (ms_area_clean$BLANKS_TOTAL > 0 &
                                               (cluster$rtMean >= ms_area_clean$rtMean - rt_tol*4 & 
                                                  cluster$rtMean <= ms_area_clean$rtMean + rt_tol*4)))),]
      if (!is.null(ms_no_spectra_count)) {
        cluster_peak_ms1 <- ms_no_spectra_count[((ms_no_spectra_count$rtMean >= cluster$rtMean - rt_tol*4 & 
                                                    ms_no_spectra_count$rtMean <= cluster$rtMean + rt_tol*4) |
                                                   ((ms_no_spectra_count$BLANKS_TOTAL == 0 &
                                                       (cluster$rtMean >= ms_no_spectra_count$rtMin - rt_tol &
                                                          cluster$rtMean <= ms_no_spectra_count$rtMax + rt_tol)) | 
                                                      (ms_no_spectra_count$BLANKS_TOTAL > 0 &
                                                         (cluster$rtMean >= ms_no_spectra_count$rtMean - rt_tol*4 & 
                                                            cluster$rtMean <= ms_no_spectra_count$rtMean + rt_tol*4)))),]
      }
    } else { 
      # not blank cluster, only use rtMean for filtering blank baseline
      cluster_peak <- ms_area_clean[ms_area_clean$msclusterID > scan_num & 
                                      ((ms_area_clean$rtMean >= cluster$rtMin - rt_tol & 
                                          ms_area_clean$rtMean <= cluster$rtMax + rt_tol) |
                                         (ms_area_clean$BLANKS_TOTAL == 0 &
                                            (cluster$rtMean >= ms_area_clean$rtMin - rt_tol &
                                               cluster$rtMean <= ms_area_clean$rtMax + rt_tol)) | 
                                         (ms_area_clean$BLANKS_TOTAL > 0 &
                                            (cluster$rtMean >= ms_area_clean$rtMean - rt_tol*4 & 
                                               cluster$rtMean <= ms_area_clean$rtMean + rt_tol*4))), ]
      # if not bflag, check peak center and boundaries deviation, remove peaks not aligned
      # the spectra peak center deviation is <= 4 * rt_tol or peak boundaries deviation <= 2*rt_tol
      cluster_peak <- cluster_peak[((cluster_peak$BLANKS_TOTAL > 0) |
                                    (cluster_peak$BLANKS_TOTAL == 0 &
                                    (abs(cluster$rtMean-cluster_peak$rtMean) <= 4*rt_tol |
                                     apply(cluster_peak[,c("rtMin", "rtMax")], 1,
                                           function(x,y)  RMSE(x, y),
                                           y = c(cluster$rtMin, cluster$rtMax)) <= 2*rt_tol))),]
      if (!is.null(ms_no_spectra_count)) {
        cluster_peak_ms1 <- ms_no_spectra_count[(ms_no_spectra_count$rtMean >= cluster$rtMin - rt_tol & 
                                                   ms_no_spectra_count$rtMean <= cluster$rtMax + rt_tol) |
                                                  (ms_no_spectra_count$BLANKS_TOTAL == 0 &
                                                     (cluster$rtMean >= ms_no_spectra_count$rtMin - rt_tol &
                                                        cluster$rtMean <= ms_no_spectra_count$rtMax + rt_tol)) | 
                                                  (ms_no_spectra_count$BLANKS_TOTAL > 0 &
                                                     (cluster$rtMean >= ms_no_spectra_count$rtMean - rt_tol*4 & 
                                                        cluster$rtMean <= ms_no_spectra_count$rtMean + rt_tol*4)), ]
        cluster_peak_ms1 <- cluster_peak_ms1[((cluster_peak_ms1$BLANKS_TOTAL > 0) |
                                        (cluster_peak_ms1$BLANKS_TOTAL == 0 &
                                           (abs(cluster$rtMean-cluster_peak_ms1$rtMean) <= 4*rt_tol |
                                              apply(cluster_peak_ms1[,c("rtMin", "rtMax")], 1,
                                                    function(x,y)  RMSE(x, y),
                                                    y = c(cluster$rtMin, cluster$rtMax)) <= 2*rt_tol))),]
      }
    }
    if (nrow(cluster_peak) == 0) # no concurrent cluster
    {
      next()
    }
    
    # get the similarity values lines
    sim_i <- which_eq(pairwise_sim[[1]], scan_num, 1)
    if (length(sim_i) == 0) # read more lines to find scan_num
    {
      sim_i <- which_eq(scans_order, scan_num, 1)
      pairwise_sim <- read_csv(file.path(output_path, 
                                         "molecular_networking/similarity_tables", 
                                         paste0("similarity_table_", output_name, "_clean.csv")), 
                               n_max = table_limit_size, skip = sim_i, 
                               col_names = F, col_types = col_types)
      sim_i <- which_eq(pairwise_sim[[1]], scan_num, 1)
      if (length(sim_i) == 0) # read more lines
        stop("Error with the new similarity table, wrong scans order.")
    }
    
    if (pairwise_sim[sim_i,1] != scan_num)
      stop("Wrong order between the similarity table and the counts table. Something went wrong in the cleanning step.")
    
    # get the samples idx in which the current cluster appears
    common_samples <- ms_peaks[cluster$idx, 4][[1]]
    
    # remove from the peak the adjacent clusters that do not share a batch with the current cluster
    cluster_peak <- cluster_peak[sapply(ms_peaks[cluster_peak$idx,3], 
                                        function(x, source_batches) {
                                          any(x %in% source_batches)
                                        }, ms_peaks[cluster$idx, 3][[1]]),]
    if (!is.null(ms_no_spectra_count) && nrow(cluster_peak_ms1) > 0) {
      cluster_peak_ms1 <- cluster_peak_ms1[sapply(cluster_peak_ms1$batches, 
                                                  function(x, source_batches) {
                                                    any(x %in% source_batches)
                                                  }, ms_peaks[cluster$idx, 3][[1]]),
                                           ]
    }
    if (nrow(cluster_peak) == 0) # no concurrent cluster
    {
      next()
    }
    
    # get the mass diff between the current cluster and the others in concurrent peaks
    cluster_peak$mzDiff <- round(cluster$mzConsensus - cluster_peak$mzConsensus, 4)
    cluster_peak$mzDiffAbs <- abs(cluster_peak$mzDiff)
    
    if (!is.null(ms_no_spectra_count)) {
      cluster_peak_ms1$mzDiff <- round(cluster$mzConsensus - cluster_peak_ms1$mzConsensus, 4)
      cluster_peak_ms1$mzDiffAbs <- abs(cluster_peak_ms1$mzDiff)
    }
    # get the similary values of the clusters in concurrent peaks
    cluster_peak$sim <- unlist(pairwise_sim[sim_i,-1][match(cluster_peak$msclusterID, scans_order)], 
                               use.names = FALSE)
    
    # keep who has been annotated in the current peak
    cluster_peak$annotated <- FALSE
    
    # set mzError column as NA and mzVariant as false
    cluster_peak$mzError <- NA
    cluster_peak$mzVariant <- FALSE
    
    ### MULTI CHARGE ###
    # Multi charge annotation -> check if current cluster is a multi charge
    # store if there is a multi charge isotopic variant mz in the same peak
    isotopic_charge <- c(FALSE, FALSE) # +2, +3
    # check if there is a multi charge isotopic variant mz in the same peak
    # check isotopic distribution of double charge ions in the MS2 peak
    cluster_peak$isotopic_charge2 <- FALSE
    cluster_peak$isotopic_charge2 <- ((abs(cluster_peak$mzDiffAbs - h_e_mass/2) <= fragment_tol) &
                                      (cluster_peak$sumAreas <= 2/3*cluster$sumAreas))
    if (any(cluster_peak$isotopic_charge2))
    {
      # double charge, isotopic distribution
      isotopic_charge[1] <- TRUE
    } else if (!is.null(ms_no_spectra_count) && 
               any(abs(cluster_peak_ms1$mzDiffAbs - h_e_mass/2) <= fragment_tol &
                   cluster_peak_ms1$mzDiff < 0 &   # cluster_peak_ms1 contain the isotopic variant
                   cluster_peak_ms1$sumAreas <= 2/3*cluster$sumAreas)) {
      # check isotopic variant of the lowest mass in MS1 concurrent peak
      isotopic_charge[1] <- TRUE
    }
    # check isotopic distribution of triple charge ions in the cluster_peak MS2 peak
    cluster_peak$isotopic_charge3 <- FALSE
    cluster_peak$isotopic_charge3 <- (abs(cluster_peak$mzDiffAbs - h_e_mass/3) <= fragment_tol &
                                      cluster_peak$sumAreas <= 2/3*cluster$sumAreas)
    if (any(cluster_peak$isotopic_charge3))
    {
      # triple charge isotopic ion
      isotopic_charge[2] <- TRUE
    } else if (!is.null(ms_no_spectra_count) && 
               any(abs(cluster_peak_ms1$mzDiffAbs - iso_mass/3) <= fragment_tol &
                   cluster_peak_ms1$mzDiff < 0 &   # cluster_peak_ms1 contain the isotopic variant
                   cluster_peak_ms1$sumAreas <= 2/3*cluster$sumAreas)) {
      # check isotopic variant of the lowest mass in MS1 concurrent peak
      isotopic_charge[2] <- TRUE
    }
    # set valid isotopic charges
    isotopic_charge <- c(2,3)[isotopic_charge]
    
    if (length(isotopic_charge) > 0) {
      # check for the monocharge ion of the current cluster and possible multicharge ion [M+X]Y+
      for (multiCharge in which(rules_mod$charge %in% isotopic_charge)) {
        # ((M-ion_mode)/Y + X) - m <= mz_tol
        # ajusted formula to deal with M+1 = [M+H] and not [M]
        cluster_peak$mzError <- abs(((cluster_peak$mzConsensus - ion_mode)/rules_mod$charge[[multiCharge]] + 
                                      rules_mod$mzdiff[[multiCharge]]) - cluster$mzConsensus)
        cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol)
        if (any(cluster_peak$mzVariant))
        {
          # if multi charge isotopic pattern is present on MS2, annotate here
          if (rules_mod$charge[[multiCharge]] == 2 && any(cluster_peak$isotopic_charge2)) {
            cluster_peak[cluster_peak$isotopic_charge2,'mzError'] <- abs(cluster_peak[cluster_peak$isotopic_charge2,'mzDiffAbs'] - h_e_mass/2)
            annotation <- c(annotation, createAnnotation(paste(rules_mod$ion[[multiCharge]],"isotopic"), 
                                                       common_samples,
                                                       cluster, 
                                                       cluster_peak[cluster_peak$isotopic_charge2,]))
            # add multi charge flag to current cluster and the isotopic variants
            ms_area_clean[cluster$idx, "multicharge_ion"] <- 2
            ms_area_clean[cluster_peak$idx[cluster_peak$isotopic_charge2], "multicharge_ion"] <- 1
            # remove isotopic_charge2 flag to only create the annotation once
            cluster_peak$annotated[cluster_peak$isotopic_charge2] <- TRUE
            cluster_peak$isotopic_charge2 <- FALSE
          } else if (rules_mod$charge[[multiCharge]] == 3 && any(cluster_peak$isotopic_charge3)) {
            cluster_peak[cluster_peak$isotopic_charge3,'mzError'] <- abs(cluster_peak[cluster_peak$isotopic_charge3,'mzDiffAbs'] - h_e_mass/3)
            annotation <- c(annotation, createAnnotation(paste(rules_mod$ion[[multiCharge]],"isotopic"), 
                                                         common_samples,
                                                         cluster, cluster_peak[cluster_peak$isotopic_charge3,]))
            # add multi charge flag to current cluster and the isotopic variants
            ms_area_clean[cluster$idx, "multicharge_ion"] <- 2
            ms_area_clean[cluster_peak$idx[cluster_peak$isotopic_charge3], "multicharge_ion"] <- 1
            cluster_peak$annotated[cluster_peak$isotopic_charge3] <- TRUE
            cluster_peak$isotopic_charge3 <- FALSE
          } else {
            # add multi charge flag to current cluster if no isotopic variant in MS2, only in MS1
            ms_area_clean[cluster$idx, "multicharge_ion"] <- 2
          }
          # multi charge annotate
          for (j in which(cluster_peak$mzVariant))
          {
            # current cluster always have mz smaller than cluster_peak j
            # annotate cluster as multi charge ion of the cluster_peak j mono charged
            cluster_peak$annotated[j] <- cluster_peak$annotated[j] | cluster_peak$mzVariant[j]
            annotation <- c(annotation, createAnnotation(rules_mod$ion[[multiCharge]], common_samples,
                                                         cluster, cluster_peak[j,]))
          }
        }
      }
    }
    # if there is an multi charge isotopic variant in the current MS2 peak but no monocharge was found in MS2
    # check if the monocharge exists in the MS1 peak, if yes annotate the isotopic variants
    # if not annotate with the default multicharge [M+2H]2+
    if (any(cluster_peak$isotopic_charge2))
    {
      # multicharge isotopic variant present in MS2
      cluster_peak[cluster_peak$isotopic_charge2,'mzError'] <- abs(cluster_peak[cluster_peak$isotopic_charge2,'mzDiffAbs'] - h_e_mass/2)
      # search for the monocharge in MS1
      for (multiCharge in which(rules_mod$charge == 2)) {
        # ((M-ion_mode)/Y + X) - m <= mz_tol
        # ajusted formula to deal with M+1 = [M+H] and not [M]
        if (!is.null(ms_no_spectra_count) && any(abs(((cluster_peak_ms1$mzConsensus - ion_mode)/rules_mod$charge[[multiCharge]] + 
                  rules_mod$mzdiff[[multiCharge]]) - cluster$mzConsensus) <= mz_tol))
        {
          # add multi charge isotopic annotation of the variant present in the MS1
          annotation <- c(annotation, createAnnotation(paste(rules_mod$ion[[multiCharge]],"isotopic"), 
                                                         common_samples,
                                                         cluster, cluster_peak[cluster_peak$isotopic_charge2,]))
          # add multi charge flag to current cluster and the isotopic variants
          ms_area_clean[cluster$idx, "multicharge_ion"] <- 2
          ms_area_clean[cluster_peak$idx[cluster_peak$isotopic_charge2], "multicharge_ion"] <- 1
          cluster_peak$annotated[cluster_peak$isotopic_charge2] <- TRUE
          cluster_peak$isotopic_charge2 <- FALSE
          break()
        }
      }
      if (any(cluster_peak$isotopic_charge2)) {
        # add multi charge isotopic annotation here when the mono charge 
        # is not present in the MS2 neither in the MS1
        annotation <- c(annotation, createAnnotation(ifelse(ion_mode>0, 
                                                            "[M+2H]2+ isotopic",
                                                            "[M-2H]2- isotopic"), 
                                                     common_samples,
                                                     cluster, cluster_peak[cluster_peak$isotopic_charge2,]))
        # add multi charge flag to current cluster and the isotopic variants
        ms_area_clean[cluster$idx, "multicharge_ion"] <- 1
        ms_area_clean[cluster_peak$idx[cluster_peak$isotopic_charge2], "multicharge_ion"] <- 1
        cluster_peak$annotated[cluster_peak$isotopic_charge2] <- TRUE
        cluster_peak$isotopic_charge2 <- FALSE
      }
    } else if (ms_area_clean[[cluster$idx, "multicharge_ion"]] == 0 && 2 %in% isotopic_charge) {
      # isotopic present in MS1 and monocharge not present in MS2 nor MS1
      # set flag
      ms_area_clean[cluster$idx, "multicharge_ion"] <- 1
    }
    if (any(cluster_peak$isotopic_charge3))
    {
      # multicharge isotopic variant present in MS2
      cluster_peak[cluster_peak$isotopic_charge3,'mzError'] <- abs(cluster_peak[cluster_peak$isotopic_charge3,'mzDiffAbs'] - h_e_mass/3)
      # search for the monocharge in MS1
      for (multiCharge in which(rules_mod$charge == 3)) {
        # ((M-ion_mode)/Y + X) - m <= mz_tol
        # ajusted formula to deal with M+1 = [M+H] and not [M]
        if (!is.null(ms_no_spectra_count) && any(abs(((cluster_peak_ms1$mzConsensus - ion_mode)/rules_mod$charge[[multiCharge]] + 
                     rules_mod$mzdiff[[multiCharge]]) - cluster$mzConsensus) <= mz_tol))
        {
          # add multi charge isotopic annotation of the variant present in the MS1
          annotation <- c(annotation, createAnnotation(paste(rules_mod$ion[[multiCharge]],"isotopic"), 
                                                       common_samples,
                                                       cluster, cluster_peak[cluster_peak$isotopic_charge3,]))
          # add multi charge flag to current cluster and the isotopic variants
          ms_area_clean[cluster$idx, "multicharge_ion"] <- 2
          ms_area_clean[cluster_peak$idx[cluster_peak$isotopic_charge3], "multicharge_ion"] <- 1
          cluster_peak$annotated[cluster_peak$isotopic_charge3] <- TRUE
          cluster_peak$isotopic_charge3 <- FALSE
          break()
        }
      }
      if (any(cluster_peak$isotopic_charge3)) {
        # add multi charge isotopic annotation here when the mono charge 
        # is not present in the MS2 neither in the MS1
        annotation <- c(annotation, createAnnotation(ifelse(ion_mode>0,
                                                            "[M+3H]3+ isotopic",
                                                            "[M-3H]3- isotopic"), 
                                                     common_samples,
                                                     cluster, cluster_peak[cluster_peak$isotopic_charge3,]))
        # add multi charge flag to current cluster and the isotopic variants
        ms_area_clean[cluster$idx, "multicharge_ion"] <- 1
        ms_area_clean[cluster_peak$idx[cluster_peak$isotopic_charge3], "multicharge_ion"] <- 1
        cluster_peak$annotated[cluster_peak$isotopic_charge3] <- TRUE
        cluster_peak$isotopic_charge3 <- FALSE
      }
    } else if (ms_area_clean[[cluster$idx, "multicharge_ion"]] == 0 && 3 %in% isotopic_charge) {
      # isotopic present in MS1 and monocharge not present in MS2 nor MS1
      # set flag
      ms_area_clean[cluster$idx, "multicharge_ion"] <- 1
    }
    # update multicharge_ion flag value
    cluster$multicharge_ion <- ms_area_clean[[cluster$idx, "multicharge_ion"]]
    # test annotation
    # for (multiCharge in which(rules_mod$charge > 1)) {
    #   cat(rules_mod$ion[[multiCharge]], " mzvariant = ",
    #       (853.33089/rules_mod$charge[[multiCharge]] + rules_mod$mzdiff[[multiCharge]]),
    #       "\n")
    # }
    
    ### DIMERS/TRIMERS ###
    # check for dimers ions [2M+X]+ of the cluster
    # (2*(M-ion_mode) + X) - m <= mz_tol
    if (!cluster$bigger_peaks &&
        any(cluster_peak$mzConsensus > (cluster$mzConsensus - ion_mode)*2 + h_e_mass - h2o_mass - 2*mz_tol)) {
      # if cluster do not have bigger peaks it could be a monomer
      for (dimer in which(grepl("^\\[(2|3)M",x = rules_mod$ion))) {
        # check if it is a dimer or trimer (2 or 3) if not annotatted as multi charge
        mono_var <- as.numeric(regmatches(rules_mod$ion[[dimer]], 
                                          regexpr("(?<=^\\[)(2|3)(?=M)", 
                                                  rules_mod$ion[[dimer]], perl = T)))
        cluster_peak$mzError <- abs(((cluster$mzConsensus - ion_mode) * mono_var + 
                                         rules_mod$mzdiff[[dimer]]) - cluster_peak$mzConsensus)
        cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol) & !cluster_peak$annotated
        for (j in which(cluster_peak$mzVariant))
        {
          # check if there is no peak between the monomer mass + H2O_mass + methanol_mass and the dimer/trimer mass or 
          # no peak greater then the dimer mass in the dimer/trimer ion, or they are at least half similar
          # it is expected that the dimer/trimer fragments in the monomer, with no peaks after it
          if (!any((ms_peaks[[cluster_peak$idx[[j]]]] > cluster$mzConsensus + rules_mod$mzdiff[[dimer]] + h2o_mass + 32.02566 & 
                    ms_peaks[[cluster_peak$idx[[j]]]] < cluster_peak$mzConsensus[[j]] - iso_mass*5 - h2o_mass) | 
                   ms_peaks[[cluster_peak$idx[[j]]]] > cluster_peak$mzConsensus[[j]] + iso_mass*5) ||
              cluster_peak$sim[[j]] >= rules_mod$sim_cutoff[[dimer]]) 
          {
            cluster_peak$annotated[j] <- cluster_peak$annotated[j] | cluster_peak$mzVariant[j]
            annotation <- c(annotation, createAnnotation(rules_mod$ion[[dimer]], common_samples,
                                                         cluster, cluster_peak[j,]))
          } 
        }
        # check for the dimer with neutral loss of water when specified
        if (rules_mod$neutral_loss_h2o[[dimer]] == 1)
        {
          cluster_peak$mzError <- abs(((cluster$mzConsensus - ion_mode)*mono_var + 
                                         rules_mod$mzdiff[[dimer]] - h2o_mass) - cluster_peak$mzConsensus)
          cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol) & !cluster_peak$annotated
          for (j in which(cluster_peak$mzVariant))
          {
            # check if there is no peak between the monomer mass + H2O_mass + methanol_mass and the dimer/trimer mass or 
            # no peak greater then the dimer mass in the dimer/trimer ion, or they are at least half similar
            # it is expected that the dimer/trimer fragments in the monomer, with no peaks after it
            if (!any((ms_peaks[[cluster_peak$idx[[j]]]] > cluster$mzConsensus + rules_mod$mzdiff[[dimer]] + h2o_mass + 32.02566 & 
                      ms_peaks[[cluster_peak$idx[[j]]]] < cluster_peak$mzConsensus[[j]] - iso_mass*5 - h2o_mass) | 
                     ms_peaks[[cluster_peak$idx[[j]]]] > cluster_peak$mzConsensus[[j]] + iso_mass*5) ||
                cluster_peak$sim[[j]] >= rules_mod$sim_cutoff[[dimer]]) 
            {
              cluster_peak$annotated[j] <- cluster_peak$annotated[j] | cluster_peak$mzVariant[j]
              annotation <- c(annotation, createAnnotation(sub(ifelse(ion_mode>0, "]\\+", "]-"),
                                                               ifelse(ion_mode>0, "-H2O]+", "-H2O]-"), 
                                                               rules_mod$ion[[dimer]]), 
                                                           common_samples,
                                                           cluster, cluster_peak[j,]))
            } 
          }
        }
        # check for the dimer with neutral loss of amonia when specified
        if (rules_mod$neutral_loss_nh3[[dimer]] == 1)
        {
          cluster_peak$mzError <- abs(((cluster$mzConsensus - ion_mode)*mono_var + 
                                        rules_mod$mzdiff[[dimer]] - nh3_mass) - cluster_peak$mzConsensus )
          cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol) & !cluster_peak$annotated
          for (j in which(cluster_peak$mzVariant))
          {
            # check if there is no peak between the monomer mass + H2O_mass + methanol_mass and the dimer/trimer mass or 
            # no peak greater then the dimer mass in the dimer/trimer ion, or they are at least half similar
            # it is expected that the dimer/trimer fragments in the monomer, with no peaks after it
            if (!any((ms_peaks[[cluster_peak$idx[[j]]]] > cluster$mzConsensus + rules_mod$mzdiff[[dimer]] + h2o_mass + 32.02566 & 
                      ms_peaks[[cluster_peak$idx[[j]]]] < cluster_peak$mzConsensus[[j]] - iso_mass*5 - h2o_mass) | 
                     ms_peaks[[cluster_peak$idx[[j]]]] > cluster_peak$mzConsensus[[j]] + iso_mass*5) ||
                cluster_peak$sim[[j]] >= rules_mod$sim_cutoff[[dimer]]) 
            {
              cluster_peak$annotated[j] <- cluster_peak$annotated[j] | cluster_peak$mzVariant[j]
              annotation <- c(annotation, createAnnotation(sub(ifelse(ion_mode>0, "]\\+", "]-"),
                                                               ifelse(ion_mode>0, "-NH3]+", "-NH3]-"),
                                                               rules_mod$ion[[dimer]]), common_samples,
                                                           cluster, cluster_peak[j,]))
            } 
          }
        }
      }
    }
    # test   
    # for (dimer in which(grepl("^\\[(2|3)M",x = rules_mod$ion))) {
    #   cat(rules_mod$ion[[dimer]], " mzvariant = ",
    #     (853.33089*as.numeric(regmatches(rules_mod$ion[[dimer]], 
    #                                      regexpr("(?<=^\\[)(2|3)(?=M)", 
    #                                              rules_mod$ion[[dimer]], perl = T))) + 
    #        rules_mod$mzdiff[[dimer]]),
    #     "\n")
    # }
    
    ### FRAGMENTS ###
    # check for fragments using the fragment tolerance
    # do not annotate frags from blank samples
    # do not annotate as fragment if variant was annotated as dimer/trimer or multicharge
    ## the problem is that the monomer can be in the list of fragmented peaks of the dimer/trimer, the same could happen with multicharges
    cluster_peak$mzVariant <- !cluster_peak$annotated
    if (any_blank) {
      cluster_peak$mzVariant <- cluster_peak$mzVariant & cluster_peak$BLANKS_TOTAL == 0
    }
    for (j in which(cluster_peak$mzVariant)) {
      # cluster has bigger mz then neighbor j
      if (cluster_peak$mzDiff[[j]] > 0) {
        # look for the mass of cluster j in the fragmented peaks of the current cluster
        if (any(abs(ms_peaks[[i]] - cluster_peak$mzConsensus[[j]]) <= fragment_tol)) {
          k <- cluster_peak$idx[[j]]
          cluster_peak$sim[[j]] <- max(round(normDotProductTrim(ms_peaks[[i]], ms_peaks[[i,2]], ms_peaks[[k]], ms_peaks[[k,2]],
                                                                fragment_tol, cluster_peak$mzConsensus[[j]]-iso_mass*2), 2),
                                       cluster_peak$sim[[j]], na.rm = TRUE)
          # check if their trimmed similarity is above a fixed threshold of 0.2
          if (cluster_peak$sim[[j]] >= 0.2) {
            cluster_peak$annotated[j] <- TRUE
            cluster_peak$mzError[[j]] <- abs(ms_peaks[[i]] - cluster_peak$mzConsensus[[j]])[
              abs(ms_peaks[[i]] - cluster_peak$mzConsensus[[j]]) <= fragment_tol][[1]]
            annotation <- c(annotation, createAnnotation("fragment", common_samples,
                                                         cluster, cluster_peak[j,]))
          }
        } 
      } else {
        # cluster has smaller mz then neighbor j
        # look for the mass of the current cluster in the fragmented peaks of the cluster j
        k <- cluster_peak$idx[[j]]
        if (any(abs(ms_peaks[[k]] - cluster$mzConsensus) <= fragment_tol)) {
          cluster_peak$sim[[j]] <- max(round(normDotProductTrim(ms_peaks[[i]], ms_peaks[[i,2]], ms_peaks[[k]], ms_peaks[[k,2]],
                                                                fragment_tol, cluster$mzConsensus-iso_mass*2), 2),
                                       cluster_peak$sim[[j]], na.rm = TRUE)
          # check if their trimmed similarity is above a fixed threshold of 0.2
          if (cluster_peak$sim[[j]] >= 0.2) {
            cluster_peak$annotated[j] <- TRUE
            cluster_peak$mzError[[j]] <- abs(ms_peaks[[k]] - cluster$mzConsensus)[
              abs(ms_peaks[[k]] - cluster$mzConsensus) <= fragment_tol][[1]]
            annotation <- c(annotation, createAnnotation("fragment", common_samples,
                                                         cluster, cluster_peak[j,]))
          }
        }
      }
    }
    
    ### ADDUCTS ###
    # dist(adducts_mzs, method = "manhattan")
    # check mass diff for adducts
    # [M+X]+ or [M+X-H2O]+
    # | |M - ion_mode - m| - X| <= mz_tol
    if (min(cluster_peak$mzDiffAbs) < adduct_max + mz_tol)
    {
      # grepl adducts and not isotopes, remove neutral losses annotations as well
      for (adduct in which(grepl("^\\[M(?!\\+[0-9]\\])",
                                 x = rules_mod$ion, perl = TRUE) & 
                           rules_mod$charge == 1 &
                           rules_mod$neutral_loss == 0))
      {
        cluster_peak$mzError <- abs(cluster_peak$mzDiffAbs + ion_mode - rules_mod$mzdiff[[adduct]])
        cluster_peak$mzVariant <- ((cluster_peak$mzError <= mz_tol) & (cluster_peak$sim >= rules_mod$sim_cutoff[[adduct]]))
        if (any(cluster_peak$mzVariant))
        {
          # cat(adduct)
          cluster_peak$annotated[cluster_peak$mzVariant] <- cluster_peak$annotated[cluster_peak$mzVariant] | cluster_peak$mzVariant[cluster_peak$mzVariant]
          annotation <- c(annotation, createAnnotation(rules_mod$ion[[adduct]], 
                                                       common_samples,
                                                       cluster, cluster_peak[cluster_peak$mzVariant,]))
        } 
        # also check for adduct with neutral loss of water
        if (rules_mod$neutral_loss_h2o[[adduct]] == 1)
        {
          cluster_peak$mzError <- abs(cluster_peak$mzDiffAbs + ion_mode - 
                                        (rules_mod$mzdiff[[adduct]] - h2o_mass))
          cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol) & (cluster_peak$sim >= rules_mod$sim_cutoff[[adduct]])
          if (any(cluster_peak$mzVariant)) {
            cat(adduct)
            cluster_peak$annotated[cluster_peak$mzVariant] <- cluster_peak$annotated[cluster_peak$mzVariant] | cluster_peak$mzVariant[cluster_peak$mzVariant]
            annotation <- c(annotation, createAnnotation(sub(ifelse(ion_mode>0, "]\\+", "]-"),
                                                             ifelse(ion_mode>0, "-H2O]+", "-H2O]-"), rules_mod$ion[[adduct]]), 
                                                         common_samples,
                                                         cluster, cluster_peak[cluster_peak$mzVariant,]))
          }
        }
        # also check for adduct with neutral loss of amonia
        if (rules_mod$neutral_loss_nh3[[adduct]] == 1)
        {
          cluster_peak$mzError <- abs(cluster_peak$mzDiffAbs + ion_mode - 
                                        (rules_mod$mzdiff[[adduct]] - nh3_mass))
          cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol) & (cluster_peak$sim >= rules_mod$sim_cutoff[[adduct]])
          if (any(cluster_peak$mzVariant)) {
            cat(adduct)
            cluster_peak$annotated[cluster_peak$mzVariant] <- cluster_peak$annotated[cluster_peak$mzVariant] | cluster_peak$mzVariant[cluster_peak$mzVariant]
            annotation <- c(annotation, createAnnotation(sub(ifelse(ion_mode>0, "]\\+", "]-"),
                                                             ifelse(ion_mode>0, "-NH3]+", "-NH3]-"), rules_mod$ion[[adduct]]), 
                                                         common_samples,
                                                         cluster, cluster_peak[cluster_peak$mzVariant,]))
          }
        }
        
      }
    }
    # test   
    # for (adduct in which(grepl("^\\[M",x = rules_mod$ion) & rules_mod$charge == 1)) {
    #   cat(rules_mod$ion[[adduct]], " mzvariant = ",
    #     (853.33089 +
    #        rules_mod$mzdiff[[adduct]]),
    #     "\n")
    # }
    
    ### NEUTRAL LOSS ###
    # similarity cutoff of 0.2 by default is needed for this annotation
    # e.g. neutral loss of water = [M+ion_mode-H2O]+
    for (neutral_loss in which(rules_mod$neutral_loss == 1))
    {
      cluster_peak$mzError <- abs(cluster_peak$mzDiffAbs - rules_mod$mzdiff[[neutral_loss]])
      cluster_peak$mzVariant <- ((cluster_peak$mzError <= mz_tol) & (cluster_peak$sim >= rules_mod$sim_cutoff[[neutral_loss]]))
      if (any(cluster_peak$mzVariant))
      {
        cluster_peak$annotated[cluster_peak$mzVariant] <- cluster_peak$annotated[cluster_peak$mzVariant] | cluster_peak$mzVariant[cluster_peak$mzVariant]
        annotation <- c(annotation, createAnnotation(rules_mod$ion[[neutral_loss]],
                                                     common_samples,
                                                     cluster, cluster_peak[cluster_peak$mzVariant,]))
      }
    }
    # water neutral loss
    # cluster_peak$mzError <- abs(cluster_peak$mzDiffAbs - h2o_mass)
    # cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol) & (cluster_peak$sim >= 0.2)
    # if (any(cluster_peak$mzVariant))
    # {
    #   cluster_peak$annotated[cluster_peak$mzVariant] <- cluster_peak$annotated[cluster_peak$mzVariant] | cluster_peak$mzVariant[cluster_peak$mzVariant]
    #   annotation <- c(annotation, createAnnotation(ifelse(ion_mode>0, 
    #                                                       "[M+H-H2O]+",
    #                                                       "[M-H-H2O]-"),
    #                                                common_samples,
    #                                                cluster, cluster_peak[cluster_peak$mzVariant,]))
    # }
    # # amonia neutral loss
    # cluster_peak$mzError <- abs(cluster_peak$mzDiffAbs - nh3_mass)
    # cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol)  & (cluster_peak$sim >= 0.2)
    # if (any(cluster_peak$mzVariant))
    # {
    #   cluster_peak$annotated[cluster_peak$mzVariant] <- cluster_peak$annotated[cluster_peak$mzVariant] | cluster_peak$mzVariant[cluster_peak$mzVariant]
    #   annotation <- c(annotation, createAnnotation(ifelse(ion_mode>0, 
    #                                                       "[M+H-NH3]+",
    #                                                       "[M-H-NH3]-"),
    #                                                common_samples,
    #                                                cluster, cluster_peak[cluster_peak$mzVariant,]))
    # }
    # # water and amonia neutral loss
    # cluster_peak$mzError <- abs(cluster_peak$mzDiffAbs - nh3_mass - h2o_mass)
    # cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol) & (cluster_peak$sim >= 0.2)
    # if (any(cluster_peak$mzVariant))
    # {
    #   cluster_peak$annotated[cluster_peak$mzVariant] <- cluster_peak$annotated[cluster_peak$mzVariant] | cluster_peak$mzVariant[cluster_peak$mzVariant]
    #   annotation <- c(annotation, createAnnotation(ifelse(ion_mode>0, 
    #                                                       "[M+H-NH3-H2O]+",
    #                                                       "[M-H-NH3-H2O]-"),
    #                                                common_samples,
    #                                                cluster, cluster_peak[cluster_peak$mzVariant,]))
    # }
    
    ### ISOTOPE ###
    # check for a single C13 isotopo in the variant and high similarity value
    # e.g. [M+1]+ and [M+2]+
    for (isotope in which(grepl("^\\[M\\+[0-9]\\]", x = rules_mod$ion))) 
    {
      iso_var <- as.numeric(regmatches(rules_mod$ion[[isotope]], 
                                        regexpr("(?<=^\\[M\\+)[0-9](?=\\])", 
                                                rules_mod$ion[[isotope]], perl = T)))
      cluster_peak$mzError <- abs(cluster_peak$mzDiffAbs - iso_mass * iso_var)
      cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol & 
                                   cluster_peak$sim >= rules_mod$sim_cutoff[[isotope]])
      for (j in which(cluster_peak$mzVariant))
      {
        # add the isotope ion flag
        if (cluster$mzConsensus < cluster_peak$mzConsensus[j]) {
          ms_area_clean$isotope_ion[cluster_peak$idx[j]] <- 1
        } else {
          ms_area_clean$isotope_ion[cluster$idx] <- 1
        }
        cluster_peak$annotated[j] <- cluster_peak$annotated[j] | cluster_peak$mzVariant[j]
        annotation <- c(annotation, createAnnotation(rules_mod$ion[[isotope]], 
                                                     common_samples,
                                                     cluster, cluster_peak[j,]))
      }
    }
    # cluster_peak$mzError <- abs(cluster_peak$mzDiffAbs - iso_mass)
    # cluster_peak$mzVariant <- (cluster_peak$mzError <= mz_tol & 
    #                            cluster_peak$sim >= 0.75)
    # for (j in which(cluster_peak$mzVariant))
    # {
    #   # add the isotope ion flag
    #   if (cluster$mzConsensus < cluster_peak$mzConsensus[j]) {
    #     ms_area_clean$isotope_ion[cluster_peak$idx[j]] <- 1
    #   } else {
    #     ms_area_clean$isotope_ion[cluster$idx] <- 1
    #   }
    #   cluster_peak$annotated[j] <- cluster_peak$annotated[j] | cluster_peak$mzVariant[j]
    #   annotation <- c(annotation, createAnnotation(ifelse(ion_mode>0, "[M+1]+", 
    #                                                       "[M+1]-"), 
    #                                                common_samples,
    #                                                cluster, cluster_peak[j,]))
    #   
    # }
    
    ### ANALOG ###
    # if spectra have similarity above a high threshold of 0.7 and no annotation,
    # annotate them as similar with their mz diff
    for (j in which(!cluster_peak$annotated & cluster_peak$sim >= 0.7))
    {
      cluster_peak$mzError[[j]] <- cluster_peak$mzDiffAbs[j]
      annotation <- c(annotation, 
                      createAnnotation("analog", 
                                       common_samples, 
                                       cluster, cluster_peak[j,]))
    }
    
    if (length(annotation) > 0)
      ms_area_clean$annotation[[i]] <- paste(annotation, collapse = ";")
  }
  cat("|\n")
  t0 <- Sys.time()
  cat("  * Done annotating clusters in", 
          round(t0-tf, 2), units(t0-tf), "*\n\n")
  ms_area_clean$bigger_peaks <- NULL
  ms_area_clean$idx <- NULL
  ms_area_clean$sumAreas <- NULL
  rm(ms_no_spectra_count)
  
  cat("\n** Creating the IVAMN and separating the annotations by types of ionization variants **\n")
  # create the annotation network file heading
  write_csv(data.frame("msclusterID_source","msclusterID_target","cosine","annotation",
                       "mzError", "rtError", "numCommonSamples"), 
            path = file.path(output_path, paste0("molecular_networking/", 
                                                 output_name,
                                                 "_ivamn.selfloop")),
            col_names = FALSE)
  # save annotation net before writing tables to save annotations separated by type
  annotations_cols <- bind_rows(lapply(seq_len(nrow(ms_area_clean)), 
                                       save_annotation_net, 
                                       ms_area_clean[,c("msclusterID","annotation")], 
                                       scans_order, output_path,
                                       output_name))
  # add the multichar ion flag to the annotations cols
  annotations_cols <- bind_cols(ms_area_clean[c("multicharge_ion", "isotope_ion")],
                                annotations_cols)
  # write annotation singletons
  write_annotation_net_singletons_join_duplicated(output_path, output_name, 
                                                  ms_area_clean[,c("msclusterID", 
                                                                   "mzConsensus")])
  rm(ms_area_clean)
  tf <- Sys.time()
  cat("  * Done in", round(tf-t0, 2), units(tf-t0), "*\n\n")
  
  # read clean counts and bind the annotation columns
  ms_area_clean <- suppressMessages(read_csv(path_clean_area_count, 
                                             guess_max = 5000,
                                             col_types = cols(.default="?", 
                                                              msclusterID="i", 
                                                              numSpectra="i")))
  # remove previous annotations columns if present
  ms_area_clean$isotope_ion <- NULL
  ms_area_clean$multicharge_ion <- NULL
  ms_area_clean$adducts <- NULL
  ms_area_clean$isotopes <- NULL
  ms_area_clean$dimers <- NULL
  ms_area_clean$multiCharges <- NULL
  ms_area_clean$fragments <- NULL
  ms_area_clean$analogs <- NULL
  ms_area_clean <- bind_cols(ms_area_clean, annotations_cols)
  write_csv(ms_area_clean, path = file.path(output_path, "count_tables", "clean",
                                            paste0(output_name, "_peak_area_clean_ann.csv")))
  rm(ms_area_clean)
  
  ms_spectra_count <- suppressMessages(read_csv(path_clean_spectra_count, 
                                                guess_max = 5000,
                                                col_types = cols(.default="?", 
                                                                 msclusterID="i", 
                                                                 numSpectra="i")))
  # remove previous annotations columns if present
  ms_spectra_count$isotope_ion <- NULL
  ms_spectra_count$multicharge_ion <- NULL
  ms_spectra_count$adducts <- NULL
  ms_spectra_count$isotopes <- NULL
  ms_spectra_count$dimers <- NULL
  ms_spectra_count$multiCharges <- NULL
  ms_spectra_count$fragments <- NULL
  ms_spectra_count$analogs <- NULL
  ms_spectra_count <- bind_cols(ms_spectra_count, annotations_cols)
  write_csv(ms_spectra_count, path = file.path(output_path, "count_tables", "clean", 
                                               paste0(output_name, "_spectra_clean_ann.csv")))
  
  # remove intermediary clean count files when exists
  if (file.exists(file.path(output_path, "count_tables", "clean", paste0(output_name, "_peak_area_clean.csv")))) {
    invisible(file.remove(file.path(output_path, "count_tables", "clean", paste0(output_name, "_peak_area_clean.csv")),
                        file.path(output_path, "count_tables", "clean", paste0(output_name, "_spectra_clean.csv"))))
  }
}
