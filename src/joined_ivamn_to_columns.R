##
# Step 7 - Get the annotations present in the IVAMN and convert them to 
# the annotation columns to be added to the clean counts tables
# use reverse engineering from the original annotation flow
##

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("One argument is necessary to write to the clean quantifications the joined annotations present in the joined IVAMN:\n", 
       " 1 - output_path: Path to the final output folder of the joined job result, inside the 'outs' dir.\n",
       call.=FALSE)
} else {
  output_path <- file.path(args[[1]])
  if (!dir.exists(output_path)) {
    stop("The provided output_path folder '", output_path,
         "' does not exists. ",
         "Provide a valid output path, inside the 'outs' folder where the final ",
         "result of the joined job is located.")
  }
}

cat("Loading packages readr, dplyr...\n")
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

# format the annotations present in one edge of the ivamn, separate them by
# ionization variant type and return a list with each type filled with the 
# respective formatted annotation
annFormat <- function(ann, sim, mzError, rtError, variant_ID, num_common_samples) {
  #mzError <- round(as.numeric(mzError), 3)
  #variant_ID <- as.integer(trimws(variant_ID))
  # if more than one annotation between the source and target node, split them
  # and categorize each one of them
  ann_split <- strsplit(ann, split=";")[[1]]
  
  annotations <- list(adducts = NA,
                      isotopes = NA,
                      dimers = NA,
                      multiCharges = NA,
                      fragments = NA)
  # retrieve and format each annotation and set it in the list of annotations
  for (i in seq_along(ann_split)) {
    # format the annotation
    ann_format <- paste0(ann_split[[i]], " (sim ", sim, " - mzE ", 
                         mzError, " - rtE ",rtError,")", "[", variant_ID, "]", 
                         "{", num_common_samples, "}")
    
    # fill annotations by type of variant, concatenate them if needed
    if (grepl(pattern = "\\[(2M|3M)", ann_split[[i]])) {
      annotations$dimers <- ifelse(is.na(annotations$dimers), 
                                   ann_format, 
                                   paste(annotations$dimers, ann_format, sep = ";"))
    } else if (grepl(pattern = "\\](2\\+|3\\+)", ann_split[[i]])) {
      annotations$multiCharges <- ifelse(is.na(annotations$multiCharges), 
                                         ann_format, 
                                         paste(annotations$multiCharges, ann_format, sep = ";")) 
    } else if (grepl(pattern = "fragment", ann_split[[i]])) {
      annotations$fragments <- ifelse(is.na(annotations$fragments), 
                                      ann_format, 
                                      paste(annotations$fragments, ann_format, sep = ";"))
    } else if (grepl(pattern = "\\[M\\+[0-9]+\\]", ann_split[[i]])) {
      annotations$isotopes <- ifelse(is.na(annotations$isotopes), 
                                     ann_format, 
                                     paste(annotations$isotopes, ann_format, sep = ";"))
    } else { # adduct and neutral loss
      annotations$adducts <- ifelse(is.na(annotations$adducts), 
                                    ann_format, 
                                    paste(annotations$adducts, ann_format, sep = ";")) 
    }
  }
  return(annotations)
}

convert_ivamn_annotation_columns <- function(output_path) {
  output_path <- normalizePath(output_path)
  output_name <- basename(output_path)
  path_clean_area_count <- file.path(output_path, "count_tables", "clean",
                                     paste0(output_name,"_peak_area_clean_ann.csv"))
  path_clean_spectra_count <- file.path(output_path, "count_tables", "clean",
                                        paste0(output_name,"_spectra_clean_ann.csv"))
  if (!file.exists(path_clean_area_count) || 
      !file.exists(path_clean_spectra_count)) {
    stop("The clean counts table files '", path_clean_area_count,
         "' and '", path_clean_spectra_count,"' do not exists. ",
         "Provide a valid output path to where the csv files with the counts ",
         "of peak area and spectra are located.")
  }
  path_ivamn <- file.path(output_path, "molecular_networking", 
                          paste0(output_name,"_ivamn.selfloop"))
  if (!file.exists(path_ivamn)) {
    stop("The joined IVAMN file '", path_ivamn,
         "' inside the provided output_path does not exists. ",
         "Provide a valid output path to where the IVAMN with the ionization ",
         "variants annotations is located.")
  }
  
  ivamn <- suppressMessages(read_csv(path_ivamn, guess_max = 5000,
                                     col_types = cols(.default="?", 
                                                      msclusterID_source="i", 
                                                      msclusterID_target="i")))
  # remove selfloops - NA annotations
  ivamn <- ivamn[!is.na(ivamn$annotation),]
  
  # if no annotations are present, skip conversion
  # set the columns as empty cols
  if (nrow(ivamn) > 0) {
    # for each target, extract the income edges annotations
    # and convert them to the annotations columns
    ivamn <- ivamn %>% arrange(msclusterID_target)
    
    # create a table to store the formated annotation by ionization variant type
    # the analogs are not added here
    mspectra_annotation <- data.frame(list(msclusterID = unique(ivamn$msclusterID_target), 
                                           adducts = NA,
                                           isotopes = NA,
                                           dimers = NA,
                                           multiCharges = NA,
                                           fragments = NA))
    # for each msclusterID in the target list, retrive its annotations from
    # income edges and format them to be added to the counts table
    for (i in seq_along(mspectra_annotation$msclusterID)) {
      # retrieve the income edges of msID, separate its annotations by ;
      # for each annotation convert them to the annFormat
      # and then, separate them by type
      msID <- mspectra_annotation$msclusterID[[i]]
      edges_msID <- ivamn[ivamn$msclusterID_target == msID,]
      for (j in seq_len(nrow(edges_msID))) {
          edge_annFormat <- annFormat(edges_msID$annotation[[j]], sim=edges_msID$cosine[[j]], 
                    mzError=edges_msID$mzError[[j]], rtError=edges_msID$rtError[[j]], 
                    variant_ID=edges_msID$msclusterID_source[[j]], 
                    num_common_samples=edges_msID$numCommonSamples[[j]])
          # store the valid annotations by ionizaion variant type
          for (ion_type in names(edge_annFormat)[!is.na(edge_annFormat)]) {
            mspectra_annotation[i, ion_type] <- ifelse(
              is.na(mspectra_annotation[i, ion_type]), 
              edge_annFormat[[ion_type]],
              paste(mspectra_annotation[i, ion_type], 
                    edge_annFormat[[ion_type]], sep = ";"))
          }
      }
    }
  } 
  # append the annotation columns to the clean quantification files
  # count by peak area
  ms_clean_count <- suppressMessages(read_csv(path_clean_area_count, 
                                             guess_max=5000,
                                             col_types = cols(.default="?", 
                                                              msclusterID="i", 
                                                              numSpectra="i")))
  # remove annotation columns if already exists (previous run)
  ms_clean_count[,c("adducts","isotopes","dimers","multiCharges","fragments")] <- NULL
  if (nrow(ivamn) > 0) {
    ms_clean_count <- left_join(ms_clean_count, mspectra_annotation,
                                by="msclusterID")
  } else {
    ms_clean_count[,c("adducts","isotopes","dimers","multiCharges","fragments")] <- NA
  }
  write_csv(ms_clean_count, path = path_clean_area_count)
  # count by number of spectra
  ms_clean_count <- suppressMessages(read_csv(path_clean_spectra_count, 
                                              guess_max=5000,
                                              col_types = cols(.default="?", 
                                                               msclusterID="i", 
                                                               numSpectra="i")))
  # remove annotation columns if already exists (previous run)
  ms_clean_count[,c("adducts","isotopes","dimers","multiCharges","fragments")] <- NULL
  if (nrow(ivamn) > 0) {
    ms_clean_count <- left_join(ms_clean_count, mspectra_annotation, 
                                by="msclusterID")
  } else {
    ms_clean_count[,c("adducts","isotopes","dimers","multiCharges","fragments")] <- NA
  }
  write_csv(ms_clean_count, path = path_clean_spectra_count)
  
  # TODO avaliate if analogs should be retrieved again
  # remove intermediary clean count files when exists
  if (file.exists(file.path(output_path, "count_tables", "clean", paste0(output_name, "_peak_area_clean.csv")))) {
    invisible(file.remove(file.path(output_path, "count_tables", "clean", paste0(output_name, "_peak_area_clean.csv")),
                          file.path(output_path, "count_tables", "clean", paste0(output_name, "_spectra_clean.csv"))))
  }
}

# call the function to convert the annotations present in the ivamn to columns 
# in the clean counts tables
convert_ivamn_annotation_columns(output_path = output_path)
