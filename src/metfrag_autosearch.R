cat("Loading packages metfRag, dlls...\n")
suppressPackageStartupMessages(library(metfRag))
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
Rcpp::sourceCpp(file.path(script_path(), 'read_mgf_peak_list_R.cpp'))

try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf)
{
  y <- try({setTimeLimit(cpu, elapsed); expr}, silent = TRUE) 
  if (inherits(y, "try-error")) NULL else y 
}

# output_path <- "/home/crisfbazz/Documents/CNPEM/ms_entity_annotation/Data_Collections/Bra346/Bra_test_metfrag/outs/Bra_test_metfrag"
mode_run <- 1
ppm_tol <- 5
fragment_tol <- 0.003
corr_method <- "spearman"
scale_factor <- 0.5

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Six arguments must be supplied to create the specification files for MSCLuster:\n", 
       " 1 - Path to the job result folder inside the 'outs' directory;\n",
       " 2 - The correlation method used;\n",
       " 3 - The ion adduct type: 1 = [M+H]+ or -1 = [M-H]-;\n",
       " 4 - mz tolerance in PPM;\n",
       " 5 - mz tolerance for fragment peaks in Daltons;\n",
       " 6 - The spectra fragmented peaks scaling method: sqrt, ln or NULL.",
       call.=FALSE)
} else {
  output_path <- args[[1]]
  if (!dir.exists(output_path))
  {
    stop("The job output folder '", output_path, 
         "' do not exists. Provide a valid path to where the job final result is located.")
  }
  output_name <- basename(output_path)
  
  path_mgf_sample <- file.path(output_path, "mgf", paste0(output_name, "_all.mgf"))
  if (!file.exists(path_mgf_sample))
  {
    stop("The mgf file '", path_mgf_sample,
         "' do not exists. Provide a valid output path to where the clustered mgf exists.")
  }
  
  path_area_count <- file.path(output_path, "count_tables", "clean", 
                               paste0(output_name, "_peak_area_clean_ann.csv"))
  if (!file.exists(path_area_count))
  {
    path_area_count <- file.path(output_path, "count_tables",
                                 paste0(output_name, "_peak_area.csv"))
    if (!file.exists(path_area_count))
    {
      stop("The count file '", path_area_count,
           "' do not exists. Provide a valid output path to where the CSV file with the count ",
           "of peak area is located.")
    }
  }
  
  corr_method <- args[[2]]
  path_count_corr <- file.path(output_path, "count_tables", "clean", 
                               paste0(output_name, 
                                      paste0("_peak_area_clean_annotated_corr_", 
                                             corr_method, ".csv")))
  if (!file.exists(path_count_corr))
  {
    path_count_corr <- file.path(output_path, "count_tables",
                                 paste0(output_name, 
                                        paste0("_peak_area_corr_", 
                                               corr_method, ".csv")))
    if (!file.exists(path_count_corr))
    {
      warning("The count file '", path_count_corr,
              "' with the correlation scores do not exists. Filter disabled, ",
              "identification will be perfomed in all not blank m/zs", 
              call. = FALSE)
      path_count_corr <- "FALSE"
    }
  }
  
  mode_run <- as.numeric(args[[3]])
  if (!(mode_run %in% c(1, -1)))
  {
    warning("Wrong precursor ion mode provided '", mode_run,
            "'. It should be one of 1 = [M+H]+ or -1 = [M-H]-.",
            "Ion adduct type set to 1.", call. = FALSE)
    mode_run <- 1
  }
  
  if (length(args) == 6) {
    ppm_tol <- as.numeric(args[[4]])
    fragment_tol <- as.numeric(args[[5]])
    scale_factor <- as.numeric(args[[6]])
    if (is.na(scale_factor) || is.null(scale_factor) || scale_factor < 0)
    {
      warning("Invalid scale factor provided, it must be a numeric value greater or equal to 0: \n", 
              "  - 0 : ln natural logarithm of the intensities\n",
              "  - 1 : no scale\n",
              "  - x : x > 0 pow of the intensities to x. (e.g. x = 0.5 square root)\n",
              "Scale factor equals 0.5 (sqrt) will be selected by default.", call. = FALSE)
      scale_factor <- 0.5
    }
  }
}

# first define the settings object
settings_object <- list()

# set polarity mode and accepted adducts
# Ion adduct type (1 = [M+H]⁺, 0 = [M]⁺⁻, 18 = [M+NH4]⁺, 23 = [M+Na]⁺, 39 = [M+K]⁺,
# -1 = [M-H]⁻, 35 = [M+Cl]⁻, 45 = [M+HCOO]⁻, 59 = [M+CH3COO]⁻)  
if (mode_run == 1)
{
  settings_object[["IsPositiveIonMode"]] <- TRUE
  settings_object[["PrecursorIonMode"]]  <- mode_run # [M+H]+ | [M]⁺⁻
  mz_diff_adduct_search <- -1.007284 * mode_run
  # mz_diff_adduct_search <- c("[M]+" = 0, "[M+H]+" = 1)
} else {
  settings_object[["IsPositiveIonMode"]] <- FALSE
  settings_object[["PrecursorIonMode"]]  <- mode_run # [M-H]-
  mz_diff_adduct_search <- 1.007284 
}

ms_area_count <- read.csv(path_area_count, stringsAsFactors = FALSE,
                          comment.char = "", strip.white = TRUE)

# read mgf and get scan index
msn_sample <- readMgfPeaksList(path_mgf_sample, bin_size = fragment_tol, 
                               trim_mz = FALSE, scale_factor = scale_factor,
                               join_isotopic_peaks = 0)

# check if number of spectra in the count file match with mgf
if (length(msn_sample$SCANS) < length(unique(c(ms_area_count$msclusterID, 
                                             unlist(strsplit(ms_area_count$joinedIDs[
                                               !is.na(ms_area_count$joinedIDs)], ";"))))))
  stop("The number of spectra in the mgf do not match with the number of ", 
       "spectra in the count files. Check if the files are correct and retry.")

# create the identifications output directory if it does not exists
output_path <- file.path(output_path, "identifications")
if (!dir.exists(output_path))
  dir.create(output_path, mode = "777")

if ("BLANKS_TOTAL" %in% names(ms_area_count)) {
  ms_area_count <- ms_area_count[ms_area_count$BLANKS_TOTAL == 0,]
}
if ("BEDS_TOTAL" %in% names(ms_area_count)) {
  ms_area_count <- ms_area_count[ms_area_count$BEDS_TOTAL == 0,]
}

# select max 20 not flagged blank m/zs and with any corr value >= 0.5 & with any corr >= max(corr) - 0.3
if (path_count_corr != "FALSE")
{
  count_corr <- suppressMessages(readr::read_csv(path_count_corr, skip = 1, guess_max = 2000))
  if ("BLANKS_TOTAL" %in% names(count_corr)) {
    count_corr <- count_corr[count_corr$BLANKS_TOTAL == 0,]
  }
  if ("BEDS_TOTAL" %in% names(count_corr)) {
    count_corr <- count_corr[count_corr$BEDS_TOTAL == 0,]
  }
  count_corr <- apply(count_corr[, startsWith(names(count_corr), "COR_")], 1, 
                      function(x)
                      {
                        x[which(x == "CTE")] <- 0
                        x <- suppressWarnings(as.numeric(x)) # suppress NA introduced warnings
                        x[is.na(x)] <- -1
                        
                        max(x, na.rm = TRUE)
                      })
  # valid_search_idx <- which(ms_area_count$BLANKS_TOTAL == 0 & 
  #                           rowSums(count_corr >= corr_cutoff, 
  #                                   na.rm = TRUE) > 0)
  valid_search_idx <- which(count_corr >= max(count_corr, na.rm = T)*0.7)
  
  # get the top 30 correlated
  valid_search_idx <- valid_search_idx[order(count_corr[valid_search_idx], decreasing = TRUE)]
  if (length(valid_search_idx) > 30)
    valid_search_idx <- valid_search_idx[1:30]
  
  rm(count_corr)
} else {
  valid_search_idx <- which(ms_area_count$msclusterID >= 0) # all mzs
}

# select the most intense spectra, the one that have more peaks
ms_area_count <- ms_area_count[valid_search_idx, c("msclusterID", "joinedIDs")]
ms_area_count$index <- sapply(1:nrow(ms_area_count), function(x) {
  ids <- as.numeric(unique(c(ms_area_count[x,1], strsplit(ms_area_count[x,2], ";")[[1]])))
  ids <- match(ids[!is.na(ids)], msn_sample$SCANS)
  ids <- ids[which.max(sapply(msn_sample$MZS[ids], length))]
  ids
})

settings_object[["DatabaseSearchRelativeMassDeviation"]]    <- ppm_tol/2
settings_object[["FragmentPeakMatchAbsoluteMassDeviation"]] <- fragment_tol
settings_object[["FragmentPeakMatchRelativeMassDeviation"]] <- ppm_tol
settings_object[["MetFragDatabaseType"]] <- "PubChem"
settings_object[["FragmenterScore"]] <- "1.0"

# define pre and post process filters 
# filter non-connected compounds (e.g. salts)
settings_object[["MetFragPreProcessingCandidateFilter"]] <- c("UnconnectedCompoundFilter", "IsotopeFilter")
# filter stereoisomers by comparing first part of compounds' InChIKeys
# only the best-scored candidate remains in the result list
settings_object[["MetFragPostProcessingCandidateFilter"]] <- c("InChIKeyFilter")

n_search <- nrow(ms_area_count)
cat("** Running MetFrag Identification with PubChem Database **\n")
ti <- Sys.time()
# write the result table header
write.table(list("msclusterID", "ionMode_PubChem", "fragmenterScore_PubChem", "numExplPeaks_PubChem",
                 "identifier_PubChem", "SMILES_PubChem", "molecularFormula_PubChem"), 
            file = file.path(output_path, "metfrag_pubchem_best_results.csv"), 
            row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
# search the best correlated idxs
search_best_results <- lapply(seq_len(n_search), function(x)
{
  gc()
  cat("  \n* Searching msclusterID", ms_area_count$msclusterID[[x]], 
          "(",x,"/", n_search, ") *\n")
  i <- ms_area_count$index[[x]]
  # spectrumq <- msn_sample[[ms_area_count$index[[x]]]]
  
  # TODO: check all mz_diff_adduct_search mzs, add adduct column and mz column
  # set precursor mass
  settings_object[["NeutralPrecursorMass"]] <- msn_sample$PREC_MZ[[i]] + mz_diff_adduct_search
  cat("\nNeutral precursor mass search: ", settings_object[["NeutralPrecursorMass"]], "\n")
  # define the peaklist as 2-dimensional matrix
  # get the top 30 peaks
  settings_object[["PeakList"]] <- matrix(c(msn_sample$MZS[[i]], 
                                            msn_sample$INTS[[i]]), 
                                          ncol=2)
  
  # run MetFrag
  t1 <- Sys.time()
  # max time to 20 min
  scored_candidates <- try_with_time_limit(run.metfrag(settings_object), 1200, 1200)
  t2 <- Sys.time()
  cat("\n Done seaching in PubChem: ", round(t2-t1, 2), " ", units(t2-t1), "\n")
  
  if (is.null(scored_candidates)) 
  {
    message("\n!!! Connection with MetFrag server took too long and was interrupted. =X !!!")
    search_result <- list(msclusterID = ms_area_count$msclusterID[[x]],
                          ionMode_PubChem = as.character(mode_run),
                          fragmenterScore_PubChem = "timeout",
                          numExplPeaks_PubChem = "timeout",
                          identifier_PubChem = "timeout",
                          SMILES_PubChem = "timeout",
                          molecularFormula_PubChem = "timeout")
  } else if (nrow(scored_candidates) == 0) {
    search_result <- list(msclusterID = ms_area_count$msclusterID[[x]],
                          ionMode_PubChem = as.character(mode_run),
                          fragmenterScore_PubChem = "no results",
                          numExplPeaks_PubChem = "no results",
                          identifier_PubChem = "no results",
                          SMILES_PubChem = "no results",
                          molecularFormula_PubChem = "no results")
  } else {
    # get top 30 candidates
    scored_candidates <- scored_candidates[1:min(nrow(scored_candidates), 30),]
    # save all results to csv
    write.csv(scored_candidates, 
              file.path(output_path, 
                        paste0(ms_area_count$msclusterID[[x]], "_", mode_run,  "_PubChem.csv")),
              row.names = FALSE)
    
    # return the best 10 results
    scored_candidates$numExplPeaks <- paste0(scored_candidates[["NoExplPeaks"]], "/", 
                                             scored_candidates[["NumberPeaksUsed"]])
    scored_candidates <- apply(scored_candidates[1:min(nrow(scored_candidates), 10), ], 
                               2, paste, collapse=";")
    search_result <- list(msclusterID = ms_area_count$msclusterID[[x]],
                          ionMode_PubChem = as.character(mode_run),
                          fragmenterScore_PubChem = as.character(scored_candidates[["FragmenterScore"]]),
                          numExplPeaks_PubChem = as.character(scored_candidates[["numExplPeaks"]]),
                          identifier_PubChem = as.character(scored_candidates[["Identifier"]]),
                          SMILES_PubChem = as.character(gsub(";", ",", scored_candidates[["SMILES"]])),
                          molecularFormula_PubChem = as.character(scored_candidates[["MolecularFormula"]]))
  }
  
  write.table(search_result, 
              file = file.path(output_path, "metfrag_pubchem_best_results.csv"), 
              row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE, quote = TRUE)
})
tf <- Sys.time()
cat("\n** Done seaching for ", n_search, 
        " PubChem identifications in ", 
        round(tf-ti, 2), " ", units(tf-ti), " **\n")