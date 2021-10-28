## ----load-libs, message = FALSE--------------------------------------------
message("Loading packages metfRag, MSnbase...")
suppressPackageStartupMessages(library(metfRag))
suppressPackageStartupMessages(library(MSnbase))

.Last <- function() {
  close(filestdin)
  rm(msn_sample)
  
  # remove duplicated searchs and keep the last one
  if (file.exists(file.path(output_path, "metfrag_pubchem_best_results.csv")))
  {
    metfrag_best_results <- read.csv(file.path(output_path, "metfrag_pubchem_best_results.csv"),
                                     stringsAsFactors = FALSE, comment.char = "", 
                                     strip.white = TRUE)
    
    metfrag_best_results <- metfrag_best_results[!duplicated(
      metfrag_best_results[,c("msclusterID", "ionMode_PubChem")], fromLast = T),]
    write.csv(metfrag_best_results, file.path(output_path, "metfrag_pubchem_best_results.csv"),
              row.names = FALSE)
  }
}

# path_mgf <- "../Data_Collections/Bra346/Bra_RT1_2_small/outs/Bra_RT1_2_small/mgf/Bra_RT1_2_small_all_clean.mgf"
# output_path <- "../Data_Collections/Bra346/Bra_RT1_2_small/outs/Bra_RT1_2_small/identifications/"

filestdin <- file("stdin")

run_metfrag <- function(i, all_run = FALSE)
{
  gc()
  x <- match(i, scan_indexes)
  message("  \n* Searching msclusterID ", i, ifelse(all_run, 
                                                    paste0(" (",x,"/", n_search, ") *"),
                                                    "*"))
  spectrumq <- msn_sample[[x]]
  
  # TODO: check all mz_diff_adduct_search mzs, add adduct column and mz column
  # set precursor mass
  settings_object[["NeutralPrecursorMass"]] <- spectrumq@precursorMz + mz_diff_adduct_search
  
  # define the peaklist as 2-dimensional matrix
  # get the top 30 peaks
  peaks_idx <- sort(order(spectrumq@intensity, decreasing = TRUE)[1:30])
  settings_object[["PeakList"]] <- matrix(c(spectrumq@mz[peaks_idx], 
                                            spectrumq@intensity[peaks_idx]), 
                                          ncol=2)
  # print settings of single search
  if (!all_run)
    print(settings_object)
  # run MetFrag
  t1 <- Sys.time()
  scored_candidates <- tryCatch(run.metfrag(settings_object), error = function(e) {
    print(e)
    return(NULL)
  })
  t2 <- Sys.time()
  message("\n Done seaching in PubChem: ", round(t2-t1, 2), " ", units(t2-t1), "")
  
  if (!exists("scored_candidates") || is.null(scored_candidates) || is.na(scored_candidates)) 
  {
    search_result <- list(msclusterID = i,
                          ionMode_PubChem = as.character(ion_mode),
                          fragmenterScore_PubChem = NA,
                          numExplPeaks_PubChem = NA,
                          identifier_PubChem = NA,
                          SMILES_PubChem = NA,
                          molecularFormula_PubChem = NA)
  } else if (nrow(scored_candidates) == 0) {
    search_result <- list(msclusterID = i,
                          ionMode_PubChem = as.character(ion_mode),
                          fragmenterScore_PubChem = "no results",
                          numExplPeaks_PubChem = "no results",
                          identifier_PubChem = "no results",
                          SMILES_PubChem = "no results",
                          molecularFormula_PubChem = "no results")
  } else {
    # save all results to csv
    write.csv(scored_candidates, 
              file.path(output_path, 
                        paste0(i, "_", ion_mode, "_PubChem.csv")),
              row.names = FALSE)
    
    # return the best 2 results
    scored_candidates$numExplPeaks <- paste0(scored_candidates[["NoExplPeaks"]], "/", 
                                             scored_candidates[["NumberPeaksUsed"]])
    scored_candidates <- apply(scored_candidates[1:2, ], 2, paste, collapse=";")
    search_result <- list(msclusterID = i,
                          ionMode_PubChem = as.character(ion_mode),
                          fragmenterScore_PubChem = as.character(scored_candidates[["FragmenterScore"]]),
                          numExplPeaks_PubChem = as.character(scored_candidates[["numExplPeaks"]]),
                          identifier_PubChem = as.character(scored_candidates[["Identifier"]]),
                          SMILES_PubChem = as.character(scored_candidates[["SMILES"]]),
                          molecularFormula_PubChem = as.character(scored_candidates[["MolecularFormula"]]))
  }
  
  write.table(search_result, 
              file = file.path(output_path, "metfrag_pubchem_best_results.csv"), 
              row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
}

read_id <- function(msg="Enter id to search: ")
{
  cat(msg)
  n <- gsub(" ", "", readLines(filestdin,1))
  
  if (n %in% c("q", "0", "all"))
    return(n)
  
  # check if entry was x
  if (!grepl("^[0-9]+$",n))
  {
    message("\n* invalid input format *\n")
    return(read_id(msg))
  } 
  
  if (!(n %in% scan_indexes)) {
    message("\n* invalid ID *\n")
    return(read_id(msg))
  }
  
  return(n)
}

read_ion_mode <- function(msg="Enter the ion mode: ")
{
  cat(msg)
  n <- gsub(" ", "", readLines(filestdin,1))
  
  if (n == "")
  {
    cat("[M]+-\n")
    return("M")
  }
  
  if (n %in% c("q", "0", "1", "-1", "18", "23", "35", "39"))
  {
    return(n)
  } else {
    message("\n* invalid input *\n")
    return(read_ion_mode(msg))
  }
}

read_tol <- function(msg="Enter frag tol: ", default_tol)
{
  cat(msg)
  n <- gsub(" ", "", readLines(filestdin,1))
  
  if (n == "")
  {
    cat(default_tol, "\n")
    return(default_tol)
  }
  
  if (n %in% c("q"))
    return(n)
  
  if (!grepl("^[0-9]+\\.?[0-9]*$",n))
  {
    message("\n* invalid input format *\n")
    return(read_tol(msg))
  }
  
  return(as.numeric(n))
}

read_countfile <- function(msg="Enter filename: ")
{
  cat(msg)
  filename <- readLines(filestdin,1)
  
  if (filename %in% c("q", "0"))
    return(filename)
  
  if (filename == "")
  {
    cat("all\n")
    return("all")
  }
    
  if (!file.exists(filename))
  {
    message("\n* invalid path *\n")
    return(read_countfile(msg))
  }
  
  return(filename)
}


# ids prompt before first plot
id_prompt <- "\n|====  | \nEnter the spectra SCANS number (msclusterID) to be searched in PubChem or 'all' to search all the missing SCANS identifications in the provided output folder:
q : quit
0 : back
all : search all the missing SCANS in the output folder
<X> : a numeric equals the SCANS number (msclusterID) present in the MGF file
spectra ID = "

tol_prompt <- "\n|==     | \nEnter the m/z tolerance in Daltons for fragment peak matching:
q : quit
ENTER : set default equals '0.025'
<X> : a numeric equals the fragment peak tolerance
tolerance = "

rel_prompt <- "\n|=    | \nEnter the m/z tolerance in PPM for the precursor mass matching:
q : quit
ENTER : set default equals '3'
<X> : a numeric equals the PPM tolerance
tolerance = "

ion_mode_prompt <- "\n|===   | \nEnter the precursor ion mode equals to the ion adduct type to be searched. The adduct mass will be subtracted or added to the precurson ion mass:
q : quit
0 : back
ENTER : set default equals [M]+- for neutral precursor mass (no adduct)
1 : [M+H]+
18 : [M+NH4]+
23 : [M+Na]+
39 : [M+K]+
-1 : [M-H]-
35 : [M+Cl]-
ion mode = "

file_prompt <- "\n|======| \nEnter the count file path corresponding to the provided MGF to exclude the search of all spectra that appeared in at least one blank sample (BLANKS_TOTAL > 0):
q : quit
0 : back
ENTER : search all missing SCANS
count file path = "

# defaults
mz_tol <- 0.05
ppm_tol <- 15

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Two arguments must be supplied to run the MetFrag search in the PubChem DB:\n",
       " 1 - Path to the output folder where the search results will be saved. ",
       "Preferably the diretory 'identifications' inside the final clustering result folder ",
       "inside the 'outs' directory;\n",
       " 2 - Path to the MGF containing the clustering result.", call. = FALSE)
} else {
  path_mgf <- file.path(args[[1]])
  if (!file.exists(path_mgf))
  {
    stop("The MGF file '", path_mgf,
         "' do not exists. Provide a valid path to where the MGF file is located.")
  }
  path_mgf <- normalizePath(path_mgf)
  
  output_path <- file.path(args[[2]])
  if (!dir.exists(output_path))
  {
    stop("The output path '", output_path,
         "' do not exists. Provide a valid path to where the search results will be saved.")
  }
  output_path <- normalizePath(output_path)
}

message("Reading MGF...\n")

## Create a phenodata data.frame
pd <- data.frame(sample_name = "search",
                 sample_group = 1,
                 class = "sample",
                 stringsAsFactors = FALSE) 

# read mgf and get scan index
msn_sample <- readMgfData(path_mgf, new("NAnnotatedDataFrame", pd), 
                          verbose = TRUE)
scan_indexes <- scanIndex(msn_sample)

# first define the settings object
settings_object <- list()
settings_object[["MetFragDatabaseType"]] <- "PubChem"
settings_object[["FragmenterScore"]] <- "1.0"
# define pre and post process filters 
# filter non-connected compounds (e.g. salts)
# TODO check IsotopeFilter
settings_object[["MetFragPreProcessingCandidateFilter"]] <- c("UnconnectedCompoundFilter")
# filter stereoisomers by comparing first part of compounds' InChIKeys
# only the best-scored candidate remains in the result list
settings_object[["MetFragPostProcessingCandidateFilter"]] <- c("InChIKeyFilter")

repeat
{
  ppm_tol <- read_tol(rel_prompt, default_tol = 3)
  
  if (ppm_tol == "q")
    q("no")
  
  settings_object[["DatabaseSearchRelativeMassDeviation"]]    <- ppm_tol
  settings_object[["FragmentPeakMatchRelativeMassDeviation"]] <- ppm_tol
  
  repeat
  {
    mz_tol <- read_tol(tol_prompt, default_tol = 0.025)
    
    if (mz_tol == "0")
      break()
    else if (mz_tol == "q")
      q("no")
    
    settings_object[["FragmentPeakMatchAbsoluteMassDeviation"]] <- mz_tol
    
    repeat
    {
      ion_mode <- read_ion_mode(ion_mode_prompt)
      
      if (ion_mode == "0")
        break()
      else if (ion_mode == "q")
        q("no")
      else if (ion_mode == "M")
        ion_mode <- 0
      else
        ion_mode <- as.numeric(ion_mode)
      
      # set polarity mode and accepted adducts
      # Ion adduct type (1 = [M+H]⁺, 0 = [M]⁺⁻, 18 = [M+NH4]⁺, 23 = [M+Na]⁺, 39 = [M+K]⁺,
      # -1 = [M-H]⁻, 35 = [M+Cl]⁻)  
      if (ion_mode %in% c(0, 1, 18, 23, 39))
      {
        settings_object[["IsPositiveIonMode"]] <- TRUE
        settings_object[["PrecursorIonMode"]]  <- ion_mode 
        if (ion_mode == 18)
          mz_diff_adduct_search <- -18.03384
        else if (ion_mode == 23)
          mz_diff_adduct_search <- -22.98922
        else if (ion_mode == 39)
          mz_diff_adduct_search <- -38.96316
        else # ion_mode == 1 or 0
          mz_diff_adduct_search <- -1.007284 * ion_mode
        # mz_diff_adduct_search <- c("[M]+" = 0, "[M+H]+" = 1)
      } else {
        settings_object[["IsPositiveIonMode"]] <- FALSE
        settings_object[["PrecursorIonMode"]]  <- ion_mode # [M-H]-
        if (ion_mode == 35)
          mz_diff_adduct_search <- -34.9694
        else # ion_mode == -1
          mz_diff_adduct_search <- 1.007284 
      }
      
      repeat
      {
        id <- read_id(id_prompt)
        
        if (id == "0")
          break()
        else if (id == "q")
          q("no")
        
        # create the best results file if it does not exists and append the best results to it
        if (!("metfrag_pubchem_best_results.csv" %in% list.files(output_path)))
          # write the result table header
          write.table(list("msclusterID", "ionMode_PubChem", "fragmenterScore_PubChem", "numExplPeaks_PubChem",
                           "identifier_PubChem", "SMILES_PubChem", "molecularFormula_PubChem"), 
                      file = file.path(output_path, "metfrag_pubchem_best_results.csv"), 
                      row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
        if (id != "all")
        {
          # search the provided id
          run_metfrag(id)
        } else {
          # search all ids 
          # ask for count file to remove blank ids from the search
          count_file <- read_countfile(file_prompt)
          if (count_file == "0")
            break()
          else if (count_file== "q")
            q("no")
          else if (count_file == "all")
          {
            blank_ids <- c()
          } else {
            # read count file, skip first row with count file have correlation scores
            count_file <- read.csv(count_file, stringsAsFactors = FALSE,
                                   comment.char = "", strip.white = TRUE,
                                   skip = ifelse(grepl("BIOACTIVITY",
                                                       readLines(count_file, n = 1)),
                                                 1, 0))
            # check if number of spectra in the count file match with mgf
            if (length(msn_sample) != nrow(count_file))
            {
              warning("The number of spectra in the mgf do not match with the number of ", 
                      "spectra in the provided count file. Search aborted.", call. = FALSE)
              break()
            }
            if ("BLANKS_TOTAL" %in% names(count_file))
              blank_ids <- count_file[count_file$BLANKS_TOTAL > 0, "msclusterID"]
            else
              blank_ids <- c()
            rm(count_file)
          }
          
          # do not search ids that already have existing results
          existing_ids <- sapply(list.files(output_path)[endsWith(list.files(output_path), 
                                                                  paste0("_", ion_mode, "_PubChem.csv"))],
                                 function(x) as.numeric(strsplit(x, "_")[[1]][[1]]))
          
          n_search <- length(scan_indexes[!(scan_indexes %in% c(existing_ids, blank_ids))])
          message("\n** Running MetFrag Identification with PubChem Database**")
          ti <- Sys.time()
          
          # search the best correlated idxs
          search_best_results <- lapply(scan_indexes[!(scan_indexes %in% c(existing_ids, blank_ids))], 
                                        run_metfrag, TRUE)
          tf <- Sys.time()
          message("\n** Done seaching for ", n_search, 
                  " PubChem identifications in ", 
                  round(tf-ti, 2), " ", units(tf-ti), " **")
        }
      }
    }
  }
}