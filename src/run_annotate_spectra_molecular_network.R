##
# RUN annotation step
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
# default args
modification_file <- file.path(script_path(),"../rules/np3_modifications.csv")
mz_tol <- 0.025
fragment_tol <- 0.025
rt_tol <- 2
ms2_peaks_intensity_cutoff <- 15 # absolute intensity cutoff for fragmented MS2 peaks scaled from 0 to 1000
scale_factor <- 0.5 
ion_mode <- "+"  # + or -
table_limit_size <- 3000  # set max number of rows to process in a chunck

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  message("Two arguments must be supplied to annotate the clean counts of spectra and create the molecular network of annotations:\n",
       " 1 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks;\n",
       " 2 - Path to the output data folder, inside the outs directory of the clustering result folder. ",
       "It should contain the mgf folder, the peak area count CSV and the spectra count CSV. The data name will be extracted from here.;\n",
       " 3 - Path to the CSV table file with the modification rules for adducts, dimers/trimers and multi charge (default 'rules/np3_modifications.csv');\n",
       " 4 - The precursor m/z tolerance, for the numerical rules applied to the precursor m/z (default to 0.025 Da);\n",
       " 5 - The fragments m/z tolerance, for matching with the m/z of MS2 peaks in the fragments rule (default to 0.025 Da);\n",
       " 6 - The retention time tolerance to enlarge the peak limits (rtMin and rtMax) when searching for other concurrent peaks (default to 2s);\n",
       " 7 - The absolute intensity cutoff for keeping MS2 fragmented peaks with an intensity >= x and to use them in the in-source fragment variant annotation. The MS2 fragmented peaks intensity range from 0 to 1000. (default to 15);\n",
       " 8 - The ionization mode, one of -1 or 1 (default to 1);\n",
       " 9 - The scale factor used in the clean step, to be applied in the absolute MS2 intensity cutoff. (default 0.5);\n",
       " 10 - The maximum number of spectra (rows) to be processed at a time (default to 3000 rows).")
  stop(call.=FALSE)
} else {
  path_batch_metadata <- file.path(args[[1]])
  if (!file.exists(path_batch_metadata))
  {
    stop("The CSV batch metadata file '", path_batch_metadata, 
         "' do not exists. Provide a valid path to where the metadata is located.")
  }
  
  output_path <- file.path(args[[2]])
  if (!dir.exists(output_path))
  {
    stop("The job output folder '", output_path, 
         "' do not exists. Provide a valid path to where the job final result is located.")
  }
  
  if (length(args) >= 10) {
    modification_file <- file.path(args[[3]])
    if (!file.exists(modification_file))
    {
      stop("The rules file '", modification_file,
           "' with the accepted spectra modifications do not exists. ",
           "Provide a valid path to where the rules file is located.")
    }
    mz_tol <- as.numeric(args[[4]])
    fragment_tol <- as.numeric(args[[5]]) 
    rt_tol <- as.numeric(args[[6]])
    ms2_peaks_intensity_cutoff <- as.numeric(args[[7]])
    
    ion_mode <- as.numeric(args[[8]])
    if (!(ion_mode %in% c(-1,1)))
      stop("The ion_mode arg must be a numeric value indicating the precursor ", 
           "ion mode. One of the following valid numeric values corresponding ", 
           "to a ion adduct type: '1' = positive, or '-1' = negative",  call. = FALSE)
    ion_mode <- ifelse(ion_mode > 0, "+", "-")
    
    scale_factor <- as.numeric(args[[9]])
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
    table_limit_size <- as.numeric(args[[10]])
    if (table_limit_size < 100) {
      warning("The max number of spectra (rows) to be processed at a time was too small. ",
              "Setting it to 100.", call. = FALSE)
      table_limit_size <- 100
    }
  } else if (length(args) >= 3) {
    warning("Default values used for arguments 3 to 10. ",
            "All arguments must be supplied to overwrite the defaults.\n")
  }
}

source(file.path(script_path(), "annotate_spectra_molecular_network.R"))

annotate_spectra_table_network(output_path,  
                               path_batch_metadata, 
                               modification_file=modification_file, 
                               mz_tol=mz_tol, 
                               fragment_tol=fragment_tol, 
                               rt_tol=rt_tol,
                               ms2_peaks_intensity_cutoff=ms2_peaks_intensity_cutoff,  
                               scale_factor=scale_factor, 
                               ion_mode=ion_mode,  # + or -
                               table_limit_size=table_limit_size  
)
