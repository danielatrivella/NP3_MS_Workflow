## ----load-libs, message = FALSE--------------------------------------------
cat("Loading package MSnBase...\n")
suppressPackageStartupMessages(library(MSnbase))

filenames_saved <- c() # all the files saved in one interactive section
.Last <- function() {
  graphics.off() # close devices before printing
  close(filestdin)
  
  if (length(filenames_saved) > 0)
  {
    message("\n\n---------------\n SAVED IMAGE: \n", paste(unique(filenames_saved), collapse = "\n"))
    message("---------------")
  }
}

filestdin <- file("stdin")

read_ids <- function(preids)
{
  if (preids[[1]] == "t") { # trimmed spectra
    msg <- (ids_t_pre_view_prompt)
  } else if (preids[[1]] >= 0 || preids[[1]] == "") { # not trimmed
    msg <- (ids_pre_view_prompt)
  } else { # no spectra
    msg <- (ids_prompt)
  }
  
  cat(msg)
  n <- gsub(" ", "", readLines(filestdin,1))
  
  if (preids[[1]] == "t")  { # if its trimmed
    if (n %in% c("q", "0", ""))
      return(n)
  } else if (preids[[1]] >= 0 || preids[[1]] == "") { # if its not the first plot and its not trimmed 
    if (n %in% c("q", "0", "t", ""))
      return(n)
  } else #  first start up plot 
  {
    if (n %in% c("q", "0"))
      return(n)
  } 
  
  # check if entry was x,y
  if (!grepl("^[0-9]+,[0-9]+$",n))
  {
    message("\n* invalid input format *\n")
    return(read_ids(preids))
  } 
  
  n_ids <- as.numeric(strsplit(n, split = ",")[[1]])
  if (!all(n_ids %in% scan_indexes)) {
    message("\n* invalid IDs *\n")
    return(read_ids(preids))
  }
  
  return(n_ids)
}

read_tol <- function(msg="Enter frag tol: ")
{
  cat(msg)
  n <- gsub(" ", "", readLines(filestdin,1))
  
  if (n == "")
  {
    cat(mz_tol, "\n")
    return(mz_tol)
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

read_filename <- function(msg="Enter filename: ")
{
  cat(msg)
  filename <- readLines(filestdin,1)
  
  if (filename %in% c("q", "0"))
    return(filename)
  
  if (filename == "" || !dir.exists(dirname(filename)))
  {
    message("\n* invalid path *\n")
    return(read_filename(msg))
  }
  
  return(filename)
}

read_title <- function(msg="Enter plot title: ")
{
  cat(msg)
  
  return(readLines(filestdin,1))
}

read_fun <- function(msg="Enter sim function: ")
{
  cat(msg)
  sim_fun <- gsub(" ", "", readLines(filestdin,1))
  
  if (sim_fun == "")
  {
    cat("dotproduct\n")
    return("dotproduct")
  }
  
  if (sim_fun %in% c("q", "0"))
    return(sim_fun)
  
  if (!(sim_fun %in% c("dotproduct", "common", "cor")))
  {
    message("\n* invalid comparison function *\n")
    return(read_fun(msg))
  }
  
  return(sim_fun)
}

read_preaction <- function(msg="Enter preview action: ")
{
  cat(msg)
  action <- gsub(" ", "", readLines(filestdin,1))
  
  if (action == "")
  {
    cat("save\n")
    return("s")
  }
  
  if (!(action %in% c("0", "q", "t")))
  {
    message("\n* invalid action *\n")
    return(read_preaction(msg))
  }
  
  return(action)
}

# ids prompt before first plot
ids_prompt <- "\n|====  | \nEnter the two spectra SCANS's numbers (msclusterID's) to be compared and plotted against each other:
      q : quit
      0 : back
    x,y : x and y are SCANS's numbers (msclusterID's) present in the MGF file
spectra ID's = "

# preview new plot 
# ENTER so aparece se existir um plot
ids_pre_view_prompt <- "\n|====  | \nEnter the two spectra SCANS number (msclusterID's) to be compared and plotted against each other:
      q : quit
      0 : back
      t : trim bigger spectrum based on the limits of the smaller spectrum fragmentation peaks 
  ENTER : save plot
    x,y : x and y are SCANS's numbers (msclusterID's) present in the MGF file
spectra ID's = "

# pre view trimmed plot
ids_t_pre_view_prompt <- "\n|====  | \nEnter the two spectra SCANS number (msclusterID's) to be compared and plotted against each other:
      q : quit
      0 : back
  ENTER : save plot
    x,y : x and y are SCANS's numbers (msclusterID's) present in the MGF file
spectra ID's = "

tol_prompt <- "\n|=     | \nEnter the fragments mz tolerance in Daltons to find and colour common peaks between two spectra:
    q : quit
ENTER : set default equals '0.05'
  <X> : a numeric with the fragment peaks tolerance
tolerance = "

file_prompt <- "\n|======| \nEnter the image filename to save to a PNG (complete path):
q : quit
0 : back
filename = "

title_prompt <- "\n|===== |  \nEnter the plot title:
    q : quit
    0 : back
ENTER : set empty title equals ''
title = "

fun_prompt <- "\n|===   |  \nEnter the function to be used in the spectra similarity computation:
         q : quit
         0 : back
     ENTER : set default equals 'dotproduct'
    common : the comparison is based on the number of common peaks, or
       cor : the Pearson correlation, or
dotproduct : the dot product
fun = "

compare_spectra_plot <- function(i, j, main = "")
{
  plot(mgf_data[[i]], mgf_data[[j]], tolerance = mz_tol, 
       relative = relative, main = main)
  grid()
  
  return(round(compareSpectra(mgf_data[[i]], mgf_data[[j]], 
                        sim_fun,
                        binSize = mz_tol),2))
}

compare_spectra_plot_trimmed <- function(i, j, main = "")
{
  s1 <- mgf_data[[i]]
  s2 <- mgf_data[[j]]
  
  if (s1@precursorMz < s2@precursorMz)
    s2 <- trimMz(s2, c(min(s1@mz) - mz_tol, s1@precursorMz + mz_tol*2))
  else
    s1 <- trimMz(s1, c(min(s2@mz) - mz_tol, s2@precursorMz + mz_tol*2))
  
  plot(s1, s2, tolerance = mz_tol, 
       relative = relative, main = main)
  grid()
  
  return(round(compareSpectra(s1, s2, 
                        sim_fun,
                        binSize = mz_tol),2))
}

save_spectra_plot <- function(i, j, main, filename, trim_spec)
{
  png(file= paste0(filename, ".png") ,width=1600, height=1000)
  if (trim_spec)
    compare_spectra_plot_trimmed(i, j, main)
  else
    compare_spectra_plot(i, j, main)
  dev.off()
}

# defaults
mz_tol <- 0.05
relative <- FALSE # for absolute tolerance in Daltons

path_mgf <- "../Bra-mzXML/Bra_RT1_2_area/outs/Bra_RT1_2_area/mgf/Bra_RT1_2_area_all_clean.mgf"

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Three arguments must be supplied to extract the chromatograms of raw data collections:\n",
       " 1 - The data name",
       " 2 - Path to the MGF containing the clustering result.", call. = FALSE)
} else {
  data_name <- args[[1]]
  
  path_mgf <- file.path(args[[2]])
  if (!file.exists(path_mgf))
  {
    stop("The MGF file '", path_mgf,
         "' do not exists. Provide a valid path to where the MGF file is located.")
  }
  path_mgf <- normalizePath(path_mgf)
}

message("Reading ", data_name, " MGF...\n")

## Create a phenodata data.frame
pd <- data.frame(sample_name = data_name,
                 sample_group = 1,
                 class = "sample",
                 stringsAsFactors = FALSE) 

# read mgf and get scan index
mgf_data <- readMgfData(path_mgf, new("NAnnotatedDataFrame", pd), 
                         verbose = T)
scan_indexes <- scanIndex(mgf_data)

# create preview device
X11(width = 6, height = 5, bg = "white", title = "spectra comparision preview")
startup <- TRUE

repeat
{
  mz_tol <- read_tol(tol_prompt)
  
  if (mz_tol == "q")
    q("no")
    
  repeat
  {
    sim_fun <- read_fun(fun_prompt)
    
    if (sim_fun == "0")
      break()
    else if (sim_fun == "q")
      q("no")
    
    ids <- -1 # reset plots
    repeat 
    {
      newplot <- FALSE
      ids <- read_ids(ids)
      
      # new plot
      if (length(ids) > 1)
      {
        # trim spectra and update preview
        i <- match(ids[[1]], scan_indexes)
        j <- match(ids[[2]], scan_indexes)
        title_plot <- paste("SCANS:", ids[[1]], "x", ids[[2]])
        
        message("\n------------------------", 
                "\nspectra similarity = ", 
                compare_spectra_plot(i,j, main = title_plot), 
                " (", sim_fun, ")", "\n------------------------")
        trim_spec <- FALSE
        
        next() # loop preview
      } else if (ids == "t") {
        message("\n------------------------", 
                "\nspectra similarity trimmed = ", 
                compare_spectra_plot_trimmed(i,j, main = paste("TRIMMED", title_plot)),
                " (", sim_fun, ")", "\n------------------------")
        trim_spec <- TRUE
        
        next()
      } else if (ids == "0") {
        break()
      } else if (ids == "q") {
        q("no")
      }
      
      repeat
      {
        plot_main <- read_title(title_prompt)
        
        if (plot_main == "0")
          break()
        else if (plot_main == "q")
          q("no")
        
        repeat
        {
          filename <- read_filename(file_prompt)
          
          if (filename == "0")
            break()
          else if (filename == "q")
            q("no")
          
          save_spectra_plot(i, j, plot_main, filename, trim_spec)
          
          filename <- paste0(filename, ".png")
          message("\n\n---------------\n SAVED IMAGE: \n", filename)
          message("---------------\n\n**** NEW PLOT ****")
          filenames_saved <- c(filenames_saved, filename)
          
          newplot <- TRUE
          break()
        }
        if (newplot) break() # return to the ids selection
      }
    }
  }
}