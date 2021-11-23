## ----load-libs, message = FALSE--------------------------------------------
cat("Loading package XCMS...\n")
suppressPackageStartupMessages(library(xcms))
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

.Last <- function() {
  graphics.off() # close devices before printing
  close(filestdin)
}

filestdin <- file("stdin")

read_aggfun <- function(msg="Enter aggregation function: ")
{
  cat(msg)
  agg_fun <- gsub(" ", "", readLines(filestdin,1))
  
  if (agg_fun == "")
  {
    cat("max\n")
    return("max")
  }
  
  if (agg_fun %in% c("q", "0"))
    return(agg_fun)
  
  if (!(agg_fun %in% c("sum", "max", "mean", "min")))
  {
    message("\n* invalid aggregation function *\n")
    return(read_aggfun(msg))
  }
  
  return(agg_fun)
}

read_range <- function(msg="Enter range: ")
{
  cat(msg)
  n <- gsub(" ", "", readLines(filestdin,1))
  
  if (n == "")
  {
    cat("all\n")
    return("all")
  }
  
  if (n %in% c("q", "0", "all"))
    return(n)
  
  if (!grepl("^[0-9]+\\.?[0-9]*,[0-9]+\\.?[0-9]*$",n))
  {
    message("\n* invalid input *\n")
    return(read_range(msg))
  }
  
  return(type.convert(strsplit(n, split = ",")[[1]],
                      numerals = "warn.loss"))
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

read_view <- function(msg="Enter chr view: ")
{
  cat(msg)
  chr_view <- gsub(" ", "", readLines(filestdin,1))
  
  if (chr_view == "")
  {
    cat("all\n")
    return("all")
  }
  
  if (chr_view %in% c("q", "0"))
    return(chr_view)
  
  if (!(chr_view %in% c("batch", "single", "all")))
  {
    chr_view <- strsplit(chr_view, split = ",")[[1]]
    
    # check if they exist in the metadata file, error if not
    if (!all(chr_view %in% batch_metadata$SAMPLE_CODE))
    {
      message("\n* invalid samples codes ",
              paste(chr_view[ which(!(chr_view %in% batch_metadata$SAMPLE_CODE))], collapse = ", "),
              " *\n")
      
      return(read_view(msg))
    }
  }
  
  return(chr_view)
}

read_legendpos <- function(msg="Enter legend pos: ")
{
  cat(msg)
  pos <- gsub(" ", "", readLines(filestdin,1))
  
  if (pos == "")
  {
    cat("topright\n")
    return("topright")
  }
  
  if (pos %in% c("q", "0"))
    return(pos)
  
  if (!(pos %in% c('bottomright', 'bottom', 'bottomleft', 'left', 'topleft',
                   'top', 'topright', 'right', 'center')))
  {
    message("\n* invalid file type *\n")
    return(read_legendpos(msg))
  }
  
  return(pos)
}

agg_prompt <- "\n[*      ] \nEnter the function to be used to aggregate intensity values across the mz value range for the same retention time:
      q | 0 : quit
      ENTER : set default equals 'max'
        sum : total ion chromatogram (TIC)
        max : base peak chromatogram (BPC)
        min : minimum intensity
       mean : mean intensity
aggregationFun = "

file_prompt <- "\n[*******] \nEnter the image filename to save to a PNG (complete path):
      q : quit
      0 : back
filename = "

title_prompt <- "\n[*****  ] \nEnter the plot title:
      q : quit
      0 : back
  ENTER : set empty title equals ''
title = "

view_prompt <- "\n[****   ] \nEnter the chromatogram view determining what need to be plotted and how:
          q : quit
          0 : back
      ENTER : set default equals 'all'
      batch : plot all samples chromatograms grouped by data collection batch, as specified in the metadata file, one image per collection
     single : plot each sample chromatogram individually
        all : plot all samples chromatograms together
  s1,s2,... : plot together the chromatograms of the given samples codes, as specified in the metadata file, separated by comma
view = "

legend_prompt <- "\n[****** ] \nEnter the relative position of the legend to be placed inside of the chromatogram plot image. Allowed values are: 'bottomright', 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right' and 'center':
      q : quit
      0 : back
  ENTER : set default equals 'topright'
pos = "

save_chr_PNG <- function(chr, main, legend, legend_pos, filename, color_bio)
{
  mardef <- par("mar")
  png(file= paste0(filename, ".png"), width=1500, height=800, res = 100)
  par(mar = c(4.5, 4.5, 4, 2) + 0.1)
  plot(chr, col = color_bio, main = main, cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.9)
  grid()
  if (legend[[1]] != "")
  {
    legend(x=legend_pos, legend=legend, fill = color_bio, bty = "n", cex = 0.8)
  }
  
  par(mar = mardef)
  dev.off()
}

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Three arguments must be supplied to extract the chromatograms of raw data collections:\n",
       " 1 - The data name",
       " 2 - Path to the CSV batch metadata file containing filenames, sample codes and data collection batches;\n",
       " 3 - Path to the raw data folder;\n", call. = FALSE)
} else {
  data_name <- args[[1]]
  
  path_batch_metadata <- file.path(args[[2]])
  # validate input
  if (!file.exists(path_batch_metadata))
  {
    stop("The CSV batch metadata file '", path_batch_metadata,
         "' do not exists. Provide a valid path to where the metadata is located.")
  }
  path_batch_metadata <- normalizePath(path_batch_metadata)
  
  path_raw_data <- file.path(args[[3]])
  if (!dir.exists(path_raw_data))
  {
    stop("The raw data file folder '", path_raw_data,
         "' do not exists. Provide a valid path to where the raw data is located.")
  }
  path_raw_data <- normalizePath(path_raw_data)
  
  if (length(args) > 3) {
    list_samples_codes <- strsplit(args[[4]],',')[[1]]
  } else {
    list_samples_codes <- character(0)
  }
}

batch_metadata <- readMetadataTable(path_batch_metadata)
if (length(list_samples_codes) > 0) {
  # filter only the samples in the provided list
  batch_metadata <- batch_metadata[batch_metadata$SAMPLE_CODE %in% list_samples_codes,]
}
# read.csv(path_batch_metadata, stringsAsFactors = FALSE,
#                            comment.char = "", strip.white = TRUE)
# if (!all(c("FILENAME", "SAMPLE_CODE", "DATA_COLLECTION_BATCH", "SAMPLE_TYPE") %in% names(batch_metadata)))
# {
#   stop("Wrong batch metadata file format. It should have at least 4 columns named as follow:\n",
#        "- FILENAME: must contain the raw MS/MS file name\n",
#        "- SAMPLE_CODE: must contain a unique code identifying each file\n",
#        "- DATA_COLLECTION_BATCH: must contain a numeric index indicading which files should be grouped in a distinct batch",
#        "- SAMPLE_TYPE: must be 'blank' if the file is a blank, 'bed' or 'control' if the file is a control and 'sample' otherwise",
#        call. = F)
# }
# batch_metadata$SAMPLE_TYPE <- tolower(batch_metadata$SAMPLE_TYPE)
batch_metadata$FILENAME <- file.path(path_raw_data, batch_metadata$FILENAME)

# set the samples colors: all blank samples will be grey and not blank will be
# collered with its respective bioactivity score
# check if the bioactivity was an inhibition or an activation
if ("INHIBITION" %in% names(batch_metadata))
{
  bioactivity <- "INHIBITION"
} else {
  bioactivity <- "ACTIVATION"
}
batch_metadata$color <- "#00000060" 
if (!is.null(batch_metadata[[bioactivity]])) # has bioactivity score
{
  batch_metadata$color[batch_metadata$SAMPLE_TYPE != "blank"] <- rainbow(length(batch_metadata$SAMPLE_TYPE != "blank"))[
    order(order(batch_metadata[[bioactivity]][batch_metadata$SAMPLE_TYPE != "blank"]))]
} else # no bioactivity score, color the ramaining samples uniformaly
{
  batch_metadata$color[batch_metadata$SAMPLE_TYPE != "blank"] <- rainbow(sum(batch_metadata$SAMPLE_TYPE != "blank"))
}

message("Reading ", data_name, " raw data...\n")
# Create a phenodata data.frame
pd <- data.frame(sample_name = batch_metadata$SAMPLE_CODE,
                 sample_group = batch_metadata$DATA_COLLECTION_BATCH,
                 class = ifelse(batch_metadata$SAMPLE_TYPE == "blank", "blank", "sample"),
                 stringsAsFactors = FALSE)

raw_data <- readMSData(files = batch_metadata$FILENAME,
                       pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk")

# get mz range and rt range
rt_range <- range(raw_data@featureData@data$retentionTime)
rt_range <- paste(c(floor(rt_range[1]), ceiling(rt_range[2])), collapse = ",")
mz_range <- paste(c(round(min(raw_data@featureData@data$lowMZ), 4), 
                    round(max(raw_data@featureData@data$highMZ), 4)), 
                  collapse = ",")
mzrPrompt <- paste0("\n[**     ] \nEnter mz range [", mz_range,"] for the MS data slice:
      q : quit
      0 : back
  ENTER : set default equals 'all'
    x,y : numeric lower mz, upper mz values
mzr = ")

rtrPrompt <- paste0("\n[***    ] \nEnter rt range [", rt_range, "] for the retention time slice:
      q : quit
      0 : back
  ENTER : set default equals 'all'
    x,y : numeric lower rt, upper rt boundaries or 'all'
rtr = ")

# create preview device
X11(width = 8, height = 5, bg = "white", title = "chromatogram preview - view default")

repeat
{
  agg_fun <- read_aggfun(agg_prompt)
  
  if (agg_fun %in% c("q", "0"))
    q("no")
  
  repeat
  {
    mzr <- read_range(mzrPrompt)
    
    # TODO trycatch chromatogram
    
    if (length(mzr) == 1)
    {
      if (mzr[[1]] == "0")
        break()
      else if (mzr[[1]] == "q")
        q("no")
      else # mzr == all
      {
        message("\nExtracting chromatogram...")
        chr <- chromatogram(raw_data, aggregationFun = agg_fun)
        main <- paste0("Chromatogram ", agg_fun)
      }
    } else { # a range was provided
      message("\nExtracting chromatogram...")
      chr <- chromatogram(raw_data, aggregationFun = agg_fun, mz = mzr)
      main <- paste0("Chromatogram ", agg_fun, " mz range [", mzr[[1]],
                     ",", mzr[[2]], "]")
    }
    
    # plot chr preview
    plot(chr,
         col = batch_metadata$color,
         main = main)
    legend(x="topright", legend=paste(batch_metadata$SAMPLE_CODE, ifelse(is.null(batch_metadata[[bioactivity]]), "",
                                                                         paste("%", batch_metadata[[bioactivity]]))),
           fill = batch_metadata$color, bty = "n",
           cex = 0.5)
    grid()
    
    repeat
    {
      rtr <- read_range(rtrPrompt)
      
      if (length(rtr) == 1) {
        if (rtr == "0")
          break()
        else if (rtr == "q")
          q("no")
        else # rtr == "all"
          chr_rt <- chr
      } else { # a range was provided
        message("\nFiltering chromatogram...")
        chr_rt <- Chromatograms(lapply(chr, filterRt, rt = rtr), 
                                nrow = 1, phenoData = pd)
      }
      
      plot(chr_rt,
           col = batch_metadata$color,
           main = main)
      legend(x="topright", 
             legend=paste(batch_metadata$SAMPLE_CODE, ifelse(is.null(batch_metadata[[bioactivity]]), "",
                                                             paste("%", batch_metadata[[bioactivity]]))),
             fill = batch_metadata$color, bty = "n",
             cex = 0.5)
      grid()
      
      repeat
      {
        newplot <- FALSE
        chr_view <- read_view(view_prompt)
        
        if (chr_view[[1]] == "0")
          break()
        else if (chr_view[[1]] == "q")
          q("no")
        
        # TODO update preview?
        
        repeat
        {
          plot_main <- read_title(title_prompt)
          
          if (plot_main == "0")
            break()
          else if (plot_main == "q")
            q("no")
          
          repeat
          {
            legend_pos <- read_legendpos(legend_prompt)
            
            if (legend_pos == "0")
              break()
            else if (legend_pos == "q")
              q("no")
            
            repeat
            {
              filename <- read_filename(file_prompt)
              
              if (filename == "0")
                break()
              else if (filename == "q")
                q("no")
              
              message("\n\n---------------\n SAVED IMAGES: \n")
              if (chr_view == "all") # plot all selected chr together
              {
                save_chr_PNG(chr_rt[1,], plot_main,
                             paste(batch_metadata$SAMPLE_CODE, ifelse(is.null(batch_metadata[[bioactivity]]), "",
                                                                      paste("%", batch_metadata[[bioactivity]]))),
                             legend_pos, filename,
                             batch_metadata$color)
                message(" - ", filename, ".png")
              } else if (chr_view == "batch") {
                n_batches <- unique(batch_metadata$DATA_COLLECTION_BATCH)
                
                for (i in n_batches)
                {
                  chr_samples <- batch_metadata[batch_metadata$DATA_COLLECTION_BATCH == i, "SAMPLE_CODE"]
                  chr_i <- which(batch_metadata$SAMPLE_CODE %in% chr_samples)
                  save_chr_PNG(chr_rt[1, chr_i], plot_main, paste(chr_samples, 
                                                                  ifelse(is.null(batch_metadata[[bioactivity]]), "",
                                                                         paste("%", batch_metadata[[bioactivity]][chr_i]))),
                               legend_pos, paste0(filename, "_batch_", i),
                               batch_metadata$color[chr_i])
                  message(" - ", paste0(filename, "_batch_", i, ".png"))
                }
              } else if (chr_view == "single") {
                for (i in seq_len(nrow(batch_metadata)))
                {
                  save_chr_PNG(chr_rt[1, i], plot_main, "", legend_pos,
                               paste0(filename, "_sample_", batch_metadata$SAMPLE_CODE[i]),
                               "#00000060")
                  message(" - ", paste0(filename, "_sample_", batch_metadata$SAMPLE_CODE[i], ".png"))
                }
              } else { # selected samples
                chr_i <- which(batch_metadata$SAMPLE_CODE %in% chr_view)
                save_chr_PNG(chr_rt[1,chr_i], plot_main,
                             paste(batch_metadata$SAMPLE_CODE[chr_i], 
                                   ifelse(is.null(batch_metadata[[bioactivity]]), "",
                                          paste("%", batch_metadata[[bioactivity]][chr_i]))),
                             legend_pos, filename,
                             batch_metadata$color[chr_i])
                message(" - ", filename, ".png")
              }
              message("\n---------------\n\n**** NEW PLOT ****")
              
              newplot <- TRUE
              break()
            }
            if (newplot) break() # return to the view selection
          }
          if (newplot) break() # return to the view selection
        }
      }
    }
  }
}