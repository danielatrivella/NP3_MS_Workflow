# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Three arguments must be supplied to join the GNPS library identification results to the NP3 count files (clustering and clean counts):\n", 
       " 1 - Path to the file of the GNPS library identification output located inside the folder named 'clusterinfo';\n", 
       " 2 - Path to the file of the GNPS library identification output located inside the folder named 'result_specnets_DB';\n", 
       " 3 - Path to any of the count tables (peak_area or spectra) resulting ",
       "from the NP3 clustering or clean steps. If the peak_area is informed and ",
       "the spectra file exists in the same path (or the opposite), it will merge the GNPS results to both files.\n",
       call.=FALSE)
} else {
  cluster_info_path <- file.path(args[[1]])
  if (cluster_info_path != "" && !file.exists(cluster_info_path))
  {
    stop("The clusterinfo file '", cluster_info_path, 
         "' from the GNPS result does not exists. Provide a valid path to where it is located.")
  }
  
  lib_idres_path <- file.path(args[[2]])
  if (!file.exists(lib_idres_path))
  {
    stop("The result specnets DB file '", lib_idres_path, 
         "' from the GNPS result does not exists. Provide a valid path to where it is located.")
  }
  
  ms_count_path <- file.path(args[[3]])
  if (!file.exists(ms_count_path))
  {
    stop("The count file '", ms_count_path, 
         "' from the NP3 results does not exists. Provide a valid path to where it is located.")
  }
}

join_gnps_result_count_tables <- function(ms_count_path, gnps_columns_keep, lib_idres_scans) 
{
  # bind the gnps result columns to the count tables
  ms_count <- read.csv(ms_count_path, stringsAsFactors = FALSE,
                       comment.char = "", 
                       strip.white = TRUE)
  ms_count[gnps_columns_keep[-1]] <- NULL # remove previous columns
  ms_count <- merge(ms_count, lib_idres_scans, all.x = TRUE, by = "msclusterID")
  gnps_columns_keep = gnps_columns_keep[-1] # remove msclusterID
  ms_count <- ms_count[,c(names(ms_count)[!(names(ms_count) %in% gnps_columns_keep)], 
                         gnps_columns_keep)]
  if (sum(!is.na(unique(ms_count$gnps_SpectrumID))) != length(unique(lib_idres_scans$gnps_SpectrumID)))
  {
    stop("ERROR. An inconsistency happen when joining the GNPS results to the ",
         "NP3 count table, probably missing IDs.")
  }
  # write the output
  write.csv(ms_count, ms_count_path,
            row.names = FALSE, quote = TRUE)
}

suppressPackageStartupMessages(library(dplyr))

cat("\n*** Joining the GNPS library indentification results to the NP3 count files ***\n\n")
t0 <- Sys.time()

# read gnps results spec net file with the libraries identifications results
lib_idres <- read.table(lib_idres_path,
                        header = T, sep = "\t", stringsAsFactors = F, 
                        comment.char = "", quote="\"")
if (cluster_info_path != "")
{
  # this is a result from molecular networking
  # read gnps clusterinfo file
  cluster_info <- read.table(cluster_info_path,
                             header = T, sep = "\t", stringsAsFactors = F, 
                             comment.char = "")
  names(cluster_info) <- sub("X.","", names(cluster_info))
  # fix the gnps result when dealing with zero scans value, it gives a value of 1 to the scans equals 0
  if (cluster_info$Scan[cluster_info$SpecIdx == 0] == 1 &&
      cluster_info$Scan[cluster_info$SpecIdx == 1] == 1)
  {
    # SpecIdx == 0 will have scan == 0 in the np3 mgf files, 
    # if the msclusterID = 0 was not removed by the noise filtering
    cluster_info$Scan[cluster_info$SpecIdx == 0] <- 0 
  }
  # the Scan column contains the clusterIdx
  # join the clusterinfo and the results spec tables using as key the ClusterIdx 
  # column in the first and the #Scan column in the second
  # rename the #Scan colum in the results spec table to ClusterIdx and join
  names(lib_idres) <- sub("X.Scan.", "ClusterIdx", names(lib_idres))
  lib_idres_scans  <- merge(lib_idres, cluster_info, 
                            all.x = TRUE, by = "ClusterIdx")
  rm(cluster_info)
} else {
  # this is a result from Library Search 
  # the Scan column contains the spectra scans
  names(lib_idres) <- sub("X.Scan.", "Scan", names(lib_idres))
  lib_idres_scans  <- lib_idres
}
rm(lib_idres)
# use the Scan column of the clusterinfo table to match with the 
# msclusterIDs of the count table
lib_idres_scans$msclusterID <- lib_idres_scans$Scan
# if it is the clean count tables match with the msclusterIDs present in the 
# column joinedIDs of the clean count tables
if (grepl("clean", ms_count_path))
{
  ms_count <- read.csv(ms_count_path, stringsAsFactors = FALSE,
                       comment.char = "",
                       strip.white = TRUE)[, c("msclusterID", "joinedIDs")]
  num_msclusterid <- nrow(ms_count)
  # get the joined ids
  ms_count <- ms_count[!is.na(ms_count$joinedIDs),]
  
  if (nrow(ms_count) > 0) {
    invisible(apply(ms_count, 1, function(x) 
    {
      joinedIDs <- as.numeric(strsplit(x[[2]], ";")[[1]])
      
      lib_idres_scans$msclusterID[
        lib_idres_scans$Scan %in% joinedIDs] <<- as.numeric(x[[1]])
    }))
  }
} else {
  ms_count <- read.csv(ms_count_path, stringsAsFactors = FALSE,
                       comment.char = "",
                       strip.white = TRUE)[, c("msclusterID")]
  num_msclusterid <- nrow(ms_count)
}
rm(ms_count)

# order by the msclusterID ascending and descending by the MQScore, to put the 
# best results by spectra on top
lib_idres_scans <- dplyr::arrange(lib_idres_scans, msclusterID, desc(MQScore))

# filter only the wanted columns
# merge the ion source with the instrument
gnps_columns_keep <- c('msclusterID', 'SpectrumID','Adduct','Smiles', 'InChIKey',
                       'CAS_Number','Compound_Name','LibMZ','MZErrorPPM',
                       'MQScore','LibraryQualityString','SharedPeaks','Organism',
                       'npclassifier_superclass', 'npclassifier_class', 'npclassifier_pathway', # old gnps columns for npclassifier result 'superclass','class','subclass',
                       'Ion_Source_Instrument')

#verify if there is a missing column and exclude it from the filter
gnps_columns_keep <- gnps_columns_keep[gnps_columns_keep %in% names(lib_idres_scans)]


lib_idres_scans$Ion_Source_Instrument <- paste(lib_idres_scans$Ion_Source, 
                                               lib_idres_scans$Instrument, 
                                               sep =  "_")
lib_idres_scans <- lib_idres_scans[gnps_columns_keep]
# remove duplicates, keep first occurence with higher MQScore
lib_idres_scans <- lib_idres_scans[!duplicated(lib_idres_scans[,c('msclusterID', 'SpectrumID')]),]

# remove first column with the msclusterID and
# rename columns adding a prefix equals gnps_
gnps_columns_keep[-1] <- paste0("gnps_",gnps_columns_keep[-1]) 
names(lib_idres_scans) <- gnps_columns_keep 

# aggregate the rows based on the same msclusterID using a ; as sep
n <- nrow(lib_idres_scans)
n_unique <- length(unique(lib_idres_scans$gnps_Compound_Name))
lib_idres_scans <- aggregate(lib_idres_scans, 
                             by = list(lib_idres_scans$msclusterID), 
                             "paste0", collapse=";")
# extract the final msclusterID, remove the repeated ones
lib_idres_scans$msclusterID <- sapply(lib_idres_scans$msclusterID, 
                                      function(x) strsplit(x, ";")[[1]][[1]])
lib_idres_scans <- lib_idres_scans[,gnps_columns_keep] # remove the group column
# in the smiles columns substitute the ; separator by a comma , for better visualization in cytoscape
lib_idres_scans$gnps_Smiles <- gsub(";", ",", lib_idres_scans$gnps_Smiles)

# bind the gnps result columns to the count tables
join_gnps_result_count_tables(ms_count_path, gnps_columns_keep, lib_idres_scans)
# write to the other count file if present
if (grepl("peak_area", ms_count_path) && file.exists(sub("peak_area", "spectra", 
                                                         ms_count_path))) 
{
  ms_count_path <- sub("peak_area", "spectra", ms_count_path)
  join_gnps_result_count_tables(ms_count_path, gnps_columns_keep, lib_idres_scans)
} else if (grepl("spectra", ms_count_path) && file.exists(sub("spectra","peak_area", 
                                                              ms_count_path)))
{
  ms_count_path <- sub("spectra", "peak_area",  ms_count_path)
  join_gnps_result_count_tables(ms_count_path, gnps_columns_keep, lib_idres_scans)
}

tf <- Sys.time()
cat("  * Done joining a total of", n, "GNPS identifications from",
        n_unique, "unique identifications to", nrow(lib_idres_scans),
        "consensus spectra in", 
        round(tf-t0, 2), units(tf-t0), "*\n\n")
cat("  * Spectra identification rate =", 
        round(nrow(lib_idres_scans)/num_msclusterid*100, 2), "%\n\n")
