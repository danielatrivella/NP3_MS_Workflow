suppressPackageStartupMessages(library(dplyr))
Rcpp::sourceCpp('src/read_mgf_peak_list_R.cpp')
# compute_peak_area(processed_data_path, 
#                   ms_spectra_count$msclusterID,
#                   lapply(ms_spectra_count$scans, function(x) strsplit(x, ";")[[1]]),
#                   lapply(ms_spectra_count$peakIds, function(x) strsplit(x, ";")[[1]]),
#                   batch_metadata)

# retrieve the peak areas and max base peak int for each msclusterID based on
# the scans_counts and check if the peakIds are correct
compute_peak_area <- function(processed_data_path, msclusterIDs, scans_count, 
                              peakIds_count, metadata)
{
  peak_areas <- Reduce(bind_cols, lapply(metadata$SAMPLE_CODE, function(x)
  {
    mgf_data <- readMgfHeader(file.path(processed_data_path, paste0(x, '_peak_info.mgf')))
    mgf_data$scans <- paste0(mgf_data$scans, "_", x)
    
    peak_area_x <- bind_rows(lapply(seq_len(length(scans_count)), function(i)
    {
      scans <- scans_count[[i]]
      scans_x <- grepl(pattern = paste0("$", x, "$"), fixed = TRUE, 
                       x = paste0(sub(pattern = "_", replacement = "$", scans), "$")) # add _ at the end to prevent matching with a sample with the same prefix
      # at least one scan appears in this sample, check consistency 
      if (any(scans_x)) {
        scans_header <- mgf_data[match(scans[scans_x], mgf_data$scans), ]
        peakIds <- peakIds_count[[i]][grepl(pattern = paste0("$", x, "$"), fixed = T, 
                                              x = paste0(sub(pattern = "_", 
                                                             replacement = "$", 
                                                             peakIds_count[[i]]), "$"))]
        
        # test if all peakIds are correct for this sample
        if (!all(peakIds %in% scans_header$peak_id) & 
            !all(scans_header$peak_id  %in% peakIds)) {
          stop("Wrong peakIds for msclusterId ", msclusterIDs[[i]], 
                  " in sample ", x, 
                  ". The following peaks are not present in the original MGF: ",
                  paste(peakIds[!(peakIds %in% scans_header$peak_id)],collapse = ","))
        }
        
        # sum the unique peak areas (fake peaks got integrated here because they have different peak_areas)
        real_peak_area_x <- sum(unique(trunc(scans_header$peak_area)))
        # return the maximum base peak found
        base_peak_int = max(scans_header$base_peak_int)
        
        list(x = real_peak_area_x, base_peak_int = base_peak_int)
      } else {
        list(x = 0, base_peak_int = 0)
      }
    }))
    names(peak_area_x) <- c(x, paste0('base_peak_int_',x))
    
    list(peak_area_x)
  }))
  if (length(metadata$SAMPLE_CODE) == 1) {
    peak_areas <- peak_areas[[1]]
  }
  base_peak_int <- apply(peak_areas[,startsWith(names(peak_areas), "base_peak_int")], 1, max)
  # test if all base peak ints were retrieved correctly
  if (any(base_peak_int == 0)) {
    stop("Wrong base peak intensities retrieved for msclusterId's ", 
         paste(msclusterIDs[(base_peak_int == 0)],collapse = ","))
  }
  peak_areas <- peak_areas[,!startsWith(names(peak_areas), "base_peak_int")]
  names(peak_areas) <- paste0(names(peak_areas), "_area")
  peak_areas$basePeakInt <- base_peak_int
  return(peak_areas)
}