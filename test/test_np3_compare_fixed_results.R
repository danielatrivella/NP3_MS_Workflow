##
# functions to execute equality tests for the NP3 results from the run, join_jobs and pre process commands
# - compare new results with fixed results - check the consistency in the steps
# - print any mismatching result detected and their differences by groups:
#   - tested groups in the Run and join_jobs command: Count_tables,MGFs,Similarity_matches_tables,Molecular_network_tables,Tremolo_identifications,GNPS_identifications,Quantification_report,Chemical_reports,Molecular_networking_reports"
#   - tested groups in the Pre process command: pre processed MGFs, MS1 tables
##

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

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
Rcpp::sourceCpp(file.path(script_path(), '../src/read_mgf_peak_list_R.cpp'))
source(file.path(script_path(), "../src/read_metadata_table.R"))

## args
output_path_test <- "/tmp/np3_ms_workflow/"
fixed_result_path <- "./test/L754_bacs/fixed_results/L754_bacs_all"
new_result_path <- "./test/L754_bacs/L754_bacs_all"
sim_w_cutoff="06"
topk="15"
maxComponentSize=200
minMachedPeaks=6
gnps_library_search_tool="gnps"
tremolo_exec=TRUE
##

## function to compare two tables 
## Where the fixed_table_path is the path to a fixed result of a count table;
# new_table_path is the path to a new result of the same count table of a similar job;
# table_type is the type source of the tables, could be clustering, clean, merge, similarity, quantification_report results;
# the job_name is the output_name of the provided fixed table; 
# the output_path_test_diff is the path to where the difference table will be stored in case the fixed table does
# not match with the new table; 
# the anti_join_by is a list of column names to be used in the anti_join to 
#   extract the different cells from the fixed and new table, 
#   this must be the set of columns in common between the two tables. 
#   This set of columns will be used to filter the two tables (columns expected to be equal). 
#   If this value is NULL all columns with a common name are used instead (default);
# sep is the separator character of the provided tables, used to parse and read them;
## Returns TRUE if the fixed table matches with the new table (all equal); otherwise
# returns a list with the differences in terms of columns, rows, cells and the
# output path to the difference table (cells present in the fixed table and absent in the new table, if any)
compare_tables <- function(fixed_table_path, new_table_path, 
                           table_type, job_name, output_path_test_diff, 
                           anti_join_by = NULL, sep=",") {
  if (sep == ",") {
    # interpret as csv
    fixed_table <- suppressMessages(read_csv(fixed_table_path, col_types = cols(.default = "c"), guess_max = 5000))
    new_table <- suppressMessages(read_csv(new_table_path, , col_types = cols(.default = "c"), guess_max = 5000))
  } else if (sep == "\t") {
    # interpret as tsv
    fixed_table <- suppressMessages(read_tsv(fixed_table_path, guess_max = 5000))
    new_table <- suppressMessages(read_tsv(new_table_path, guess_max = 5000))
  } else {
    return("Error in the separator type, not supported.")
  }
  if (!is.null(anti_join_by)) {
    # filter only the columns expected to be equal
    if (!all(anti_join_by %in% names(fixed_table))) {
      return(paste0("Error in the columns selected for the anti_join_by, the following are not present in the fixed_table: ",
             paste0(anti_join_by[!(anti_join_by %in% names(fixed_table))], collapse=","),"."))
    }
    fixed_table <- fixed_table[,anti_join_by]
    new_table <- new_table[,anti_join_by]
  }
  
  # runs equality test
  test_tables_equal <- all.equal(fixed_table, new_table)
  if (length(test_tables_equal) > 1 || test_tables_equal != TRUE)
  {
    test_tables_equal <- c(paste("  Differences in the fixed table (x) and the new table (y) of", table_type),
                           test_tables_equal)
    # creates df with differences present in fixed_table and missing/diff in new_table
    # suppress msg when anti_join_by is NULL - uses all columns for complete equality
    table_diff <- suppressMessages(anti_join(fixed_table, new_table, by=anti_join_by))
    if (is.null(anti_join_by))
    {
      anti_join_by <- "all common columns"
    }
    # computes number of different rows
    test_tables_equal <- c(test_tables_equal, 
                           paste0("  Number of rows in x but not in y (by=",paste(anti_join_by,collapse=","),"): ", 
                                  nrow(table_diff)))
    # store difference table if there is different rows/cells
    if (nrow(table_diff) > 0) {
      table_diff_out_path <- file.path(output_path_test_diff, 
                                       paste0(job_name, "_", table_type,"_diff.csv"))
      test_tables_equal <- c(test_tables_equal, 
                             paste0("  Difference table (cells present in the fixed table and absent in the new table) stored at: ", 
                                    table_diff_out_path))
      write_csv(table_diff, path=table_diff_out_path)
    }
  }
  
  if (length(test_tables_equal) > 1) {
    # concatenate the error messages
    test_tables_equal <- paste0(test_tables_equal, collapse = "\n")
  }
  
  return(test_tables_equal)
}

## Function to compare two MGFs files
# start by comparing the MGFs headers' table, for the matching ions, also 
# compare the MGFs peak lists - read the full peak list without filters and 
# check for equality in terms of the m/zs and intensities
## Where the fixed_mgf_path is the path to a MGF used as reference;
# new_mgf_path is the path to a MGF to be compared with;
# mgf_type is the step from which this MGF is from, one of preprocess, clustering or clean;
# output_path_test_diff is the path to where the difference mgf header table will be stored in case the fixed mgf header does
# not match with the new mgf header table; 
## Returns TRUE if both the mgf headers and peak lists match all equal; otherwise
# return an error message describing the difference found in the MGFs headers concatenated
# with any difference found the in the peak lists of the ions with matching headers
compare_mgfs <- function(fixed_mgf_path, new_mgf_path, 
                         mgf_type, job_name, output_path_test_diff) {
  # start reading the MGFs header tables and comparing them
  fixed_mgf_path <- normalizePath(fixed_mgf_path)
  new_mgf_path <- normalizePath(new_mgf_path)
  fixed_mgf_header <- readMgfHeader(fixed_mgf_path)
  new_mgf_header <- readMgfHeader(new_mgf_path)
  
  # compare the mgf header table
  # runs equality test
  test_mgf_header_tables_equal <- all.equal(fixed_mgf_header, new_mgf_header)
  if (length(test_mgf_header_tables_equal) > 1 || test_mgf_header_tables_equal != TRUE)
  {
    test_mgf_header_tables_equal <- c(paste("Differences in the MGF header fixed table (x) and the new table (y) of", 
                                            mgf_type),
                                      test_mgf_header_tables_equal)
    # creates df with differences present in fixed_table and missing/diff in new_table, using all columns
    table_diff <- suppressMessages(anti_join(fixed_mgf_header, new_mgf_header))
    
    # store difference table if there is different rows/cells
    if (nrow(table_diff) > 0) {
      # computes number of different rows
      test_mgf_header_tables_equal <- c(test_mgf_header_tables_equal, 
                                        paste0("Rows in x but not in y (by=all): ", 
                                               nrow(table_diff)))
      # store the diff table to the output test path
      table_diff_out_path <- file.path(output_path_test_diff, 
                                       paste0(job_name, "_mgf_header_", mgf_type,"_diff.csv"))
      test_mgf_header_tables_equal <- c(test_mgf_header_tables_equal, 
                             paste0("Difference MGF header table (cells present in the fixed table and absent in the new table) stored at: ", 
                                    table_diff_out_path))
      write_csv(table_diff, path=table_diff_out_path)
    }
    if (length(test_mgf_header_tables_equal) > 1) {
      # concatenate the error messages
      test_mgf_header_tables_equal <- paste0(test_mgf_header_tables_equal, collapse = "\n")
    }
    # extract the row names index of the scans with a common header
    idx_common_headers <- row.names(fixed_mgf_header[!(fixed_mgf_header$scans %in% table_diff$scans),])
  } else {
    idx_common_headers <- row.names(fixed_mgf_header)
  }
  
  # if there is any header in common, read the mgfs peaks list, 
  # filter the common header and check if their peaks list are equal
  if (length(idx_common_headers) > 0) {
    fixed_mgf_peaks_data <- as_tibble(readMgfPeaksList(fixed_mgf_path, 
                                                       bin_size=0.05, trim_mz=F, 
                                                       scale_factor=0.05, 
                                                       join_isotopic_peaks=F))
    new_mgf_peaks_data <- as_tibble(readMgfPeaksList(new_mgf_path, bin_size=0.05, 
                                                     trim_mz=F, scale_factor=0.05, 
                                                     join_isotopic_peaks=F))
    
    # check if the peaks list MZS and INTS are equal across the common ions
    test_mgf_peaks_equal <- sapply(idx_common_headers, function(i) {
      fixed_peaks_mzs <- fixed_mgf_peaks_data[[i,"MZS"]]
      new_peaks_mzs <- new_mgf_peaks_data[[i,"MZS"]]
      fixed_peaks_ints <- fixed_mgf_peaks_data[[i,"INTS"]]
      new_peaks_ints <- new_mgf_peaks_data[[i,"INTS"]]
      test_mzs <- all.equal(fixed_peaks_mzs,new_peaks_mzs)
      test_ints <- all.equal(fixed_peaks_ints,new_peaks_ints)
      if (test_mzs == TRUE && test_ints == TRUE) {
        test_peaks <- TRUE
      } else {
        test_peaks <- paste0("MZs = ",test_mzs, "; INTs = ",test_ints)
      }
      return(test_peaks)
    })
    # check if all peaks list matches and reduce their results, if there is
    # at least one mismatch concatenate the error messages and count the number
    # of not matching lists
    check_test_peaks_match <- (test_mgf_peaks_equal == TRUE)
    if (all(check_test_peaks_match)) {
      test_mgf_peaks <- TRUE
    } else {
      number_mismatch_peaks <- sum(!check_test_peaks_match)
      scans_mismatch <- fixed_mgf_header[idx_common_headers[which(!check_test_peaks_match)],"scans"]
      test_mgf_peaks <- paste0("A total of ", number_mismatch_peaks, 
                                     " ions peaks list mismatched, their scans number and differences are highlighted next: ", 
                                     paste0("(SCANS ",scans_mismatch,": ", 
                                            test_mgf_peaks_equal[!check_test_peaks_match], 
                                            ")", collapse=";"))
    }
  } else {
    test_mgf_peaks <- paste("No ions header match in the MGFs ", mgf_type,
                            ". MGFs peaks list comparison aborted.")
  }
  
  if (test_mgf_header_tables_equal == TRUE && test_mgf_peaks == TRUE) {
    return(TRUE)
  } else {
    return(paste0("MGF header = ", test_mgf_header_tables_equal, 
                  ";\nMGF peaks = ", test_mgf_peaks))
  }
}

## Function to compare the results of two np3 jobs from the run or the join_jobs commands - a fixed result (reference) and a new result ##
# It compares the clustering, clean and merger count tables; the MGFs from clustering and from clean;
# the similarity and matches tables; the molecular networks tables and the final report tables;
## Where fixed_result_path is the path to a fixed result of a np3 job to be used as reference -
# the job name is extracted from here;
# new_result_path is the path to a new result of the same np3 job to be checked for consistency -
# the new job name is extracted from here;
# output_path_test is the path to the output_path from the test command parameter. 
# This directory path concatenated with the sub folders '/test/compare_np3_fixed_results/' is where
# the difference tables will be stored, if any;
# sim_w_cutoff is the similarity cutoff parm of the SSMN in character as written in the output file;
# topk is the top k neigbours parm of the SSMN of the jobs;
# maxComponentSize is the maximum component size parm of the SSMN of the jobs;
# minMatchedPeaks is the minimum number of matched peaks parm of the SSMN of the jobs;
# gnps_library_search_tool is the gnps library search tool used or empty string "" if it was disabled
# tremolo_exec is a TRUE or FALSE indicating if the tremolo search tool was executed
## Results print at the screen the number of correct tests and subtests by groups
# if any mismatch is found, also print the number of failures and the
# error messages and detected differences by groups
compare_two_np3_run_results <- function(fixed_result_path, new_result_path, output_path_test,
                                        sim_w_cutoff="06",
                                        topk="15",
                                        maxComponentSize=200,
                                        minMachedPeaks=6,
                                        gnps_library_search_tool="gnps",
                                        tremolo_exec=TRUE)
    {
  cat("* Testing the equality of the new NP3 results against the fixed results *\n")
  # creates the output dir to store the test results if not present yet
  output_path_test_diff <- file.path(output_path_test, "/test/compare_np3_fixed_results/")
  if (!dir.exists(output_path_test_diff)) {
    dir.create(output_path_test_diff, recursive = TRUE)
  }
  output_path_test_diff <- normalizePath(output_path_test_diff)
  
  # check the results paths
  if (!dir.exists(fixed_result_path))
  {
    stop("The fixed NP3 result path '", fixed_result_path,
         "' does not exists. Provide a valid path to where the fixed results are located.")
  }
  if (!dir.exists(new_result_path))
  {
    stop("The new NP3 result path '", new_result_path,
         "' does not exists. Provide a valid path to where the new results are located.")
  }
  
  ## set jobs names
  job_name <- basename(fixed_result_path)
  new_job_name <- basename(new_result_path)
  
  # check the results outs paths
  if (!dir.exists(file.path(fixed_result_path, "outs", job_name)))
  {
    stop("The fixed NP3 outs result path '", file.path(fixed_result_path, "outs", job_name),
         "' does not exists. Provide a valid path to where the fixed results are located.")
  }
  if (!dir.exists(file.path(new_result_path, "outs", new_job_name)))
  {
    stop("The new NP3 outs result path '", file.path(new_result_path, "outs", new_job_name),
         "' does not exists. Provide a valid path to where the new results are located.")
  }
  
  ## Compare count tables ##
  
  count_table_names <- c("clustering", "clean", "merge")
  count_tables_suffix <- c("_peak_area.csv", "_peak_area_clean_ann.csv", "_peak_area_merged_ann.csv")
  count_table_subfolder <- c("", "clean", "merge")
  count_table_anti_join_by <- c(c("msclusterID","mzConsensus","rtMean"), 
                                c("msclusterID","mzConsensus","rtMean","joinedIDs"),
                                c("msclusterID","mzConsensus","rtMean","mergedIDs_all"))
  # exec comparison for each count table for clustering, clean and merge
  test_count_tables <- lapply(seq_along(count_table_names), function(i) {
    fixed_count_table_path <- file.path(fixed_result_path, "outs", job_name, 
                                        "count_tables", count_table_subfolder[[i]],
                                        paste0(job_name, count_tables_suffix[[i]]))
    new_count_table_path <- file.path(new_result_path, "outs", new_job_name, 
                                      "count_tables", count_table_subfolder[[i]],
                                      paste0(new_job_name, count_tables_suffix[[i]]))
    if (!file.exists(fixed_count_table_path)) {
      # if the fixed count table does not exists, this result is not expected
      test_count <- FALSE
    } else if (!file.exists(new_count_table_path)) {
      test_count <- paste("The ", count_table_names[[i]],
                          " count table of the new job was not created - missing file.")
    } else {
      test_count <- compare_tables(fixed_table_path=fixed_count_table_path, 
                                  new_table_path=new_count_table_path,
                                  table_type=paste0(count_table_names[[i]],
                                                    "_count"), 
                                  job_name=job_name, 
                                  output_path_test_diff=output_path_test_diff)
    }
    return(test_count)
  })
  names(test_count_tables) <- count_table_names
  # remove any missing data from the fixed result, which had test result = FALSE
  if (any(unlist(test_count_tables) == FALSE)) {
    test_count_tables <- test_count_tables[unlist(test_count_tables) != FALSE]
  }
  
  ## Compare MGFs ##
  
  # compare final MGFs by reading them and checking for equality in the 
  # header fields and the fragmented peaks list
  mgfs_type <- c("clustering", "clean")
  mgfs_suffix <- c("all", "clean")
  
  # exec comparison for each MGF file
  test_mgfs <- lapply(seq_along(mgfs_type), function(i) {
    fixed_mgf_path <- file.path(fixed_result_path, "outs", job_name, 
                                "mgf", paste0(job_name, "_", mgfs_suffix[[i]],".mgf"))
    new_mgf_path <- file.path(new_result_path, "outs", new_job_name, 
                              "mgf", paste0(new_job_name, "_", mgfs_suffix[[i]],".mgf"))
    if (!file.exists(new_mgf_path)) {
      test_mgfs_ <- paste("The ", mgfs_type[[i]],
                          " MGF file of the new job was not created - missing file.")
    } else {
      test_mgfs_ <- compare_mgfs(fixed_mgf_path=fixed_mgf_path, 
                                 new_mgf_path=new_mgf_path, 
                                 mgf_type=mgfs_type[[i]], job_name=job_name, 
                                 output_path_test_diff=output_path_test_diff)
    }
    return(test_mgfs_)
  })
  names(test_mgfs) <- paste0("mgf_",mgfs_type)
  
  ## Compare pairwise similarity and matches tables ##
  
  sim_table_names <- c("clustering", "clean_aggMax", "clean","clustering", "clean")
  sim_tables_type <- c("similarity", "similarity", "similarity", "matches","matches")
  sim_tables_suffix <- c(".csv", "_aggMax.csv", "_clean.csv", ".csv", "_clean.csv")
  # exec comparison for each count table for clustering, clean and merge
  test_sim_matches_tables <- lapply(seq_along(sim_table_names), function(i) {
    fixed_sim_table_path <- file.path(fixed_result_path, "outs", job_name, 
                                      "molecular_networking", "similarity_tables", 
                                      paste0(ifelse(sim_tables_type[[i]] == "similarity", 
                                                    "similarity_table_",
                                                    "similarity_table_matches_") , 
                                             job_name, sim_tables_suffix[[i]]))
    new_sim_table_path <- file.path(new_result_path, "outs", new_job_name, 
                                    "molecular_networking", "similarity_tables", 
                                    paste0(ifelse(sim_tables_type[[i]] == "similarity", 
                                                  "similarity_table_",
                                                  "similarity_table_matches_") , 
                                           new_job_name, sim_tables_suffix[[i]]))
    if (!file.exists(new_sim_table_path)) {
      test_sim_matches <- paste0("The ", sim_table_names[[i]]," pairwise ",
                                 sim_tables_type[[i]], 
                                 " table of the new job was not created - missing file:",new_sim_table_path,".")
    } else {
      test_sim_matches <- compare_tables(fixed_table_path=fixed_sim_table_path, 
                                   new_table_path=new_sim_table_path,
                                   table_type=paste0("pairwise_",sim_tables_type[[i]],"_",sim_table_names[[i]]), 
                                   job_name=job_name, 
                                   output_path_test_diff=output_path_test_diff)
    }
    return(test_sim_matches)
  })
  names(test_sim_matches_tables) <- paste0(sim_tables_type, "_", sim_table_names)
  
  
  ## Compare molecular networks tables ## 
  
  # test ssmn, ssmn filtered, ssmn protonated, ivamn, ivamn protonated and ivamn att tables
  mn_table_names <- c("SSMN", "filtered SSMN", "protonated SSMN","IVAMN", 
                      "protonated IVAMN", "IVAMN attributes")
  mn_tables_type <- c("complete_ssmn", "filtered_ssmn", "protonated_ssmn", 
                      "complete_ivamn","protonated_ivamn", "ivamn_attributes")
  mn_tables_suffix <- c(paste0("_ssmn_w_",sim_w_cutoff,".selfloop"), 
                        paste0("_ssmn_w_",sim_w_cutoff,"_mmp_", 
                               minMachedPeaks,"_k_", topk,"_x_",maxComponentSize,".selfloop"),
                        paste0("_ssmn_w_",sim_w_cutoff,"_protonated_mmp_", 
                               minMachedPeaks,"_k_", topk,"_x_",maxComponentSize,".selfloop"),
                        "_ivamn.selfloop", "_ivamn_protonated.selfloop", 
                        "_ivamn_attributes.csv")
  # exec comparison for each count table for clustering, clean and merge
  test_mn_tables <- lapply(seq_along(mn_table_names), function(i) {
    fixed_mn_path <- file.path(fixed_result_path, "outs", job_name, 
                                 "molecular_networking",
                                 paste0(job_name, mn_tables_suffix[[i]]))
    new_mn_path <- file.path(new_result_path, "outs", new_job_name, 
                               "molecular_networking", 
                               paste0(new_job_name, mn_tables_suffix[[i]]))
    if (!file.exists(new_mn_path)) {
      test_mn <- paste0("The ", mn_table_names[[i]],
                        " of the new job was not created - missing file:",
                        new_mn_path,".")
    } else {
      test_mn <- compare_tables(fixed_table_path=fixed_mn_path, 
                                  new_table_path=new_mn_path,
                                  table_type=mn_tables_type[[i]], 
                                  job_name=job_name, 
                                  output_path_test_diff=output_path_test_diff)
    }
    return(test_mn)
  })
  names(test_mn_tables) <- mn_tables_type
  
  ## compare Tremolo-UNPD and GNPS identifications complete results ##
  # test tremolo complete result
  if (tremolo_exec) {
    fixed_tremolo_path <- file.path(fixed_result_path, "outs", job_name, 
                                    "identifications", "tremolo_results.csv")
    new_tremolo_path <- file.path(new_result_path, "outs", new_job_name, 
                                  "identifications", "tremolo_results.csv")
    if (!file.exists(new_tremolo_path)) {
      test_tremolo <- "The complete Tremolo-UNPD identification table of the new job was not created - missing file." 
    } else {
      test_tremolo <- compare_tables(fixed_table_path=fixed_tremolo_path, 
                                     new_table_path=new_tremolo_path,
                                     table_type="tremolo_identifications", 
                                     job_name=job_name, 
                                     output_path_test_diff=output_path_test_diff,
                                     anti_join_by = c("msclusterID", "dbIndex", "Charge",
                                            "MQScore", "CompoundName", "mzErrorPPM",
                                            "SMILES", "LibSearchSharedPeaks", 
                                            "chemicalNames", "molecularFormula",
                                            "molecularWeight", "CAS", "PARENTMASS",
                                            "InChIKey", "NPClassifier_superclass",
                                            "ClassyFire_subclass", "NPAtlas_id"))
    }
  } else {
    # tremolo was not executed, return TRUE
    test_tremolo <- TRUE
  }
  # test GNPS result
  if (gnps_library_search_tool != "") {
    fixed_gnps_path <- file.path(fixed_result_path, "outs", job_name, 
                                    "identifications", 
                                    paste0(job_name, "_library_search_",
                                           gnps_library_search_tool,".tsv"))
    new_gnps_path <- file.path(new_result_path, "outs", new_job_name, 
                                  "identifications", 
                                  paste0(new_job_name, "_library_search_",
                                         gnps_library_search_tool,".tsv"))
    if (!file.exists(new_gnps_path)) {
      test_gnps <- "The complete GNPS identification table of the new job was not created - missing file." 
    } else {
      test_gnps <- compare_tables(fixed_table_path=fixed_gnps_path, 
                                   new_table_path=new_gnps_path,
                                   table_type="gnps_identifications", 
                                   job_name=job_name, 
                                   output_path_test_diff=output_path_test_diff, sep="\t")
    }
  } else {
    # gnps was not executed, return TRUE
    test_gnps <- TRUE
  }
  
  ## Compare final report table by equality ##
 
  # test quantification report
  fixed_quant_stats_path <- file.path(fixed_result_path, "outs", job_name, 
                               "final_reports", "quantification_report", 
                               "quantification_statistics.csv")
  new_quant_stats_path <- file.path(new_result_path, "outs", new_job_name, 
                                    "final_reports", "quantification_report", 
                                    "quantification_statistics.csv")
  if (!file.exists(new_quant_stats_path)) {
    test_quant_report <- "The quantification statistics report table of the new job was not created - missing file." 
  } else {
    test_quant_report <- compare_tables(fixed_table_path=fixed_quant_stats_path, 
                                new_table_path=new_quant_stats_path,
                                table_type="quantification_report_statistics", 
                                job_name=job_name, 
                                output_path_test_diff=output_path_test_diff)
  }
  
  # test chemical reports
  if (gnps_library_search_tool != "") {
    chemical_reports_filenames <- c("chemical_statistics_UNPD.csv",
                                    "chemical_identification_statistics_GNPS.csv",
                                    "chemical_identification_statistics_UNPDxGNPS_best_origin.csv")
    chemical_reports_names <- c("UNPD",
                                "GNPS",
                                "UNPDxGNPS")
  } else {
    chemical_reports_filenames <- c("chemical_statistics_UNPD.csv")
    chemical_reports_names <- c("UNPD")
  }
  # exec comparison for each chemical report table
  test_chem_reports <- lapply(seq_along(chemical_reports_names), function(i) {
    fixed_report_path <- file.path(fixed_result_path, "outs", job_name, 
                                   "final_reports", "chemical_report",
                                   chemical_reports_filenames[[i]])
    new_report_path <- file.path(new_result_path, "outs", new_job_name,
                                 "final_reports", "chemical_report",
                                 chemical_reports_filenames[[i]])
    if (!file.exists(new_report_path)) {
      test_report_stats <- paste("The ", chemical_reports_names[[i]],
                                 " chemical report table of the new job was not created - missing file.")
    } else {
      test_report_stats <- compare_tables(fixed_table_path=fixed_report_path, 
                                          new_table_path=new_report_path,
                                          table_type=paste0("chem_report_",
                                                            chemical_reports_names[[i]]), 
                                          job_name=job_name, 
                                          output_path_test_diff=output_path_test_diff)
    }
    return(test_report_stats)
  })
  names(test_chem_reports) <- chemical_reports_names
  
  # test molecular_networking reports
  mn_reports_suffix <- c("ivamn", "ivamn_protonated", 
                         paste0("ssmn_w_",sim_w_cutoff), 
                         paste0("ssmn_w_",sim_w_cutoff,"_mmp_", 
                                minMachedPeaks,"_k_", topk,"_x_",maxComponentSize),
                         paste0("ssmn_w_",sim_w_cutoff,"_protonated_mmp_", 
                                minMachedPeaks,"_k_", topk,"_x_",maxComponentSize))
  mn_reports_names <- c("ivamn", "ivamn_protonated", "ssmn", "ssmn_filtered", 
                        "ssmn_protonated")
  # exec comparison for each report table
  test_mn_reports <- lapply(seq_along(mn_reports_names), function(i) {
    fixed_report_path <- file.path(fixed_result_path, "outs", job_name, 
                                 "final_reports", "molecular_networking_report",
                                 paste0("molecular_networking_statistics_", 
                                        job_name, "_", 
                                        mn_reports_suffix[[i]],".csv"))
    new_report_path <- file.path(new_result_path, "outs", new_job_name,
                                 "final_reports", "molecular_networking_report",
                                 paste0("molecular_networking_statistics_", 
                                        new_job_name, "_", 
                                        mn_reports_suffix[[i]],".csv"))
    if (!file.exists(fixed_report_path)) {
      # if the fixed report table does not exists, this result is not expected
      test_report_stats <- FALSE
    } else if (!file.exists(new_report_path)) {
      test_report_stats <- paste("The ", mn_reports_names[[i]],
                         " molecular networking report table of the new job was not created - missing file: ",
                         new_report_path,".")
    } else {
      test_report_stats <- compare_tables(fixed_table_path=fixed_report_path, 
                                  new_table_path=new_report_path,
                                  table_type=paste0("mn_report_",
                                                    mn_reports_names[[i]]), 
                                  job_name=job_name, 
                                  output_path_test_diff=output_path_test_diff)
    }
    return(test_report_stats)
  })
  names(test_mn_reports) <- mn_reports_names
  if (any(unlist(test_mn_reports) == FALSE)) {
    test_mn_reports <- test_mn_reports[unlist(test_mn_reports) != FALSE]
  }
  ## Reduce all the tests results and print the error messages with the differences if any, 
  # using the tag ERROR if any mismatch was found
  # test_count_tables : contains list of 3 results for "clustering","clean","merge"
  # test_mgfs : contains list of 2 results for "mgf_clustering", "mgf_clean"
  # test_sim_matches_tables : contains list of 5 results for similarity_clustering,similarity_clean_aggMax,similarity_clean,matches_clustering,matches_clean
  # test_mn_tables : contains list of 6 results for "complete_ssmn,filtered_ssmn,protonated_ssmn,complete_ivamn,protonated_ivamn,ivamn_attributes"
  # test_tremolo : contains 1 result 
  # test_gnps : contains 1 result  
  # test_quant_report : contains 1 result
  # test_chem_reports : contains list of 3 results for "UNPD,GNPS,UNPDxGNPS"
  # test_mn_reports : contains list of 5 results for "ivamn,ivamn_protonated,ssmn,ssmn_filtered,ssmn_protonated" 
  list_tests <- list(Count_tables = test_count_tables, MGFs=test_mgfs, 
                     Similarity_matches_tables=test_sim_matches_tables,
                     Molecular_network_tables=test_mn_tables, 
                     Tremolo_identifications=test_tremolo, 
                     GNPS_identifications=test_gnps, 
                     Quantification_report=test_quant_report,
                     Chemical_reports=test_chem_reports, 
                     Molecular_networking_reports=test_mn_reports)
  number_tests <- length(list_tests)
  number_subtests <- sum(sapply(list_tests, length))
  #print(list_tests)
  # check the correct groups of tests and all subtests
  correct_tests <- sapply(list_tests, function (x) all(unlist(x) == TRUE))
  correct_subtests <- unlist(sapply(list_tests, function (x) (unlist(x) == TRUE)))
  cat("\nCorrect equality tests by groups = ", 
      sum(correct_tests),"/",number_tests)
  cat("\nCorrect equality subtests within groups = ", 
      sum(correct_subtests),"/",number_subtests)
  
  if (sum(correct_tests) == number_tests && sum(correct_subtests) == number_subtests) {
    cat("\nDone! :)\n")
  } else {
    cat("\n\nTotal tests failures:", sum(!correct_tests), ":(",
        "\nTotal subtests failures:", sum(!correct_subtests), ":(\n",
        "\nThe following ERRORS were detected from mismatching results:\n")
    
    for (i in which(!correct_tests)) {
      cat("\n- Error in the *",names(list_tests[i]), "* group, check the detected mismatches:\n")
      if (length(list_tests[[i]]) > 1) {
        failed_subtests <- which(unlist(list_tests[[i]]) != TRUE)
        for (j in failed_subtests) {
          cat("    - Difference in the subgroup *", names(list_tests[[i]][j]), "* =", list_tests[[i]][[j]], "\n")
        }
      } else {
        cat("    - Difference =", list_tests[[i]][[1]], "\n")
      }
    }
  }
}


# make a comparison for pre process command result between a fixed and a new job based on the 
# pre processed MGFs and MS1 tables
# compare the pre processed MGFs of all sample codes present in the fixed metadata
# and also compare the tables created in the pre process:"MS1_no_MS2", "MS1_with_MS2", "log_MS2_no_MS1_fake_peaks"
## Where job_name is the name of the job from which this pre process is from;
# metadata_fixed_path is the path to the fixed metadata used to compare the pre process results;
# fixed_pp_result_path is the path to a fixed result of a np3 pre process job to be used as reference;
# new_pp_result_path is the path to a new result of the same np3 pre process job to be checked for consistency;
# output_path_test is the path to the output_path from the test command parameter. 
# This directory path concatenated with the sub folders '/test/compare_np3_fixed_results/' is where
# the difference tables will be stored, if any;
## Results print at the screen the number of correct tests and subtests by groups
# if any mismatch is found, also print the number of failures and the
# error messages and detected differences by groups
compare_two_np3_preprocess_results <- function(job_name,
                                               metadata_fixed_path,
                                               fixed_pp_result_path, 
                                               new_pp_result_path, 
                                               output_path_test) 
{
  cat("* Testing the equality of the new NP3 pre process results against the fixed results *\n")
  # creates the output dir to store the test results if not present yet
  output_path_test_diff <- file.path(output_path_test, "/test/compare_np3_fixed_results/")
  if (!dir.exists(output_path_test_diff)) {
    dir.create(output_path_test_diff, recursive = TRUE)
  }
  output_path_test_diff <- normalizePath(output_path_test_diff)
  
  # for each pre processed mgf present in the metadata, 
  # compare the fixed with the new result
  metadata_fixed <- readMetadataTable(metadata_fixed_path)
  
  # compare the pre processed MGFs by reading them and checking for equality in the 
  # header fields and the fragmented peaks list
  mgfs_type <- c("pre_process")
  mgfs_sample_code_prefix <- metadata_fixed$SAMPLE_CODE
  
  # exec comparison for each pp MGF file
  test_pp_mgfs <- lapply(seq_along(mgfs_sample_code_prefix), function(i) {
    fixed_mgf_path <- file.path(fixed_pp_result_path,  
                                paste0(mgfs_sample_code_prefix[[i]], "_peak_info.mgf"))
    new_mgf_path <- file.path(new_pp_result_path, 
                              paste0(mgfs_sample_code_prefix[[i]], "_peak_info.mgf"))
    if (!file.exists(new_mgf_path)) {
      test_mgfs_ <- paste("The ", mgfs_type,
                          " MGF file of sample code ",
                          mgfs_sample_code_prefix[[i]],
                          " of the new job was not created - missing file:", 
                          new_mgf_path,".")
    } else {
      test_mgfs_ <- compare_mgfs(fixed_mgf_path=fixed_mgf_path, 
                                 new_mgf_path=new_mgf_path, 
                                 mgf_type=paste0(mgfs_type,"_",mgfs_sample_code_prefix[[i]]), 
                                 job_name=job_name, 
                                 output_path_test_diff=output_path_test_diff)
    }
    return(test_mgfs_)
  })
  names(test_pp_mgfs) <- paste0("mgf_",mgfs_type,"_",mgfs_sample_code_prefix)
  
  # compare the pre process MS1 table lists and MS2 no match log table
  ms1_tables_filenames <- c("MS1_list_no_MS2.csv", "MS1_list_with_MS2.csv",
                            "log_MS2_no_MS1peak_match.csv")
  ms1_tables_names <- c("MS1_no_MS2", "MS1_with_MS2", "MS2_no_MS1_fake_peaks")
  # exec comparison for each pre process table
  test_pp_tables <- lapply(seq_along(ms1_tables_names), function(i) {
    fixed_pp_table_path <- file.path(fixed_pp_result_path, 
                                     ms1_tables_filenames[[i]])
    new_pp_table_path <- file.path(new_pp_result_path, 
                                   ms1_tables_filenames[[i]])
    if (!file.exists(new_pp_table_path)) {
      test_pp_table <- paste("The ", ms1_tables_names[[i]],
                             " pre processing table of the new job was not created - missing file: ",
                             new_pp_table_path,".")
    } else {
      test_pp_table <- compare_tables(fixed_table_path=fixed_pp_table_path, 
                                          new_table_path=new_pp_table_path,
                                          table_type=paste0("pp_table_",
                                                            ms1_tables_names[[i]]), 
                                          job_name=job_name, 
                                          output_path_test_diff=output_path_test_diff)
    }
    return(test_pp_table)
  })
  names(test_pp_tables) <- paste0("pre_process_table_",ms1_tables_names)
  
  ## Reduce all the tests results and print the error messages with the differences if any, 
  # using the tag ERROR if any mismatch was found
  # test_pp_mgfs : contains list of n results, where n is the number of rows in the fixed metadata
  # MS1_tables : contains list of 3 results for MS1_no_MS2 and MS1_w_MS2 and log_MS2_no_MS1
  list_tests <- list(Pre_process_mgfs = test_pp_mgfs, 
                     Pre_process_tables=test_pp_tables)
  number_tests <- length(list_tests)
  number_subtests <- sum(sapply(list_tests, length))
  # check the correct groups of tests and all subtests
  correct_tests <- sapply(list_tests, function (x) all(unlist(x) == TRUE))
  correct_subtests <- unlist(sapply(list_tests, function (x) (unlist(x) == TRUE)))
  cat("\nCorrect equality tests by groups = ", 
      sum(correct_tests),"/",number_tests)
  cat("\nCorrect equality subtests within groups = ", 
      sum(correct_subtests),"/",number_subtests)
  
  if (sum(correct_tests) == number_tests) {
    cat("\nDone! :)\n")
  } else {
    cat("\n\nTotal tests failures:", sum(!correct_tests), ":(",
        "\nTotal subtests failures:", sum(!correct_subtests), ":(\n",
        "\nThe following ERRORS were detected from mismatching results:\n")
    
    for (i in which(!correct_tests)) {
      cat("\n- Error in the *",names(list_tests[i]), "* group, check the detected mismatches:\n")
      if (length(list_tests[[i]]) > 1) {
        failed_subtests <- which(unlist(list_tests[[i]]) != TRUE)
        for (j in failed_subtests) {
          cat("    - Difference in the subgroup *", names(list_tests[[i]][j]), "* =", list_tests[[i]][[j]], "\n")
        }
      } else {
        cat("    - Difference =", list_tests[[i]][[1]], "\n")
      }
    }
  }
}

