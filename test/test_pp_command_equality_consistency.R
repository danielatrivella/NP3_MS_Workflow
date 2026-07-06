##
# run the equality test for the NP3 pre_process command results
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

source(file.path(script_path(),"test_np3_compare_fixed_results.R"))

## args
output_path_test <- "/tmp/np3_ms_workflow/"
job_name <- "L754_bacs_all"
fixed_pp_result_path <- "./test/L754_bacs/fixed_results/preprocess/processed_data_all/"
new_pp_result_path <- "./test/L754_bacs/mzxml/processed_data/"
metadata_fixed_path <- "./test/L754_bacs/marine_bacteria_library_L754_metadata.csv"
##

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  cat(paste0("Five arguments must be supplied to test the pre_process commands results for equality consistency:\n", 
             " 1 - job_name: is the name of the job from which this pre process results are from;\n",
             " 2 - metadata_fixed_path: is the path to the fixed metadata used to compare the pre process results;\n",
             " 3 - fixed_pp_result_path: is the path to a fixed result of a np3 pre process job to be used as reference;\n",
             " 4 - new_pp_result_path: is the path to a new result of the same np3 pre process job to be checked for consistency;\n",
             " 5 - output_path_test: the path to the output_path from the test command parameter. This directory path concatenated with the sub folders '/test/compare_np3_fixed_results/' is where the difference tables will be stored, if any.\n"))
  stop("ERROR",call.=FALSE)
} else {
  #print(args)
  job_name <- args[[1]]
  metadata_fixed_path <- normalizePath(file.path(args[[2]]))
  fixed_pp_result_path <- normalizePath(file.path(args[[3]]))
  new_pp_result_path <- normalizePath(file.path(args[[4]]))
  output_path_test <- file.path(args[[5]])
}

## runs the function to compare the results of two np3 jobs from run or join_jobs commands
# - compares a fixed result (reference) and a new result ##
# at the end it prints the testing results and the errors found, if any
compare_two_np3_preprocess_results(job_name=job_name,
                                   metadata_fixed_path=metadata_fixed_path,
                                   fixed_pp_result_path=fixed_pp_result_path, 
                                   new_pp_result_path=new_pp_result_path, 
                                   output_path_test=output_path_test)
