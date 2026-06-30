##
# run the equality test for the NP3 run or join_jobs commands results
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
fixed_result_path <- "./test/L754_bacs/fixed_results/L754_bacs_all"
new_result_path <- "./test/L754_bacs/L754_bacs_all"
sim_w_cutoff="06"
topk="15"
maxComponentSize=200
minMachedPeaks=6
gnps_library_search_tool="gnps"
tremolo_exec=TRUE
##

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 9) {
  stop("Nine arguments must be supplied to test the run or the join_jobs commands results for equality consistency:\n", 
       " 1 - output_path_test: the path to the output_path from the test command parameter. This directory path concatenated with the sub folders '/test/compare_np3_fixed_results/' is where the difference tables will be stored, if any;\n",
       " 2 - fixed_result_path: the path to a fixed result of a np3 run or join_jobs result to be used as reference - the job name is extracted from here;\n", 
       " 3 - new_result_path: the path to a new result of the same np3 job to be checked for consistency - the new job name is extracted from here;\n", 
       " 4 - sim_w_cutoff: the similarity cutoff parm of the SSMN in character as written in the output file;\n",
       " 5 - topk: the top k neigbours parm of the SSMN of the jobs;\n",
       " 6 - maxComponentSize: the maximum component size parm of the SSMN of the jobs;\n",
       " 7 - minMachedPeaks: the minimum number of matched peaks parm of the SSMN of the jobs;\n",
       " 8 - gnps_library_search_tool: the name of the gnps library search tool used or an empty string '' if it was disabled;\n",
       " 9 - tremolo_exec: a boolean TRUE or FALSE indicating if the tremolo search tool was executed;\n",
       call.=FALSE)
} else {
  #print(args)
  output_path_test <- file.path(args[[1]])
  fixed_result_path <- normalizePath(file.path(args[[2]]))
  new_result_path <- normalizePath(file.path(args[[3]]))
  
  sim_w_cutoff <- as.character(args[[4]])
  topk <- as.numeric(args[[5]])
  maxComponentSize <- as.numeric(args[[6]])
  minMachedPeaks <- as.numeric(args[[7]])
  gnps_library_search_tool <- as.character(args[[8]])
  tremolo_exec <- as.logical(args[[9]])
  if (is.na(tremolo_exec)) {
    stop("The tremolo_exec must be a TRUE or FALSE value indicating if the ",
         "tremolo tools was executed. A wrong value was provided '",
         args[[9]],"' and converted to NA.")
  }
}

## runs the function to compare the results of two np3 jobs from run or join_jobs commands
# - compares a fixed result (reference) and a new result ##
# at the end it prints the testing results and the errors found, if any
compare_two_np3_run_results(fixed_result_path=fixed_result_path, 
                            new_result_path=new_result_path, 
                            output_path_test=output_path_test,
                            sim_w_cutoff=sim_w_cutoff,
                            topk=topk,
                            maxComponentSize=maxComponentSize,
                            minMachedPeaks=minMachedPeaks,
                            gnps_library_search_tool=gnps_library_search_tool,
                            tremolo_exec=tremolo_exec)
