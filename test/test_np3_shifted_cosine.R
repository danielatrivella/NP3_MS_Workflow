##
# functionality test for the NP3 shifted cosine similarity function results
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
Rcpp::sourceCpp(file.path(script_path(),
                          '../src/norm_dot_product.cpp'))

cat("\n*** Testing the NP3 shifted cosine similarity function ***\n")

test_no_shift <- list()
# test equal no shift
test_no_shift$all_equal <- all.equal(normDotProductShift(c(1,2,3,4,5), c(1,2,3,4,5), 
                                                         c(1,2,3,4,5), c(1,2,3,4,5), 0.0, 1.0),
                                    1.0)
# test no intersection
test_no_shift$all_diff <- all.equal(
  c(normDotProductShift(c(3), c(5),c(1,2), c(1,2), 0.025, 0),normDotProductShift(c(3), c(5),c(1,2), c(1,2), 0.025, 0)),
  c(normDotProductShift(c(3,4,5), c(3,4,5),c(1,2), c(1,2), 0.025, 0), 0))
# test one intersection  no shift
test_no_shift$one_equal <- all.equal(c(normDotProductShift(c(1,2,3), c(1,2,3), c(3), c(3),0.025, 0),
                                       normDotProductShift(c(1,2,3), c(1,2,3), c(3), c(3),0.025, 0)),
                                     c(normDotProductShift(c(3), c(3),c(1,2,3), c(1,2,3), 0.025, 0),
                                       3*3/sqrt((1+2*2+3*3)*(3*3))))

test_w_shift <- list()
# test equal with shift in the end
test_w_shift$all_equal_shiftEnd <- all.equal(normDotProductShift(c(1,2,3,4,5), c(1,2,3,4,5), c(1,2,3,4,6), c(1,2,3,4,5), 0.0, -1.0),
                                             1.0)
# test equal with shift in the end and in the beginning
test_w_shift$all_equal_shiftBeginEnd <- all.equal(normDotProductShift(c(1,3,4,5,7), c(1,2,3,4,5), c(2,3,4,5,8), c(1,2,3,4,5), 0.0, -1.0), 
                                                  1.0)
# middle intersection with shift
test_w_shift$middle_intersection <- all.equal(normDotProductShift(c(1,2,3,4), c(1,2,3,4), c(2.5,3), c(2,3), 0.025, -0.5),
                                              (2*2+3*3)/sqrt((2*2+3*3)*(1+2*2+3*3+4*4)))
# test all equal with shift 
test_w_shift$all_equal <- all.equal(c(normDotProductShift(c(3,6), c(3,2), c(5,8), c(3,2), 0.025, -2),
                                      normDotProductShift(c(3,6), c(3,2), c(5,8), c(3,2), 0.025, -2)),
                                    c(1,normDotProductShift(c(3,6), c(3,2), c(3,6), c(3,2), 0.025, 0)))
# test all shifted
test_w_shift$all_equal2 <-all.equal(normDotProductShift(c(3,6), c(5,3), c(1,4), c(5,3), 0.025, 2),
                                    1.0)
# test equality and inverted peaks with shift
test_w_shift$inverted_ABBA <- all.equal(c(normDotProductShift(c(0,4,5,7), c(7,5,3,1), c(1,2,5,7), c(7,5,3,1), 0.025, 1),
                                          normDotProductShift(c(1,2,5,7), c(7,5,3,1),c(0,4,5,7), c(7,5,3,1), 0.025, -1)),
                                        c(10/(7*7+5*5+9+1),10/(7*7+5*5+9+1)))
# test inverted peaks A,B e B,A with shift
test_w_shift$inverted_equal <- all.equal(normDotProductShift(c(3,8,10), c(3,2,1), c(3, 7,9), c(3,2,1),0.025, 1),
                                         normDotProductShift(c(3, 7,9), c(3,2,1), c(3,8,10), c(3,2,1), 0.025, -1))
# test match equal first without the shift
test_w_shift$match_first_no_shift <- all.equal(c(normDotProductShift(c(5,6), c(3,5), c(5,7), c(3,5), 0.025, 2),
                                                 normDotProductShift(c(5,6), c(3,5), c(5,7), c(3,5), 0.025, 2)),
                                               c((3*3)/sqrt((3*3+5*5)*(3*3+5*5)),
                                                 normDotProductShift(c(3,6), c(3,5), c(3,5), c(3,5), 0.025, -2)))

test_no_shift <- unlist(test_no_shift)
cat("\nCorrect tests without shift = ", 
    sum(test_no_shift == TRUE, na.rm = TRUE),"/",length(test_no_shift))
if (!all(test_no_shift == TRUE)) {
  cat("\nTest ERRORS: \n", paste(names(test_no_shift)[test_no_shift != TRUE], 
                               test_no_shift[test_no_shift != TRUE], sep = ' - ',
                               collapse='\n'), sep = '')
}
test_w_shift <- unlist(test_w_shift)
cat("\nCorrect tests with shift = ", 
    sum(test_w_shift == TRUE, na.rm = TRUE),"/",length(test_w_shift))
if (!all(test_w_shift == TRUE)) {
  cat("\nTest ERRORS: \n", paste(names(test_w_shift)[test_w_shift != TRUE], 
                               test_w_shift[test_w_shift != TRUE], sep = ' - ',
                               collapse='\n'), sep = '')
}

if (all(test_no_shift == TRUE) && all(test_w_shift == TRUE)) {
  cat("\nDone! :)\n")
} else {
  cat("\nTotal tests failures:", sum(c((test_no_shift != TRUE),(test_w_shift != TRUE))), 
      "of", sum(length(test_no_shift),length(test_w_shift)), ":(\n")
}