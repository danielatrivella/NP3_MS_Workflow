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
                                   c(1.0,5))
# test no intersection
test_no_shift$all_diff <- all.equal(
  c(normDotProductShift(c(3), c(5),c(1,2), c(1,2), 0.025, 0),
    normDotProductShift(c(3), c(5),c(1,2), c(1,2), 0.025, 0)),
  c(normDotProductShift(c(3,4,5), c(3,4,5),c(1,2), c(1,2), 0.025, 0), c(0,0)))
# test one intersection  no shift
test_no_shift$one_equal <- all.equal(c(normDotProductShift(c(1,2,3), c(1,2,3), c(3), c(3),0.025, 0),
                                       normDotProductShift(c(1,2,3), c(1,2,3), c(3), c(3),0.025, 0)),
                                     c(normDotProductShift(c(3), c(3),c(1,2,3), c(1,2,3), 0.025, 0),
                                      c(3*3/sqrt((1+2*2+3*3)*(3*3)),1)))

test_w_shift <- list()
# test equal with shift in the end
test_w_shift$all_equal_shiftEnd <- all.equal(normDotProductShift(c(1,2,3,4,5), c(1,2,3,4,5), c(1,2,3,4,6), c(1,2,3,4,5), 0.0, -1.0),
                                            c(1.0,5))
# test equal with shift in the end and in the beginning
test_w_shift$all_equal_shiftBeginEnd <- all.equal(normDotProductShift(c(1,3,4,5,7), c(1,2,3,4,5), c(2,3,4,5,8), c(1,2,3,4,5), 0.0, -1.0), 
                                                 c(1.0,5))
# middle intersection with shift
test_w_shift$middle_intersection <- all.equal(normDotProductShift(c(1,2,3,4), c(1,2,3,4), c(2.5,3), c(2,3), 0.025, -0.5),
                                             c((2*2+3*3)/sqrt((2*2+3*3)*(1+2*2+3*3+4*4)),2))
# test all equal with shift 
test_w_shift$all_equal <- all.equal(c(normDotProductShift(c(3,6), c(3,2), c(5,8), c(3,2), 0.025, -2),
                                      normDotProductShift(c(3,6), c(3,2), c(5,8), c(3,2), 0.025, -2)),
                                    c(c(1,2),normDotProductShift(c(3,6), c(3,2), c(3,6), c(3,2), 0.025, 0)))
# test all shifted
test_w_shift$all_equal2 <- all.equal(normDotProductShift(c(3,6), c(5,3), c(1,4), c(5,3), 0.025, 2),
                                     c(1.0,2))
# test equality and inverted peaks with shift
test_w_shift$inverted_ABBA <- all.equal(c(normDotProductShift(c(0,4,5,7), c(7,5,3,1), c(1,2,5,7), c(7,5,3,1), 0.025, 1),
                                          normDotProductShift(c(1,2,5,7), c(7,5,3,1),c(0,4,5,7), c(7,5,3,1), 0.025, -1)),
                                       c(c(10/(7*7+5*5+9+1),2),c(10/(7*7+5*5+9+1),2)))
# test inverted peaks A,B e B,A with shift
test_w_shift$inverted_equal <- all.equal(normDotProductShift(c(3,8,10), c(3,2,1), c(3, 7,9), c(3,2,1),0.025, 1),
                                         normDotProductShift(c(3, 7,9), c(3,2,1), c(3,8,10), c(3,2,1), 0.025, -1))
# test match equal first without the shift
test_w_shift$match_first_no_shift <- all.equal(c(normDotProductShift(c(5,6), c(3,5), c(5,7), c(3,5), 0.025, 2),
                                                 normDotProductShift(c(5,6), c(3,5), c(5,7), c(3,5), 0.025, 2)),
                                               c(c((3*3)/sqrt((3*3+5*5)*(3*3+5*5)),1),
                                                 normDotProductShift(c(3,6), c(3,5), c(3,5), c(3,5), 0.025, -2)))

# test match spec B fully intersect spec A, but spec A have more peaks -> sim should not be 1.0
test_w_shift$match_full_intersect_B <- all.equal(normDotProductShift(
  c(308.102234), c(31.62278),
   c(269.145386, 287.221497, 288.2575687), c(3.87944576234677,29.1911304206282,11.525094591494),
  0.005,
  326.1-287.14),
  c(31.62278*3.87944576234677/sqrt(31.62278^2*sum( c(3.87944576234677,29.1911304206282,11.525094591494)^2)),1))

# test match 2 shifted m/zs peak (first and third) in spec B with one not matched peak (second) between them
test_w_shift$match_2_shifts_spaced <- all.equal(normDotProductShift(c(77.038330078125,79.054443359375,81.0339050292969,89.0368041992188,91.0538482666016,92.0575714111328,95.0485000610352,105.04426574707,106.048477172852,107.04906463623,108.05281829834,118.039901733398,119.049049377441,129.04150390625,135.043365478516,163.038772583008,164.042190551758), 
                                                                    c(6.54952979244544,3.9730417514932,2.4677802183366,2.28829201282441,4.40948950078457,1.78493015786444,7.43108582222807,13.0018173753784,2.40045724834929,17.3090629174236,3.23141965622511,1.8187707210105,6.82311363182329,1.89471345234566,4.74020167336629,16.6688514848295,3.67610117832939),
                                                                    c(56.0500030517578,94.0288162231445,112.039352416992,112.049842834473),
                                                                    c(3.18129157686205,8.31729792964629,30.0262137313114,4.3736058366176),
                                                                    0.005, 163.039-112.039),
                                               c((17.309062917423*3.1812915768620+16.668851484829*30.026213731311)/1000,2))

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