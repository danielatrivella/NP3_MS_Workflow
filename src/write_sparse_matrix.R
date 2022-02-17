#' @param x A sparse matrix from the Matrix package.
#' @param file A filename that ends in ".gz".
#' @param total_written_spec current total number of written spectra, to correctly number the rows.
writeMMgz_upperTri <- function(x, file, total_written_spec_x, total_written_spec_y) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  n_written_values <- length(x@x)
  x <- summary(x)
  # skip written spectra
  x[,1] <- x[,1] + total_written_spec_x
  x[,2] <- x[,2] + total_written_spec_y
  # convert to a upper triangular matrix, invert i and j (x is symmetric)
  names(x) <- c("j","i","x")
  x <- x[, c("i","j","x")]
  
  readr::write_delim(
    x = x,
    path = file,
    append = TRUE,
    delim = " ",
    col_names = FALSE
  )
  # return the number of not 0 values
  n_written_values
}

#' @param x A sparse matrix from the Matrix package.
#' @param file A filename that ends in ".gz".
#' @param mnrows number of rows of the sparse matrix
#' @param mncols number of cols of the sparse matrix
#' @param nvalues number of non zero values in the sparse matrix
writeMMgzHeader <- function(x, file, mnrows, mncols, nvalues) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
      sprintf("%s %s %s", mnrows, mncols, nvalues)
    ),
    gzfile(file)
  )
}

# TODO make read function to read sparse matrix and the header and assign rownames and colnames

# For converting base matrix to sparse matrix
# sparsematrix <- as(BaseMatrix, "sparseMatrix")
# as(as.matrix(bind_cols(comp_row_sim_matches[1,], list(c(rep(0.000,n_scans-1), 1.000)))), "sparseMatrix")

writeMMgz(x=as(as.matrix(bind_cols(comp_row_sim_matches[1,], list(c(rep(0.000,n_scans-1), 1.000)))), "sparseMatrix"),
          file=file.path(output_path, 
                    paste0("similarity_table_", data_name, "_tmp.csv")))


writeMMgz(x=as(as.matrix(bind_cols(comp_row_sim_matches[2,],
                                   list(c(rep(0.000,n_scans-1), length(ms2_sample$MZS[[n_scans]]))))), 
          "sparseMatrix"),
     file = file.path(output_path, 
                      paste0("similarity_table_matches_", data_name, "_tmp.csv")))