#include <RcppArmadillo.h>

using namespace std;

// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(RcppArmadillo)]]

// make the lower triangular matrix between cols [m,n) equals zero
// [[Rcpp::export]]
void zero_matrix_tri_down(Rcpp::NumericMatrix x_, int m, int n) {
  arma::mat M( x_.begin(), x_.nrow(), x_.ncol(), false ) ;
  // indices start from 0 
  m = m - 1;
  //n = n - 1 ;
  for (int i = 0; i < n - m; i++) {
    for (int j = i + 1; j < n - m; j++) {
      M(j, i + m) = 0.0;
    }
  }
}

// mirror an upper matrix between columns [m,n)
// [[Rcpp::export]]
void mirror_matrix_tri_upper(Rcpp::NumericMatrix x_, int m, int n) {
  arma::mat M( x_.begin(), x_.nrow(), x_.ncol(), false ) ;
  // indices start from 0 
  m = m - 1;
  //n = n - 1 ;
  for (int i = 0; i < n - m; i ++) {
    for (int j = i + 1; j < n - m; j++) {
      M(j, i + m) = M(i, j + m);
    }
  }
}

// [[Rcpp::export]]
vector<int> which_ge(vector<double>  x, double value, int gap) {
  vector<int> idx;
  
    for (unsigned int i = 0; i < x.size(); i++)
  {
    if (x[i] >= value)
      idx.push_back(i+1+gap);
  }
  
  return(idx);
}

// [[Rcpp::export]]
vector<int> which_eq(vector<double>  x, double value, int stop) {
  vector<int> idx;
  
  for (unsigned int i = 0; i < x.size(); i++)
  {
    if (x[i] == value) 
    {
      idx.push_back(i+1);
      if (stop)
        break;
    }
  }
  
  return(idx);
}

// [[Rcpp::export]]
vector<int> which_not_NA(vector<string>  x) {
  vector<int> idx;
  
  for (unsigned int i = 0; i < x.size(); i++)
  {
    //Rcpp::Rcout << x[i] << endl;
    if (x[i] != "NA")
      idx.push_back(i+1);
  }
  
  return(idx);
}


// // [[Rcpp::export]]
// void transpose_matrix(double *src, double *dst, const int N, const int M) {
//   
//   for(int n = 0; n<N*M; n++) {
//     int i = n/N;
//     int j = n%N;
//     dst[n] = src[M*j + i];
//   }
// }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// # /*** R
// # # timesTwo(data.frame(c(1,2), c(3,4)), 2)
// # # a <- matrix(1:25000000, nrow = 5000)
// # a <- matrix(as.double(1:15), 3, 5)
// # a
// # zero_matrix_tri_down(a, 2, 4)
// # a
// # mirror_matrix_tri_upper(a, 2, 4)
// # a
// # */