#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat pmat(double x) {
  arma::mat Ans(x,x);
  Ans.fill(-1/x);
  Ans.diag() += 1;
  return Ans;
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat projectionmatrix(arma::vec x) {
  int lenx = x.size();
  arma::mat Ans = arma::zeros(x(lenx - 1),x(lenx - 1));
  Ans.submat(0, 0, x(0) - 1, x(0) - 1) = pmat(x(0));
  for(int i=1; i<lenx; i++){
    Ans.submat(x(i - 1), x(i - 1), x(i) - 1, x(i) - 1) = pmat(x(i) - x(i - 1));
  }
  return Ans;
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sumdatamatrix(arma::mat X, int N) {
  int C = X.n_cols / N;
  arma::mat Ans = X.cols(0, C-1);
  for(int i=2; i<N+1; i++){
    Ans += X.cols((i - 1) * C, i * C - 1);
  }
  return Ans;
}




#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat centerdatamatrix(arma::mat X, int N) {
  int C = X.n_cols/N;
  arma::mat Ans = X;
  arma::mat X_mean = sumdatamatrix(X, N)/N;
  for(int i=1; i<N+1; i++){
    Ans.cols((i - 1) * C, i * C - 1) -= X_mean;
  }
  return Ans;
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat transposedatamatrix(arma::mat X, int N) {
  int C = X.n_cols/N;
  int R = X.n_rows;
  arma::mat Ans(C, R * N);
  for(int i=1; i<N+1; i++){
    Ans.cols((i - 1) * R, i * R - 1) = X.cols((i - 1) * C, i * C - 1).t();
  }
  return Ans;
}



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat tcrossprodcpp(arma::mat X) {
  arma::mat ans = X * trans(X) ;  
  return ans;  
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat tcrossprod2cpp(arma::mat X, arma::mat Y) {
  arma::mat ans = X * trans(Y) ;  
  return ans;  
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat crossprodcpp(arma::mat X) {
  arma::mat ans = trans(X) * X ;  
  return ans;  
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat crossprod2cpp(arma::mat X, arma::mat Y) {
  arma::mat ans = trans(X) * Y ;  
  return ans;  
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sampleSigmaR(arma::mat X, int N) {
  int C = X.n_cols / N;
  int R = X.n_rows; 
  arma::mat Ans = arma::zeros(R, R);
  for(int i=1; i<N+1; i++){
    Ans += X.cols((i - 1) * C, i * C - 1);
  }
  return Ans;
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec statistics_centered(arma::mat X, double N) {
  int C = X.n_cols / N;
  int R = X.n_rows;
  arma::mat X_i_tcrossprod = arma::zeros(R, R);
  arma::mat X_sq_sum = arma::zeros(R, N);
  double sum_X_i_had_X_j = 0;
  double Y2N = 0;
  arma::mat X_i = X.cols((N-1) * C, N * C - 1);
  arma::mat X_tcrossprod = tcrossprodcpp(X_i);
  arma::mat X_sq = pow(X, 2);
  X_sq_sum.col(N - 1) = sum(X_sq.cols((N - 1) * C, N * C - 1), 1);
  for(int i=N-1; i>0; i--){
    X_i = X.cols((i - 1) * C, i * C - 1);
    X_i_tcrossprod = tcrossprodcpp(X_i);
    Y2N += accu(X_i_tcrossprod % X_tcrossprod);
    X_tcrossprod += X_i_tcrossprod;
    X_sq_sum.col(i - 1) = sum(X_sq.cols((i - 1) * C, i * C - 1), 1);
    for(int j=i+1; j<N+1; j++){
      sum_X_i_had_X_j += pow(accu(X_i % X.cols((j - 1) * C, j * C - 1)), 2);
    }
  }
  arma::mat X_sq_sum_cum = fliplr(X_sq_sum);
  X_sq_sum_cum = cumsum(X_sq_sum_cum, 1);
  X_sq_sum_cum = fliplr(X_sq_sum_cum);
  arma::vec ans = arma::vec(4);
  ans[0] = accu(X_sq) / (N * C);
  ans[1] = (2 * Y2N) / (N * (N - 1) * pow(C, 2));  
  ans[2] = (2 * accu(X_sq_sum.cols(0, N - 2) % X_sq_sum_cum.cols(1, N - 1))) / 
    (N * (N - 1) * pow(C, 2));
  ans[3] = (sum_X_i_had_X_j * 2) / (N * (N - 1));
  return ans;
}





#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec statistics_trans_centered(arma::mat X, double N) {
  int C = X.n_cols / N;
  int R = X.n_rows;
  arma::mat X_i = arma::zeros(R, C);
  arma::mat X_j = arma::zeros(R, C);
  arma::mat X_sq_sum = arma::zeros(R, N);
  double sum_X_i_had_X_j = 0;
  double Y2N = 0;
  arma::mat X_sq = pow(X, 2);
  X_sq_sum.col(N - 1) = sum(X_sq.cols((N - 1) * C, N * C - 1), 1);
  for(int i=N-1; i>0; i--){
    X_i = X.cols((i - 1) * C, i * C - 1);
    X_sq_sum.col(i - 1) = sum(X_sq.cols((i - 1) * C, i * C - 1), 1);
    for(int j=i+1; j<N+1; j++){
      X_j = X.cols((j - 1) * C, j * C - 1);
      Y2N += accu(pow(crossprod2cpp(X_i, X_j), 2));
      sum_X_i_had_X_j += pow(accu(X_i % X_j), 2);
    }
  }
  arma::vec ans = arma::vec(4);
  ans[0] = accu(X_sq) / (N * C);
  ans[1] = (2 * Y2N) / (N * (N - 1) * pow(C, 2));  
  arma::mat X_sq_sum_cum = fliplr(X_sq_sum);
  X_sq_sum_cum = cumsum(X_sq_sum_cum, 1);
  X_sq_sum_cum = fliplr(X_sq_sum_cum);
  ans[2] = (2 * accu(X_sq_sum.cols(0, N - 2) % X_sq_sum_cum.cols(1, N - 1))) / 
    (N * (N - 1) * pow(C, 2));
  ans[3] = (sum_X_i_had_X_j * 2) / ( N *(N-1));
  return ans;
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec statistics(arma::mat X, double N) {
  int C = X.n_cols / N;
  int R = X.n_rows;
  arma::mat X_i_X_j_tcrossprod = arma::zeros(R, R);
  arma::mat X_cen_i_sq = arma::zeros(R, C);
  arma::mat X_cen_sq_colsums = arma::zeros(N, 1);
  arma::mat X_cen_sq_sum = arma::zeros(R, N);
  arma::mat X_sq_sum = arma::zeros(R, N);
  arma::mat X_j = arma::zeros(R, C);
  double sum_X_i_had_X_j = 0;
  double Y2N = 0;
  double Y4N = 0;
  double Y62N = 0;
  double Y63N = 0;
  double Y82N = 0;
  arma::mat X_i = X.cols((N - 1) * C, N * C - 1);
  arma::mat X_i_tcrossprod = tcrossprodcpp(X_i);
  arma::mat X_cen = centerdatamatrix(X, N);
  arma::mat X_cen_sq = pow(X_cen, 2);
  arma::mat X_cen_i = X_cen.cols((N - 1) * C, N * C - 1);
  arma::mat X_cen_i_tcrossprod = tcrossprodcpp(X_cen_i);
  arma::mat X_cen_i_X_i_tcrossprod = tcrossprod2cpp(X_cen_i, X_i);
  arma::mat X_sq = pow(X, 2);
  arma::mat X_sum = sumdatamatrix(X, N);
  arma::mat X_sum_gi = X_i;
  arma::mat X_sum_ni = X_sum - X_i;
  arma::mat X_i_X_sum_li_tcrossprod = tcrossprod2cpp(X_i, X_sum_ni);
  arma::mat X_tcrossprod = X_i_tcrossprod;
  arma::mat X_times_X_cen = X % X_cen;
  X_cen_sq_sum.col(N - 1) = sum(X_cen_sq.cols((N - 1) * C, N * C - 1), 1);
  X_cen_sq_colsums[N - 1] = accu(X_cen_sq_sum.col(N - 1));
  X_sq_sum.col(N - 1) = sum(X_sq.cols((N - 1) * C, N * C - 1), 1);
  double Y1N = accu(X_sq) / N;
  double Y51N = accu(X_cen_i_tcrossprod % X_i_tcrossprod);
  double Y52N = accu(pow(X_i_tcrossprod, 2));
  double Y53N = accu(X_i_tcrossprod % X_i_X_sum_li_tcrossprod);
  double Y61N = accu(pow(X_cen_i_tcrossprod,2));
  double Y641N = accu(pow(X_cen_i_X_i_tcrossprod,2));
  double Y651N = accu(X_cen_i_X_i_tcrossprod % X_cen_i_X_i_tcrossprod.t());
  double Y73N = dot(X_sq_sum.col(N - 1), sum((X_i % X_sum_ni),1));
  double Y841N = accu(pow(sum(X_times_X_cen.cols((N - 1) * C, N * C - 1), 1), 2));
  for(int i=N-1; i>0; i--){
    X_i = X.cols((i - 1) * C, i * C - 1);
    X_i_tcrossprod = tcrossprodcpp(X_i);
    X_cen_i = X_cen.cols((i - 1) * C, i * C - 1);
    X_cen_i_sq = X_cen_sq.cols((i - 1) * C, i * C - 1);
    X_cen_i_tcrossprod = tcrossprodcpp(X_cen_i);
    X_cen_i_X_i_tcrossprod = tcrossprod2cpp(X_cen_i, X_i);
    X_cen_sq_colsums[i - 1] = accu(X_cen_i_sq);
    X_sum_ni = X_sum - X_i;
    X_i_X_sum_li_tcrossprod = tcrossprod2cpp(X_i, X_sum_ni);
    X_sq_sum.col(i - 1) = sum(X_sq.cols((i - 1) * C, i * C - 1), 1);
    X_cen_sq_sum.col(i - 1) = sum(X_cen_i_sq, 1);
    Y2N += accu(X_i_tcrossprod % X_tcrossprod);
    Y4N += accu(X_i % X_sum_gi);
    Y51N += accu(X_cen_i_tcrossprod % X_i_tcrossprod);
    Y52N += accu(pow(X_i_tcrossprod, 2));
    Y53N += accu(X_i_tcrossprod % X_i_X_sum_li_tcrossprod);
    Y61N += accu(pow(X_cen_i_tcrossprod,2));
    Y641N += accu(pow(X_cen_i_X_i_tcrossprod,2));
    Y651N += accu(X_cen_i_X_i_tcrossprod % X_cen_i_X_i_tcrossprod.t());
    Y73N += dot(X_sq_sum.col(i - 1), sum((X_i % X_sum_ni),1));
    Y841N += accu(pow(sum(X_times_X_cen.cols((i - 1) * C, i * C - 1), 1), 2));
    for(int j=i+1; j<N+1; j++){
      X_j = X.cols((j - 1) * C, j * C - 1);
      X_i_X_j_tcrossprod = tcrossprod2cpp(X_i, X_j);
      Y62N += accu(pow(X_i_X_j_tcrossprod, 2));
      Y63N += accu(X_i_X_j_tcrossprod % X_i_X_j_tcrossprod.t());
      Y82N += accu(pow(sum(X_i % X_j, 1), 2));
      sum_X_i_had_X_j += pow(accu(X_cen_i % X_cen.cols((j - 1) * C, j * C - 1)), 2);
    }
    X_tcrossprod += X_i_tcrossprod;
    X_sum_gi += X_i;
  }
  Y2N = 2 * Y2N;  
  arma::mat X_sq_sum_cum = fliplr(X_sq_sum);
  X_sq_sum_cum = cumsum(X_sq_sum_cum, 1);
  X_sq_sum_cum = fliplr(X_sq_sum_cum);
  double Y3N = 2 * accu(X_sq_sum.cols(0, N - 2) % X_sq_sum_cum.cols(1, N - 1));
  Y4N = 2 * Y4N / (N * (N - 1));
  double Y5N = pow(N, 2) * Y51N - pow(N - 1, 2) * Y52N + 2 * (N - 1) * Y53N - Y2N;
  double Y64N = pow(N, 2) * Y641N + 2 * (N - 1) * Y53N - pow(N - 1, 2) * Y52N - 2 * Y62N;
  double Y65N = pow(N, 2) * Y651N + 2 * (N - 1) * Y53N - pow(N - 1, 2) * Y52N - 2 * Y63N;
  double Y6N = ((N - 1) * (pow(N, 2) - 3 * N + 3) * Y52N + 
                (2 * N - 3) * (Y2N + 2 * Y62N + 2 * Y63N) + 
                2 * (N - 3) * (Y5N + Y64N + Y65N) -
                4 * (pow(N, 2) - 3 * N + 3) * Y53N - pow(N, 3) * Y61N) / 3;
  double Y71N = accu(X_sq_sum % X_cen_sq_sum);
  double Y72N = accu(pow(X_sq_sum, 2));
  double Y7N = pow(N, 2) * Y71N - pow(N - 1, 2) * Y72N + 2 * (N - 1) * Y73N - Y3N;
  double y81N = accu(pow(X_cen_sq_sum, 2));
  double Y84N = pow(N, 2) * Y841N + 2 * (N-1) * Y73N - pow(N - 1, 2) * Y72N - 2 * Y82N;
  double Y8N = ((N - 1) * (pow(N, 2) - 3 * N + 3) * Y72N + (2 * N - 3) * (Y3N + 4 * Y82N) + 
                2 * (N - 3) * (Y7N + 2 * Y84N) - 4 * (pow(N, 2) - 3 * N + 3) * Y73N - pow(N, 3) * y81N) / 3;
  double Y9N = pow(accu(X_cen_sq_colsums) / (N - 1), 2);
  double Y11N = accu(pow(X_cen_sq_colsums, 2));
  double Y10N = (sum_X_i_had_X_j * 2 + Y11N) / pow(N - 1, 2);
  arma::vec ans = arma::vec(4);
  ans[0] = (Y1N - Y4N) / C;
  ans[1] = (Y2N / (N * (N - 1)) - 2 * Y5N / (N * (N - 1) * (N - 2)) + Y6N / (N * (N - 1) * (N - 2) * (N-3))) / pow(C, 2);
  ans[2] = (Y3N / (N * (N - 1)) - 2 * Y7N / (N * (N - 1) * (N - 2)) + Y8N / (N * (N - 1) * (N - 2) * (N-3))) / pow(C, 2);
  ans[3] = (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 1) * (N - 2) * Y10N + Y9N - N * Y11N / (N - 1));
  return ans;
}



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec statistics_trans(arma::mat X, double N) {
  int C = X.n_cols / N;
  int R = X.n_rows;
  arma::mat X_i_X_j_crossprod = arma::zeros(C, C);
  arma::mat X_cen_i_sq = arma::zeros(R, C);
  arma::mat X_cen_sq_colsums = arma::zeros(N, 1);
  arma::mat X_cen_sq_sum = arma::zeros(R, N);
  arma::mat X_sq_sum = arma::zeros(R, N);
  arma::mat X_j = arma::zeros(R, C);
  double sum_X_i_had_X_j = 0;
  double Y2N = 0;
  double Y4N = 0;
  double Y62N = 0;
  double Y63N = 0;
  double Y82N = 0;
  arma::mat X_i = X.cols((N - 1) * C, N * C - 1);
  arma::mat X_i_crossprod = crossprodcpp(X_i);
  arma::mat X_cen = centerdatamatrix(X, N);
  arma::mat X_cen_sq = pow(X_cen, 2);
  arma::mat X_cen_i = X_cen.cols((N - 1) * C, N * C - 1);
  arma::mat X_cen_i_crossprod = crossprodcpp(X_cen_i);
  arma::mat X_cen_i_X_i_crossprod = crossprod2cpp(X_cen_i, X_i);
  arma::mat X_sq = pow(X, 2);
  arma::mat X_sum = sumdatamatrix(X, N);
  arma::mat X_sum_gi = X_i;
  arma::mat X_sum_ni = X_sum - X_i;
  arma::mat X_sum_ni_X_i_crossprod = crossprod2cpp(X_sum_ni, X_i);
  arma::mat X_crossprod = X_i_crossprod;
  arma::mat X_times_X_cen = X % X_cen;
  X_cen_sq_sum.col(N - 1) = sum(X_cen_sq.cols((N - 1) * C, N * C - 1), 1);
  X_cen_sq_colsums[N - 1] = accu(X_cen_sq_sum.col(N - 1));
  X_sq_sum.col(N - 1) = sum(X_sq.cols((N - 1) * C, N * C - 1), 1);
  double Y1N = accu(X_sq) / N;
  double Y51N = accu(pow(X_cen_i_X_i_crossprod, 2));
  double Y52N = accu(pow(X_i_crossprod, 2));
  double Y53N = accu(X_i_crossprod % X_sum_ni_X_i_crossprod);
  double Y61N = accu(pow(X_cen_i_crossprod,2));
  double Y641N = accu(X_i_crossprod % X_cen_i_crossprod);
  double Y651N = accu(X_cen_i_X_i_crossprod % X_cen_i_X_i_crossprod.t());
  double Y73N = dot(X_sq_sum.col(N - 1), sum((X_i % X_sum_ni),1));
  double Y841N = accu(pow(sum(X_times_X_cen.cols((N - 1) * C, N * C - 1), 1), 2));
  for(int i=N-1; i>0; i--){
    X_i = X.cols((i - 1) * C, i * C - 1);
    X_i_crossprod = crossprodcpp(X_i);
    X_cen_i = X_cen.cols((i - 1) * C, i * C - 1);
    X_cen_i_sq = X_cen_sq.cols((i - 1) * C, i * C - 1);
    X_cen_i_crossprod = crossprodcpp(X_cen_i);
    X_cen_i_X_i_crossprod = crossprod2cpp(X_cen_i, X_i);
    X_cen_sq_colsums[i - 1] = accu(X_cen_i_sq);
    X_sum_ni = X_sum - X_i;
    X_sum_ni_X_i_crossprod = crossprod2cpp(X_sum_ni, X_i);
    X_sq_sum.col(i - 1) = sum(X_sq.cols((i - 1) * C, i * C - 1), 1);
    X_cen_sq_sum.col(i - 1) = sum(X_cen_i_sq, 1);
    Y4N += accu(X_i % X_sum_gi);
    Y51N += accu(pow(X_cen_i_X_i_crossprod, 2));
    Y52N += accu(pow(X_i_crossprod, 2));
    Y53N += accu(X_i_crossprod % X_sum_ni_X_i_crossprod);
    Y61N += accu(pow(X_cen_i_crossprod,2));
    Y641N += accu(X_i_crossprod % X_cen_i_crossprod);
    Y651N += accu(X_cen_i_X_i_crossprod % X_cen_i_X_i_crossprod.t());
    Y62N += accu(X_i_crossprod % X_crossprod);
    Y73N += dot(X_sq_sum.col(i - 1), sum((X_i % X_sum_ni),1));
    Y841N += accu(pow(sum(X_times_X_cen.cols((i - 1) * C, i * C - 1), 1), 2));
    for(int j=i+1; j<N+1; j++){
      X_j = X.cols((j - 1) * C, j * C - 1);
      X_i_X_j_crossprod = crossprod2cpp(X_i, X_j);
      Y2N += accu(pow(X_i_X_j_crossprod, 2));
      Y63N += accu(X_i_X_j_crossprod % X_i_X_j_crossprod.t());
      Y82N += accu(pow(sum(X_i % X_j, 1), 2));
      sum_X_i_had_X_j += pow(accu(X_cen_i % X_cen.cols((j - 1) * C, j * C - 1)), 2);
    }
    X_crossprod += X_i_crossprod;
    X_sum_gi += X_i;
  }
  Y2N = 2 * Y2N;  
  arma::mat X_sq_sum_cum = fliplr(X_sq_sum);
  X_sq_sum_cum = cumsum(X_sq_sum_cum, 1);
  X_sq_sum_cum = fliplr(X_sq_sum_cum);
  double Y3N = 2 * accu(X_sq_sum.cols(0, N - 2) % X_sq_sum_cum.cols(1, N - 1));
  Y4N = 2 * Y4N / (N * (N - 1));
  double Y5N = pow(N, 2) * Y51N - pow(N - 1, 2) * Y52N + 2 * (N - 1) * Y53N - Y2N;
  double Y64N = pow(N, 2) * Y641N + 2 * (N - 1) * Y53N - pow(N - 1, 2) * Y52N - 2 * Y62N;
  double Y65N = pow(N, 2) * Y651N + 2 * (N - 1) * Y53N - pow(N - 1, 2) * Y52N - 2 * Y63N;
  double Y6N = ((N - 1) * (pow(N, 2) - 3 * N + 3) * Y52N + 
                (2 * N - 3) * (Y2N + 2 * Y62N + 2 * Y63N) + 
                2 * (N - 3) * (Y5N + Y64N + Y65N) -
                4 * (pow(N, 2) - 3 * N + 3) * Y53N - pow(N, 3) * Y61N) / 3;
  double Y71N = accu(X_sq_sum % X_cen_sq_sum);
  double Y72N = accu(pow(X_sq_sum, 2));
  double Y7N = pow(N, 2) * Y71N - pow(N - 1, 2) * Y72N + 2 * (N - 1) * Y73N - Y3N;
  double y81N = accu(pow(X_cen_sq_sum, 2));
  double Y84N = pow(N, 2) * Y841N + 2 * (N-1) * Y73N - pow(N - 1, 2) * Y72N - 2 * Y82N;
  double Y8N = ((N - 1) * (pow(N, 2) - 3 * N + 3) * Y72N + (2 * N - 3) * (Y3N + 4 * Y82N) + 
                2 * (N - 3) * (Y7N + 2 * Y84N) - 4 * (pow(N, 2) - 3 * N + 3) * Y73N - pow(N, 3) * y81N) / 3;
  double Y9N = pow(accu(X_cen_sq_colsums) / (N - 1), 2);
  double Y11N = accu(pow(X_cen_sq_colsums, 2));
  double Y10N = (sum_X_i_had_X_j * 2 + Y11N) / pow(N - 1, 2);
  arma::vec ans = arma::vec(4);
  ans[0] = (Y1N - Y4N) / C;
  ans[1] = (Y2N / (N * (N - 1)) - 2 * Y5N / (N * (N - 1) * (N - 2)) + Y6N / (N * (N - 1) * (N - 2) * (N-3))) / pow(C, 2);
  ans[2] = (Y3N / (N * (N - 1)) - 2 * Y7N / (N * (N - 1) * (N - 2)) + Y8N / (N * (N - 1) * (N - 2) * (N-3))) / pow(C, 2);
  ans[3] = (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 1) * (N - 2) * Y10N + Y9N - N * Y11N / (N - 1));
  return ans;
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec covmathat_statistics_centered(arma::mat X, double N) {
  int C = X.n_cols / N;
  int R = X.n_rows;
  arma::mat X_i_tcrossprod = arma::zeros(R, R);
  double sum_X_i_had_X_j = 0;
  double Y2N = 0;
  arma::mat X_i = X.cols((N-1) * C, N * C - 1);
  arma::mat X_tcrossprod = tcrossprodcpp(X_i);
  arma::mat X_sq = pow(X, 2);
  for(int i=N-1; i>0; i--){
    X_i = X.cols((i - 1) * C, i * C - 1);
    X_i_tcrossprod = tcrossprodcpp(X_i);
    Y2N += accu(X_i_tcrossprod % X_tcrossprod);
    X_tcrossprod += X_i_tcrossprod;
    for(int j=i+1; j<N+1; j++){
      sum_X_i_had_X_j += pow(accu(X_i % X.cols((j - 1) * C, j * C - 1)), 2);
    }
  }
  arma::vec ans = arma::vec(3);
  ans[0] = accu(X_sq) / (N * C);
  ans[1] = (2 * Y2N) / (N * (N - 1) * pow(C, 2));  
  ans[2] = (sum_X_i_had_X_j * 2) / (N * (N - 1));
  return ans;
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec covmathat_statistics_trans_centered(arma::mat X, double N) {
  int C = X.n_cols / N;
  int R = X.n_rows;
  arma::mat X_i = arma::zeros(R, C);
  arma::mat X_j = arma::zeros(R, C);
  double sum_X_i_had_X_j = 0;
  double Y2N = 0;
  arma::mat X_sq = pow(X, 2);
  for(int i=N-1; i>0; i--){
    X_i = X.cols((i - 1) * C, i * C - 1);
    for(int j=i+1; j<N+1; j++){
      X_j = X.cols((j - 1) * C, j * C - 1);
      Y2N += accu(pow(crossprod2cpp(X_i, X_j), 2));
      sum_X_i_had_X_j += pow(accu(X_i % X_j), 2);
    }
  }
  arma::vec ans = arma::vec(3);
  ans[0] = accu(X_sq) / (N * C);
  ans[1] = (2 * Y2N) / (N * (N - 1) * pow(C, 2));  
  ans[2] = (sum_X_i_had_X_j * 2) / ( N *(N-1));
  return ans;
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec covmathat_statistics(arma::mat X, double N) {
  int C = X.n_cols / N;
  int R = X.n_rows;
  arma::mat X_i_X_j_tcrossprod = arma::zeros(R, R);
  arma::mat X_cen_i_sq = arma::zeros(R, C);
  arma::mat X_cen_sq_colsums = arma::zeros(N, 1);
  arma::mat X_cen_sq_sum = arma::zeros(R, N);
  arma::mat X_sq_sum = arma::zeros(R, N);
  arma::mat X_j = arma::zeros(R, C);
  double sum_X_i_had_X_j = 0;
  double Y2N = 0;
  double Y4N = 0;
  double Y62N = 0;
  double Y63N = 0;
  arma::mat X_i = X.cols((N - 1) * C, N * C - 1);
  arma::mat X_i_tcrossprod = tcrossprodcpp(X_i);
  arma::mat X_cen = centerdatamatrix(X, N);
  arma::mat X_cen_sq = pow(X_cen, 2);
  arma::mat X_cen_i = X_cen.cols((N - 1) * C, N * C - 1);
  arma::mat X_cen_i_tcrossprod = tcrossprodcpp(X_cen_i);
  arma::mat X_cen_i_X_i_tcrossprod = tcrossprod2cpp(X_cen_i, X_i);
  arma::mat X_sq = pow(X, 2);
  arma::mat X_sum = sumdatamatrix(X, N);
  arma::mat X_sum_gi = X_i;
  arma::mat X_sum_ni = X_sum - X_i;
  arma::mat X_i_X_sum_li_tcrossprod = tcrossprod2cpp(X_i, X_sum_ni);
  arma::mat X_tcrossprod = X_i_tcrossprod;
  arma::mat X_times_X_cen = X % X_cen;
  X_cen_sq_sum.col(N - 1) = sum(X_cen_sq.cols((N - 1) * C, N * C - 1), 1);
  X_cen_sq_colsums[N - 1] = accu(X_cen_sq_sum.col(N - 1));
  X_sq_sum.col(N - 1) = sum(X_sq.cols((N - 1) * C, N * C - 1), 1);
  double Y1N = accu(X_sq) / N;
  double Y51N = accu(X_cen_i_tcrossprod % X_i_tcrossprod);
  double Y52N = accu(pow(X_i_tcrossprod, 2));
  double Y53N = accu(X_i_tcrossprod % X_i_X_sum_li_tcrossprod);
  double Y61N = accu(pow(X_cen_i_tcrossprod,2));
  double Y641N = accu(pow(X_cen_i_X_i_tcrossprod,2));
  double Y651N = accu(X_cen_i_X_i_tcrossprod % X_cen_i_X_i_tcrossprod.t());
  for(int i=N-1; i>0; i--){
    X_i = X.cols((i - 1) * C, i * C - 1);
    X_i_tcrossprod = tcrossprodcpp(X_i);
    X_cen_i = X_cen.cols((i - 1) * C, i * C - 1);
    X_cen_i_sq = X_cen_sq.cols((i - 1) * C, i * C - 1);
    X_cen_i_tcrossprod = tcrossprodcpp(X_cen_i);
    X_cen_i_X_i_tcrossprod = tcrossprod2cpp(X_cen_i, X_i);
    X_cen_sq_colsums[i - 1] = accu(X_cen_i_sq);
    X_sum_ni = X_sum - X_i;
    X_i_X_sum_li_tcrossprod = tcrossprod2cpp(X_i, X_sum_ni);
    X_sq_sum.col(i - 1) = sum(X_sq.cols((i - 1) * C, i * C - 1), 1);
    X_cen_sq_sum.col(i - 1) = sum(X_cen_i_sq, 1);
    Y2N += accu(X_i_tcrossprod % X_tcrossprod);
    Y4N += accu(X_i % X_sum_gi);
    Y51N += accu(X_cen_i_tcrossprod % X_i_tcrossprod);
    Y52N += accu(pow(X_i_tcrossprod, 2));
    Y53N += accu(X_i_tcrossprod % X_i_X_sum_li_tcrossprod);
    Y61N += accu(pow(X_cen_i_tcrossprod,2));
    Y641N += accu(pow(X_cen_i_X_i_tcrossprod,2));
    Y651N += accu(X_cen_i_X_i_tcrossprod % X_cen_i_X_i_tcrossprod.t());
    for(int j=i+1; j<N+1; j++){
      X_j = X.cols((j - 1) * C, j * C - 1);
      X_i_X_j_tcrossprod = tcrossprod2cpp(X_i, X_j);
      Y62N += accu(pow(X_i_X_j_tcrossprod, 2));
      Y63N += accu(X_i_X_j_tcrossprod % X_i_X_j_tcrossprod.t());
      sum_X_i_had_X_j += pow(accu(X_cen_i % X_cen.cols((j - 1) * C, j * C - 1)), 2);
    }
    X_tcrossprod += X_i_tcrossprod;
    X_sum_gi += X_i;
  }
  Y2N = 2 * Y2N;  
  arma::mat X_sq_sum_cum = fliplr(X_sq_sum);
  X_sq_sum_cum = cumsum(X_sq_sum_cum, 1);
  X_sq_sum_cum = fliplr(X_sq_sum_cum);
  Y4N = 2 * Y4N / (N * (N - 1));
  double Y5N = pow(N, 2) * Y51N - pow(N - 1, 2) * Y52N + 2 * (N - 1) * Y53N - Y2N;
  double Y64N = pow(N, 2) * Y641N + 2 * (N - 1) * Y53N - pow(N - 1, 2) * Y52N - 2 * Y62N;
  double Y65N = pow(N, 2) * Y651N + 2 * (N - 1) * Y53N - pow(N - 1, 2) * Y52N - 2 * Y63N;
  double Y6N = ((N - 1) * (pow(N, 2) - 3 * N + 3) * Y52N + 
                (2 * N - 3) * (Y2N + 2 * Y62N + 2 * Y63N) + 
                2 * (N - 3) * (Y5N + Y64N + Y65N) -
                4 * (pow(N, 2) - 3 * N + 3) * Y53N - pow(N, 3) * Y61N) / 3;
  double Y9N = pow(accu(X_cen_sq_colsums) / (N - 1), 2);
  double Y11N = accu(pow(X_cen_sq_colsums, 2));
  double Y10N = (sum_X_i_had_X_j * 2 + Y11N) / pow(N - 1, 2);
  arma::vec ans = arma::vec(3);
  ans[0] = (Y1N - Y4N) / C;
  ans[1] = (Y2N / (N * (N - 1)) - 2 * Y5N / (N * (N - 1) * (N - 2)) + Y6N / (N * (N - 1) * (N - 2) * (N-3))) / pow(C, 2);
  ans[2] = (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 1) * (N - 2) * Y10N + Y9N - N * Y11N / (N - 1));
  return ans;
}



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec covmathat_statistics_trans(arma::mat X, double N) {
  int C = X.n_cols / N;
  int R = X.n_rows;
  arma::mat X_i_X_j_crossprod = arma::zeros(C, C);
  arma::mat X_cen_i_sq = arma::zeros(R, C);
  arma::mat X_cen_sq_colsums = arma::zeros(N, 1);
  arma::mat X_cen_sq_sum = arma::zeros(R, N);
  arma::mat X_sq_sum = arma::zeros(R, N);
  arma::mat X_j = arma::zeros(R, C);
  double sum_X_i_had_X_j = 0;
  double Y2N = 0;
  double Y4N = 0;
  double Y62N = 0;
  double Y63N = 0;
  arma::mat X_i = X.cols((N - 1) * C, N * C - 1);
  arma::mat X_i_crossprod = crossprodcpp(X_i);
  arma::mat X_cen = centerdatamatrix(X, N);
  arma::mat X_cen_sq = pow(X_cen, 2);
  arma::mat X_cen_i = X_cen.cols((N - 1) * C, N * C - 1);
  arma::mat X_cen_i_crossprod = crossprodcpp(X_cen_i);
  arma::mat X_cen_i_X_i_crossprod = crossprod2cpp(X_cen_i, X_i);
  arma::mat X_sq = pow(X, 2);
  arma::mat X_sum = sumdatamatrix(X, N);
  arma::mat X_sum_gi = X_i;
  arma::mat X_sum_ni = X_sum - X_i;
  arma::mat X_sum_ni_X_i_crossprod = crossprod2cpp(X_sum_ni, X_i);
  arma::mat X_crossprod = X_i_crossprod;
  arma::mat X_times_X_cen = X % X_cen;
  X_cen_sq_sum.col(N - 1) = sum(X_cen_sq.cols((N - 1) * C, N * C - 1), 1);
  X_cen_sq_colsums[N - 1] = accu(X_cen_sq_sum.col(N - 1));
  X_sq_sum.col(N - 1) = sum(X_sq.cols((N - 1) * C, N * C - 1), 1);
  double Y1N = accu(X_sq) / N;
  double Y51N = accu(pow(X_cen_i_X_i_crossprod, 2));
  double Y52N = accu(pow(X_i_crossprod, 2));
  double Y53N = accu(X_i_crossprod % X_sum_ni_X_i_crossprod);
  double Y61N = accu(pow(X_cen_i_crossprod,2));
  double Y641N = accu(X_i_crossprod % X_cen_i_crossprod);
  double Y651N = accu(X_cen_i_X_i_crossprod % X_cen_i_X_i_crossprod.t());
  for(int i=N-1; i>0; i--){
    X_i = X.cols((i - 1) * C, i * C - 1);
    X_i_crossprod = crossprodcpp(X_i);
    X_cen_i = X_cen.cols((i - 1) * C, i * C - 1);
    X_cen_i_sq = X_cen_sq.cols((i - 1) * C, i * C - 1);
    X_cen_i_crossprod = crossprodcpp(X_cen_i);
    X_cen_i_X_i_crossprod = crossprod2cpp(X_cen_i, X_i);
    X_cen_sq_colsums[i - 1] = accu(X_cen_i_sq);
    X_sum_ni = X_sum - X_i;
    X_sum_ni_X_i_crossprod = crossprod2cpp(X_sum_ni, X_i);
    X_sq_sum.col(i - 1) = sum(X_sq.cols((i - 1) * C, i * C - 1), 1);
    X_cen_sq_sum.col(i - 1) = sum(X_cen_i_sq, 1);
    Y4N += accu(X_i % X_sum_gi);
    Y51N += accu(pow(X_cen_i_X_i_crossprod, 2));
    Y52N += accu(pow(X_i_crossprod, 2));
    Y53N += accu(X_i_crossprod % X_sum_ni_X_i_crossprod);
    Y61N += accu(pow(X_cen_i_crossprod,2));
    Y641N += accu(X_i_crossprod % X_cen_i_crossprod);
    Y651N += accu(X_cen_i_X_i_crossprod % X_cen_i_X_i_crossprod.t());
    Y62N += accu(X_i_crossprod % X_crossprod);
    for(int j=i+1; j<N+1; j++){
      X_j = X.cols((j - 1) * C, j * C - 1);
      X_i_X_j_crossprod = crossprod2cpp(X_i, X_j);
      Y2N += accu(pow(X_i_X_j_crossprod, 2));
      Y63N += accu(X_i_X_j_crossprod % X_i_X_j_crossprod.t());
      sum_X_i_had_X_j += pow(accu(X_cen_i % X_cen.cols((j - 1) * C, j * C - 1)), 2);
    }
    X_crossprod += X_i_crossprod;
    X_sum_gi += X_i;
  }
  Y2N = 2 * Y2N;  
  arma::mat X_sq_sum_cum = fliplr(X_sq_sum);
  X_sq_sum_cum = cumsum(X_sq_sum_cum, 1);
  X_sq_sum_cum = fliplr(X_sq_sum_cum);
  Y4N = 2 * Y4N / (N * (N - 1));
  double Y5N = pow(N, 2) * Y51N - pow(N - 1, 2) * Y52N + 2 * (N - 1) * Y53N - Y2N;
  double Y64N = pow(N, 2) * Y641N + 2 * (N - 1) * Y53N - pow(N - 1, 2) * Y52N - 2 * Y62N;
  double Y65N = pow(N, 2) * Y651N + 2 * (N - 1) * Y53N - pow(N - 1, 2) * Y52N - 2 * Y63N;
  double Y6N = ((N - 1) * (pow(N, 2) - 3 * N + 3) * Y52N + 
                (2 * N - 3) * (Y2N + 2 * Y62N + 2 * Y63N) + 
                2 * (N - 3) * (Y5N + Y64N + Y65N) -
                4 * (pow(N, 2) - 3 * N + 3) * Y53N - pow(N, 3) * Y61N) / 3;
  double Y9N = pow(accu(X_cen_sq_colsums) / (N - 1), 2);
  double Y11N = accu(pow(X_cen_sq_colsums, 2));
  double Y10N = (sum_X_i_had_X_j * 2 + Y11N) / pow(N - 1, 2);
  arma::vec ans = arma::vec(3);
  ans[0] = (Y1N - Y4N) / C;
  ans[1] = (Y2N / (N * (N - 1)) - 2 * Y5N / (N * (N - 1) * (N - 2)) + Y6N / (N * (N - 1) * (N - 2) * (N-3))) / pow(C, 2);
  ans[2] = (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 1) * (N - 2) * Y10N + Y9N - N * Y11N / (N - 1));
  return ans;
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec meanmatts_statistics(arma::mat X, double N) {
  int C = X.n_cols / N;
  int R = X.n_rows;
  arma::mat X_cen_sq_colsums = arma::zeros(N, 1);
  arma::mat X_cen_sq_sum = arma::zeros(R, N);
  double sum_X_i_had_X_j = 0;
  double Y4N = 0;
  arma::mat X_i = X.cols((N - 1) * C, N * C - 1);
  arma::mat X_cen = centerdatamatrix(X, N);
  arma::mat X_cen_sq = pow(X_cen, 2);
  arma::mat X_sum_gi = X_i;
  arma::mat X_cen_i = X_cen.cols((N - 1) * C, N * C - 1);
  X_cen_sq_colsums(N-1) = accu(pow(X_cen_i, 2));
  for(int i=N-1; i>0; i--){
    X_i = X.cols((i - 1) * C, i * C - 1);
    X_cen_i = X_cen.cols((i - 1) * C, i * C - 1);
    X_cen_sq_colsums(i-1) = accu(pow(X_cen_i, 2));
    Y4N += accu(X_i % X_sum_gi);
    for(int j=i+1; j<N+1; j++){
      sum_X_i_had_X_j += pow(accu(X_cen_i % X_cen.cols((j - 1) * C, j * C - 1)), 2);
    }
    X_sum_gi += X_i;
  }
  Y4N = 2 * Y4N / (N * (N - 1));
  double Y9N = pow(accu(X_cen_sq_colsums) / (N - 1), 2);
  double Y11N = accu(pow(X_cen_sq_colsums, 2));
  double Y10N = (sum_X_i_had_X_j * 2 + Y11N) / pow(N - 1, 2);
  arma::vec ans = arma::vec(2);
  ans[0] = Y4N;
  ans[1] = (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 1) * (N - 2) * Y10N + Y9N - N * Y11N / (N - 1));
  return ans;
}