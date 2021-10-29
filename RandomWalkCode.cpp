
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
sp_mat prob_walk_net(sp_mat data, int steps) {
  steps = steps-1;
  
  mat data_2 = mat(data);
  sp_mat inv_deg = sp_mat(data);
  inv_deg.zeros();
  sp_mat deg = sum(data,1);
    sp_mat::iterator it = deg.begin();
  sp_mat::iterator it_end = deg.end();
 
  for(; it!=it_end; ++it){

    inv_deg(it.row(), it.row()) = 1 / (*it);
  }
  



  sp_mat n_sums = inv_deg * data;
  sp_mat walk_sums = n_sums;
  sp_mat prev_sums = n_sums;
  sp_mat walk_dist = n_sums;
  for(int i = 1; i <= steps; i++){
    walk_dist = walk_dist*n_sums;
    walk_sums = walk_dist%(1-walk_sums)+walk_sums;
  }
  return walk_sums;
}

// [[Rcpp::export]]
sp_mat strength_preserve_rewire(sp_mat data, mat out_in, double mean_k, double mean_s, mat out_in_deg){
  
  sp_mat::iterator it = data.begin();
  sp_mat::iterator it_end = data.end();
  
  for(; it!=it_end; ++it){
   // Rcout << it.row() <<" " <<it.col() << std::endl;
    //Rcout << (mean_k/mean_s) << " " <<out_in(it.row(),0) << " " <<out_in(it.col(), 1) << " " << out_in_deg(it.col(),0) <<  " " <<out_in_deg(it.row(),1) << std::endl;
    
     *it = (mean_k/mean_s)*out_in(it.row(),0)*out_in(it.col(), 1)/(out_in_deg(it.col(),0)*out_in_deg(it.row(),1));
  }
  
  return data;
  
  
}