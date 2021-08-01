#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>


#define REAL 0
#define IMAG 1
#define PI2 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651e+00

//using namespace Rcpp;
//using namespace std;
//using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]] 
Rcpp::IntegerVector rmultinom_1(unsigned int &size,Rcpp:: NumericVector &probs, unsigned int &N) {
    Rcpp::IntegerVector outcome(N);
    rmultinom(size, probs.begin(), N, outcome.begin());
    return outcome;
}

//[[Rcpp::export]] 
Rcpp::IntegerMatrix rmultinom_rcpp(unsigned int &n, unsigned int &size, Rcpp::NumericVector &probs) {
    unsigned int N = probs.length();
    Rcpp::IntegerMatrix sim(N, n);
    for (unsigned int i = 0; i < n; i++) {
        sim(Rcpp::_,i) = rmultinom_1(size, probs, N);
    }
    return sim;
}

//[[Rcpp::export]] 
arma::vec rpmd_arma(arma::mat pp)
{
  int mm=pp.n_cols; 
  int nn=pp.n_rows;
  int i;
  
  unsigned int n=1;
  unsigned int size=1;
  
  arma::mat tmp(nn, mm, arma::fill::zeros);
  
  for(i=0;i<nn;i++){
    Rcpp::NumericVector prob=Rcpp::wrap(pp.row(i));
    
    Rcpp::IntegerMatrix res=rmultinom_rcpp(n, size, prob);
    /*tmp.row(i)=arma::conv_to<arma::rowvec>::from(res);*/

    arma::mat xx = Rcpp::as<arma::vec>(res);
    tmp.row(i)=xx.as_row();
  
  }
  
  arma::vec finalres(mm, arma::fill::zeros);
  finalres=sum(tmp).as_col();
  
  return finalres;
}
//[[Rcpp::export]] 

double pm_simulation_arma(arma::mat pp, arma::vec x_vec, int t)
{
	/*arma::vec res(nnt, arma::fill::zeros);*/
	double res=0;
	int mm=pp.n_cols;  

  int k, u;
  double count;
  

  arma::mat sim(t,mm,arma::fill::zeros);
  for(k=0;k<t;k++){
    //arma::mat tmp(nn, mm, arma::fill::zeros);
    //rpmd(pp).as_row().print();
    sim.row(k)=rpmd_arma(pp).as_row();
    //sim.row(k).print();
  }



	count=0;
	for(u=0;u<t;u++){
		if(all(sim.row(u)==x_vec.as_row())){
		  count++;
		  }
	

	}

	res=count/t;


	return res;
}