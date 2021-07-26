#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>
#include <complex.h>
#include <fftw3.h>

//using namespace Rcpp;
//using namespace std;
//using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
int mod(int a, int n)
{
    return a - floor(a/n)*n;
}  

//*****************************************************************************//

//[[Rcpp::export]]
void l_vec_compute_arma(int k, arma::vec& l_vec, arma::vec& cn_vec, int m)
{
  int i, aa, bb;
  
  for(i=0; i<m-1; i++)
  {
    aa=mod(k, cn_vec(i));
    bb=(k-aa)/cn_vec(i);
    l_vec(i)=bb;
    k=aa;
  }

  return;
}


/*
arma::vec test(int mm, arma::vec nn_vec)
{  
  int i;
  
  int nn_vec_a[mm-1] = {0};
 
  for(i=0; i<mm-1; i++)
  {

    Rcout<<"i= "<< i<<"  value "<< nn_vec_a[i] <<endl;
    
    nn_vec_a[i]=nn_vec(i);
    
    Rcout<<"i= "<< i<<"  value "<< nn_vec_a[i] <<endl;
  
  }
  
  return nn_vec;
    
}  
*/


//[[Rcpp::export]]  
arma::vec pmn_mdfft_arma(int nnt, arma::mat pp, arma::vec nn_vec, arma::vec l_vec, arma::vec cn_vec)
{  
  arma::vec res(nnt, arma::fill::zeros);

  int mm=pp.n_cols;  
  int nn=pp.n_rows;
  
  //int nn_vec_a[mm-1] = {0};
  int* nn_vec_a = new int[mm-1];
  
  fftw_complex *in, *out;
  int i, j, k;
  int n, m, nt;
  fftw_plan p;
  double tmp, con, pij, pim, ww;
  arma::cx_double ctmp, ctmp1, ctmp2, qval, a1, a2;
  
  nt=nnt;
  n=nn;
  m=mm;
  
  
  for(i=0; i<mm-1; i++)
  {
    nn_vec_a[i]=nn_vec(i);
  }
  
  //Rprintf("nt %u, n %u, m %u \n", nt, n, m);
 
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt); 
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);  
  
  ww=2*PI/(n+1);

  //Rprintf("ww %lf \n", ww);
  
  for(k=0; k<nt; k++)
  {
    qval=0.0 + 0.0 * I;

    l_vec_compute_arma(k, l_vec, cn_vec, m);
    
    //for(ii=0; ii<m-1; ii++)
    //{
    //  Rprintf("l_vec %u \n", l_vec[ii]);
    //}
      
    for(i=0; i<n; i++)
    {
      ctmp=0.0 + 0.0 * I; 
      
      for(j=0; j<m-1; j++)
      {
        pij=pp(i,j);

        //printf("pij: %lf, l_vec[j], %u, ww, %lf, \n", pij, l_vec[j], ww);
        
        ctmp1=0.0+l_vec(j)*ww*I;
        
        //printf("ctmp1: %lf +%lf*i\n", creal(ctmp1), cimag(ctmp1));

        a1=pij+0.0*I;
        a2=exp(ctmp1);
        //printf("a1: %lf +%lf*i\n", creal(a1), cimag(a1));
        //printf("a2: %lf +%lf*i\n", creal(a2), cimag(a2));
                
        ctmp2=a1*a2;

        //printf("ctmp2: %lf +%lf*i\n\n", creal(ctmp2), cimag(ctmp2));
        
        ctmp+=ctmp2;
      }
      
      pim=pp(i, m-1);
      ctmp+=pim;
      
      ctmp=log(ctmp);
      qval+=ctmp;
    }
    
    qval=exp(qval);
    
    in[k][0]= real(qval);
    in[k][1]= imag(qval);
    
    //printf("qval: %lf +%lf*i\n", creal(qval), cimag(qval));
    //printf("in[k]: %lf +%lf*i\n", creal(in[k]), cimag(in[k]));

  
  }
  
  p=fftw_plan_dft(m-1, nn_vec_a, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  fftw_execute(p);
  
  con=pow(n+1, m-1);
  
  for(k=0; k<nt; k++)
  {
    tmp=out[k][0];
    res[k]=tmp/con;
    //printf("out[k]: %lf +%lf*i\n", creal(out[k]), cimag(out[k]));
  }
  
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);  
  
  delete[] nn_vec_a;
  
  return res;
}











