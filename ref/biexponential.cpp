#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "zeroin.h"

struct biexponential_info {
  double a, b, c,d, f, w, y;
};

double biexponential_fn(double x, void*info) {
  struct biexponential_info *p = (struct biexponential_info *)info;
  double B = p->a*exp(p->b*(x-p->w))-p->c*exp(-p->d*(x-p->w))+p->f -p-> y;
  return B;
}

struct sfun_info{
  double m, w, p, t, a, r;
};

/*
Logicle function calculation

double strans_fn(double y,void*info){
	struct sfun_info  *k = (struct sfun_info *)info;
	double tmp ;
	tmp = (y < (k->w +k->a)) ? (k->w + k->a - y) : (y - k->w - k->a) ;
	tmp =  k->t*pow(10,-1*(k->m - k->w - k->a ))*(pow(10,tmp)-k->p*k->p*pow(10,-tmp/k->p)+ k->p*k->p-1);
	return((y < (k->w +k->a)) ? -1*tmp -k->r  :tmp-k->r);
}
*/
/*10/15/09 ngopalak: Updated biexponential transform so that it now calculates
the solution of a generic biexponential function
S(x,a,b,c,d,f) = ae^(bx) -ce^(-dx) +f  instead of the logicle function
it was calculating earlier
*/

//[[Rcpp::export]]
std::vector<double> biexponential_transform(std::vector<double> input,
                                            double A, double B, double C,
                                            double D, double F, double W,
                                            double tol, int maxIt) {
  struct biexponential_info params;
  params.a = A;
  params.b = B;
  params.c = C;
  params.d = D;
  params.f = F;
  params.w = W;
  void* pParams = (void*)&params;
  int fail = 0;
  for(int i=0; i < input.size(); i++) {
    params.y = input.at(i);
    double step = 0.5;
    for(int j=0; j<maxIt; j++, step*=1.5) {
      double bi1 = biexponential_fn(step, pParams);
      double bi2 = biexponential_fn(-step, pParams);
      if (bi1 * bi2 > 0) {
        break;
      }
    }
    double tol_ = tol;
    int maxIt_ = maxIt;
    input.at(i) = R_zeroin(-step, step, biexponential_fn, pParams, &tol_, &maxIt_);
    if(MaxIt==-1){
      fail=fail+1;
    }
  }
  if(fail>0)
    Rcpp::warning("%d values have not converged.", fail);

  return input;
}
