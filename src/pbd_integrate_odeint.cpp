//' @useDynLib PBD


#define STRICT_R_HEADERS
#include "config.h"
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "odeint_helper.h"


using namespace Rcpp;


class ode_rhs
{
public:
  ode_rhs(NumericVector parsvec)
  {
    b_G  = parsvec[0];
    b_I  = parsvec[1];
    mu_G = parsvec[2];
    mu_I = parsvec[3];
    la   = parsvec[4];
  }

  void operator()(const std::vector<double>& xx, std::vector<double>& dx, double /* t */)
  {
    dx.front() = dx.back() = 0.0;
    const size_t lx = xx.size() - 1;
    dx[1] = (mu_I + la) * xx[2]
          - (mu_G + b_G) * xx[1];
	  for (size_t i = 1; i < lx/2; ++i) {
      const size_t i0 = i - 1;
      const size_t i1 = i + 1;
      dx[1 + i] = (b_G + b_I * i0) * xx[1 + i0]
                + (mu_I + la) * i1 * xx[1 + i1] * (i1 < lx/2)
                - (mu_G + b_G + (mu_I + b_I + la) * i) * xx[1 + i];
      dx[1] = dx[1] + la * i * xx[1 + i] + la * i * xx[lx/2 + i];
      dx[lx/2 + i] = b_I * i0 * xx[lx/2 + i0] * (i0 > 0)
                   + (mu_I + la) * i1 * xx[lx/2 + i1]
                   - (mu_I + b_I + la) * i * xx[lx/2 + i]
                   + mu_G * xx[1 + i];
    }
  }

private:
  double b_G;
  double b_I;
  double mu_G;
  double mu_I;
  double la;
};


// [[Rcpp::export]]
NumericVector pbd_integrate_odeint(NumericVector ry,
                                   NumericVector times,
                                   NumericVector pars,
                                   double atol,
                                   double rtol,
                                   std::string stepper)
{
  std::vector<double> y(ry.size() + 2, 0.0);    // [0,y,0]
  std::copy(ry.begin(), ry.end(), y.begin() + 1);

  auto rhs_obj = ode_rhs(pars);
  odeint_helper::integrate(stepper, std::ref(rhs_obj), y, times[0], times[1], 0.1 * (times[1] - times[0]), atol, rtol);
  return NumericVector(y.cbegin() + 1, y.cend() - 1);
}
