// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

#define twoPi 6.28318573072
using namespace Rcpp;
using namespace arma;

//' Internal function: to estimate M, A and beta initial parameters also returns residual sum of squared (RSS).
//'
//' @param alphaOmegaParameters vector of the parameters alpha and omega.
//' @param vData: data to be fitted an FMM model.
//' @param timePoints: one single period time points.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector step1FMMrcpp(const arma::colvec& alphaOmegaParameters,
                                 const arma::colvec& vData,
                                 const arma::colvec& timePoints){

  double alphaParameter = alphaOmegaParameters[0];
  double omegaParameter = alphaOmegaParameters[1];

  arma::colvec mobiusTerm = 2*atan(omegaParameter*tan((timePoints - alphaParameter)/2));
  arma::colvec tStar = alphaParameter + mobiusTerm;

  arma::colvec costStar = cos(tStar);
  arma::colvec sentstar = sin(tStar);
  arma::colvec onesVector(size(vData));
  arma::mat designMatrix = join_horiz(onesVector.ones(), costStar, sentstar);
  arma::colvec coefCosinorModel = arma::solve(designMatrix, vData);

  double mParameter = coefCosinorModel(0); // intercept
  double cosCoeff = coefCosinorModel(1); // cos coefficient
  double sinCoeff = coefCosinorModel(2); // sin coefficient
  double phiEst = atan2(-sinCoeff, cosCoeff); // acrophase (phi)
  double aParameter = sqrt(pow(cosCoeff,2) + pow(sinCoeff,2));
  double betaParameter = fmod(phiEst + alphaParameter, twoPi);

  arma::colvec  mobiusModel = mParameter + aParameter*cos(betaParameter + mobiusTerm); // Mobius regression
  double residualSS = accu(pow(vData - mobiusModel,2))/vData.size(); //residual sum of squares

  return(Rcpp::NumericVector::create(mParameter,aParameter,alphaParameter,betaParameter,
                                     omegaParameter,residualSS));
}

