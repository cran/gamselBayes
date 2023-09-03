#ifndef APPROXLOGML_H
   #define APPROXLOGML_H
double approxLogML(double muqBetaZero, double sigsqqBeta0, double logitBeta, arma::vec muqgammaBeta, 
                      double muqrecipSigsqBeta, arma::vec muqbBeta, arma::vec muqBetaTilde,
                      arma::mat SigmaqBetaTilde, double muqrecipaBeta, double lambdaqSigsqBeta,
                      double sBetaHYP, double lambdaqaBeta, double logitU, arma::vec muqgammaU,
                      arma::vec muqrecipSigsqU, arma::vec muqbU,arma::mat muqUtilde,
                      arma::mat sigsqqUtilde, arma::vec muqrecipaU, arma::vec lambdaqSigsqU,double sUHYP, 
                      arma::vec lambdaqaU, double muqrecipaEps, double muqrecipSigsqEps, double lambdaqaEps, 
                      double lambdaqSigsqEps, double sEpsHYP, int n, arma::uvec ncZvec,int familyNum, 
                      arma::vec ySign, arma::vec omega20);
#endif



