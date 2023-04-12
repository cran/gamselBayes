#ifndef APPROXLOGML_H
   #define APPROXLOGML_H
   double approxLogML(double muqBetaZero, double sigsqqBeta0, arma::vec muqgammaBeta, double AqrhoBeta,
                      double BqrhoBeta, double muqrecipSigsqBeta, arma::vec muqbBeta, arma::vec muqBetaTilde,
                      arma::mat SigmaqBetaTilde, double muqrecipaBeta, double lambdaqSigsqBeta,
                      double sBetaHYP, double lambdaqaBeta, arma::mat muqgammaU, arma::vec AqrhoU,
                      arma::vec BqrhoU, arma::vec muqrecipSigsqU, arma::vec muqbU,arma::mat muqUtilde,
                      arma::mat sigsqqUtilde, arma::vec muqrecipaU, arma::vec lambdaqSigsqU,double sUHYP, 
                      arma::vec lambdaqaU, double muqrecipaEps, double muqrecipSigsqEps, double lambdaqaEps, 
                      double lambdaqSigsqEps, double sEpsHYP, int n, arma::uvec ncZvec,int familyNum, 
                      arma::vec ySign, arma::vec omega20);
#endif



