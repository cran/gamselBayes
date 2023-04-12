/********** C++ Function gamselBayesMFVBinner **********/

/* Carries out mean field variational Bayes fitting of a 
   specified gamsel-type model. */

/* Last changed: 09 FEB 2022 */

#include <RcppArmadillo.h>
#include "printPercMsgs.h"
#include "omitVecEnt.h"
#include "zetad.h"
#include "approxLogML.h"

using namespace Rcpp;

//[[Rcpp::export]] 

List gamselBayesMFVBinner(arma::vec y, arma::mat X, arma::mat Z,int familyNum, arma::uvec ncZvec, 
                          int ncZmax, int dGeneral, arma::uvec ZsttInds, arma::uvec ZendInds,
                          arma::vec XTy, arma::mat XTX, arma::vec ZTy, arma::mat ZTX, arma::mat ZTZ,
                          double sigmaBeta0HYP, double sEpsHYP, double sBetaHYP, double sUHYP, 
                          double AbetaHYP, double BbetaHYP, double AuHYP,double BuHYP,
                          int maxIter, double toler, int msgCode) 
{
   /* Declare all non-input variables: */

   int n = y.n_elem;
   int ncX = X.n_cols;
   int ncZ = Z.n_cols;
   int percCnt; 
   arma::mat wZmat(ncZmax,dGeneral);
   arma::vec innVec(ncZmax);
   double muqBetaZero;
   double sigsqqBetaZero;
   arma::vec muqBetaTilde(ncX);
   arma::mat SigmaqBetaTildeInv(ncX,ncX);
   arma::mat SigmaqBetaTilde(ncX,ncX);
   arma::vec muqbBeta(ncX);
   arma::vec muqBeta(ncX);
   double muqrecipSigsqBeta;
   double kappaqSigsqBeta; 
   double lambdaqSigsqBeta;
   double muqrecipaBeta;
   double kappaqaBeta; 
   double lambdaqaBeta;
   double AqrhoBeta;
   double BqrhoBeta;
   double muqlogitrhoBeta;
   arma::vec muqgammaBeta(ncX);
   arma::mat muqUtilde(ncZmax,dGeneral);
   arma::mat muqu(ncZmax,dGeneral);
   arma::mat sigsqqUtilde(ncZmax,dGeneral);
   arma::vec muqbU(dGeneral);
   arma::vec muqrecipSigsqU(dGeneral);
   arma::vec kappaqSigsqU(dGeneral); 
   arma::vec lambdaqSigsqU(dGeneral);
   arma::vec muqrecipaU(dGeneral);
   double kappaqaU; 
   arma::vec lambdaqaU(dGeneral); 
   arma::vec AqrhoU(dGeneral);
   arma::vec BqrhoU(dGeneral);
   double muqlogitrhoUj;
   arma::mat muqgammaU(ncZmax,dGeneral);
   double muqrecipSigsqEps;
   double kappaqSigsqEps;
   double lambdaqSigsqEps;
   arma::mat innerMat(ncX,ncX);
   arma::mat trMat(ncX,ncX);
   arma::vec muGammauFac(ncZmax);
   double muqrecipaEps;
   double kappaqaEps;
   double lambdaqaEps;
   arma::vec muqcAux(n);
   arma::vec ySign(n);
   arma::mat Omega_gammaBeta1st(ncX,ncX); 
   arma::mat OmegaqgammaBeta(ncX,ncX);  
   arma::vec SigsqqBetaTilde(ncX); 
   bool converged;
   int itNum; 
   double omega12; 
   arma::vec omega13(ncX);
   arma::vec omega14(ncX);
   double omega15; 
   arma::vec omega16(ncZmax);
   double omega17;
   arma::vec omega18(ncZmax); 
   double omega19;
   arma::vec omega20(n);
   double etaCurr;
   double quadTerm;
   double yT1adj;
   arma::vec XTyAdj(ncX);
   arma::vec ZTyAdj(ncZ); 
   arma::mat ZTXj(ncZmax,ncX);
   arma::vec innerVec(ncX-1);
   arma::mat ZTZjjd(ncZmax,ncZmax);
   arma::mat Zj(n,ncZmax);
   arma::uvec Kminus1(dGeneral);
   arma::vec logMargLik(maxIter); 
   int stopType;
   int numIters;
   double relErr;

   /* Initialise q-density parameters: */

   muqBetaTilde = arma::zeros(ncX);
   muqbBeta = arma::ones(ncX);
   muqrecipSigsqBeta = 1.0;
   kappaqSigsqBeta = 0.5*(ncX+1);
   muqrecipaBeta = 1.0;
   kappaqaBeta = 1.0;
   muqgammaBeta = 0.5*arma::ones(ncX);
   muqrecipSigsqEps = 1.0;
   kappaqSigsqEps = 0.5*(n+1);
   muqrecipaEps = 1.0;
   kappaqaEps = 1.0;
  
   if (dGeneral>0)
   {
      muqUtilde = arma::zeros(ncZmax,dGeneral);
      sigsqqUtilde = arma::ones(ncZmax,dGeneral);
      muqbU = arma::ones(dGeneral);
      muqrecipSigsqU = arma::ones(dGeneral);
      muqrecipaU = arma::ones(dGeneral);
      kappaqaU = 1.0;
      lambdaqaU = arma::ones(dGeneral);
      AqrhoU = arma::ones(dGeneral);
      BqrhoU = arma::ones(dGeneral);
      muqgammaU = 0.5*arma::ones(ncZmax,dGeneral);
 
      for (int j=0; j<dGeneral; j++)
         kappaqSigsqU(j) = 0.5*(ncZvec(j) + 1.0);
   }

   /* Initialise the marginal log-likelihood vector and the
      number of iterations:                                   */

   logMargLik = arma::zeros(maxIter);
   numIters = maxIter;

   /* Do data matrix initialisations: */

   Kminus1 = ncZvec - 1;
   ySign = 2.0*y - 1.0;
   yT1adj = 0.0   ;   XTyAdj = XTy   ;   ZTyAdj = ZTy;

   if (dGeneral>0)
   {
      for (int j=0; j<dGeneral; j++)
      {
         wZmat.col(j) = arma::zeros(ncZmax);
         for (int k=0; k<ncZvec(j); k++)
            wZmat(k,j) = ZTZ.rows(ZsttInds(j),ZendInds(j)).cols(ZsttInds(j),ZendInds(j))(k,k);
      }
   }

   /* Perform mean field variational Bayes iterations: */

   percCnt = 0; 
   if (msgCode>0)
   {
      Rcout << "   The percentage of the maximum number of mean field variational" <<"\n";
      Rcout << "   Bayes iterations completed is:" <<"\n";
      Rcout << "   ";
   };  

   converged = FALSE;
   itNum = 0;
   while (!converged)
   { 
      itNum = itNum + 1;
    
      /* Print percentage progess message: */

      percCnt = printPercMsgs(msgCode,maxIter,itNum + 1,percCnt);

      /* Update the q(beta_0) parameters: */
    
      omega12 = yT1adj;
      sigsqqBetaZero = 1.0/(n*muqrecipSigsqEps+(1.0/(sigmaBeta0HYP*sigmaBeta0HYP)));
      muqBetaZero = sigsqqBetaZero*muqrecipSigsqEps*omega12;
    
      /* Update the q(betaTilde) parameters: */
    
      Omega_gammaBeta1st = diagmat(muqgammaBeta%(1.0-muqgammaBeta));
      OmegaqgammaBeta = Omega_gammaBeta1st + muqgammaBeta*muqgammaBeta.t();
      SigmaqBetaTildeInv = muqrecipSigsqEps*(OmegaqgammaBeta%XTX) + muqrecipSigsqBeta*diagmat(muqbBeta);
      SigmaqBetaTilde = SigmaqBetaTildeInv.i(); 
      omega13 = XTyAdj;
      if (dGeneral>0)
      {
         for (int j=0; j<dGeneral; j++)
         {
            ZTXj = ZTX.rows(ZsttInds(j),ZendInds(j));
            omega13 = omega13 - ZTXj.t()*(muqgammaU.col(j).rows(0,Kminus1(j))
                                          %muqUtilde.col(j).rows(0,Kminus1(j)));
         }
      }     

      muqBetaTilde = muqrecipSigsqEps*SigmaqBetaTilde*(muqgammaBeta%omega13);
    
      /* Update the q(b_beta) parameters: */   
    
      for (int j=0; j<ncX; j++)
         omega14(j) = SigmaqBetaTilde(j,j) + muqBetaTilde(j)*muqBetaTilde(j);
      
      muqbBeta = 1.0/sqrt(muqrecipSigsqBeta*omega14);
    
      /* Update the q(Sigsqbeta) parameters: */
    
      lambdaqSigsqBeta = muqrecipaBeta + 0.5*sum(muqbBeta%omega14);
      muqrecipSigsqBeta = kappaqSigsqBeta/lambdaqSigsqBeta;
    
      /* Update the q(a_beta) parameters: */
    
      lambdaqaBeta = muqrecipSigsqBeta + (1.0/(sBetaHYP*sBetaHYP));
      muqrecipaBeta = kappaqaBeta/lambdaqaBeta;
    
      /* Update the q(rho_beta) parameters: */
   
      AqrhoBeta = AbetaHYP + sum(muqgammaBeta);
      BqrhoBeta = BbetaHYP + ncX - sum(muqgammaBeta); 
    
      /* Update the q(gamma_beta) parameters: */
    
      muqlogitrhoBeta = R::digamma(AqrhoBeta) - R::digamma(BqrhoBeta);
      muqu = muqgammaU%muqUtilde;  
      for (int j=0 ; j<ncX ; j++)
      { 
         omega15 = XTyAdj(j);
         if (dGeneral>0)
         {
            for (int jd=0 ; jd<dGeneral ; jd++)   
               omega15 = omega15 - sum(ZTX.rows(ZsttInds(jd),ZendInds(jd)).
                                       col(j)%muqu.col(jd).rows(0,Kminus1(jd)));
         }

         omega15 = muqBetaTilde(j)*omega15; 
         if (ncX>1)
         {
            innerVec = omitVecEnt(SigmaqBetaTilde.col(j),j) + muqBetaTilde(j)*omitVecEnt(muqBetaTilde,j);
            omega15 = omega15 - sum((omitVecEnt(XTX.col(j),j)%omitVecEnt(muqgammaBeta,j))%innerVec);
         }
         quadTerm = ((muqBetaTilde(j)*muqBetaTilde(j)) + SigmaqBetaTilde(j,j))*XTX(j,j);
         etaCurr = muqlogitrhoBeta - 0.5*muqrecipSigsqEps*(quadTerm-2.0*omega15);
         muqgammaBeta(j) = 1.0/(1.0+exp(-etaCurr));
      }  
    
      if (dGeneral>0)
      {
         /* Update the q(u_Tilde) parameters: */
   
         muqu = muqgammaU%muqUtilde;  

         for (int j=0; j<dGeneral ; j++)
         {
            omega16 = ZTyAdj.rows(ZsttInds(j),ZendInds(j));
            omega16 = omega16 - ZTX.rows(ZsttInds(j),ZendInds(j))*(muqgammaBeta%muqBetaTilde);

            for (int jd=0; jd<dGeneral ; jd++)
            {
               ZTZjjd = ZTZ.rows(ZsttInds(j),ZendInds(j)).cols(ZsttInds(jd),ZendInds(jd));
               omega16 = omega16 - ZTZjjd*muqu.col(jd).rows(0,Kminus1(jd));

            }
            omega16 =  omega16 + (wZmat.rows(0,Kminus1(j)).col(j))%(muqu.rows(0,Kminus1(j)).col(j));   
            sigsqqUtilde.rows(0,Kminus1(j)).col(j) = 1.0/((muqrecipSigsqEps*
                                        (muqgammaU.rows(0,Kminus1(j)).col(j))%wZmat.rows(0,Kminus1(j)).col(j))
                                        + (muqrecipSigsqU(j)*muqbU(j)*arma::ones(ncZvec(j))));
            muqUtilde.rows(0,Kminus1(j)).col(j) = muqrecipSigsqEps*(muqgammaU.rows(0,Kminus1(j)).col(j)
                                                  %omega16)%sigsqqUtilde.rows(0,Kminus1(j)).col(j);
         }
    
         /* Update the q(b_uj), q(sigsq_uj) and q(a_uj) parameters: */
     
         for (int j=0; j<dGeneral ; j++)
         {
            omega17 = sum(muqUtilde.col(j)%muqUtilde.col(j)) + sum(sigsqqUtilde.col(j));
            muqbU(j) =  1.0/sqrt(muqrecipSigsqU(j)*omega17);
      
            lambdaqSigsqU(j) = muqrecipaU(j) + 0.5*muqbU(j)*omega17;
            muqrecipSigsqU(j) = kappaqSigsqU(j)/lambdaqSigsqU(j);
      
            lambdaqaU(j) = muqrecipSigsqU(j) + (1.0/(sUHYP*sUHYP));
            muqrecipaU(j) = kappaqaU/lambdaqaU(j);
         }
 
         /* Update the q(gamma_uj) parameters: */
    
         muqu = muqgammaU%muqUtilde;  
 
         for (int j=0; j<dGeneral ; j++)
         {
         
            omega18 = (ZTyAdj.rows(ZsttInds(j),ZendInds(j)) 
                       - ZTX.rows(ZsttInds(j),ZendInds(j))*(muqgammaBeta%muqBetaTilde));

            for (int jd=0; jd<dGeneral ; jd++)      
            {
               ZTZjjd = ZTZ.rows(ZsttInds(j),ZendInds(j)).cols(ZsttInds(jd),ZendInds(jd));
               omega18 = omega18 - ZTZjjd*muqu.rows(0,Kminus1(jd)).col(jd);
            }
            omega18 = omega18 + (wZmat.rows(0,Kminus1(j)).col(j))%(muqu.rows(0,Kminus1(j)).col(j));     

            AqrhoU(j) = AuHYP + sum(muqgammaU.rows(0,Kminus1(j)).col(j));
            BqrhoU(j) = BuHYP + ncZvec(j)- sum(muqgammaU.rows(0,Kminus1(j)).col(j));
            muqlogitrhoUj =  R::digamma(AqrhoU(j)) -  R::digamma(BqrhoU(j));

            for (int k=0; k<ncZvec(j) ; k++)      
            { 
               quadTerm = ((muqUtilde(k,j)*muqUtilde(k,j)) + sigsqqUtilde(k,j))*wZmat(k,j);
               omega19 = quadTerm - 2*muqUtilde(k,j)*omega18(k);
               etaCurr = muqlogitrhoUj - 0.5*muqrecipSigsqEps*omega19;
               muqgammaU(k,j) = 1.0/(1.0 + exp(-etaCurr));
            }   
         }
      }      

      /* Update the q(sigsq_epsilon) parameters and Albert-Chib auxiliary variable parameters: */

      omega20 = muqBetaZero*arma::ones(n) + X*(muqgammaBeta%muqBetaTilde);
      if (dGeneral>0)
      {
         for (int j=0; j<dGeneral ; j++)
         {
            Zj = Z.cols(ZsttInds(j),ZendInds(j));
            omega20 = omega20 + Zj*(muqgammaU.col(j).rows(0,Kminus1(j))
                                    %muqUtilde.col(j).rows(0,Kminus1(j)));
         }
      }         

      if (familyNum==1)
      {
         /* Update the q(sigsq_eps) parameters: */

         Omega_gammaBeta1st = diagmat(muqgammaBeta%(1.0-muqgammaBeta));
         OmegaqgammaBeta = Omega_gammaBeta1st + muqgammaBeta*muqgammaBeta.t();
         innerMat = OmegaqgammaBeta%(SigmaqBetaTilde + muqBetaTilde*muqBetaTilde.t());
         muqBeta = muqgammaBeta%muqBetaTilde;
         innerMat = innerMat - muqBeta*muqBeta.t();
         trMat = 0.5*XTX*innerMat;
            
         lambdaqSigsqEps = muqrecipaEps + 0.5*sum((y - omega20)%(y - omega20));
         lambdaqSigsqEps = lambdaqSigsqEps + 0.5*n*sigsqqBetaZero;

         for (int j=0; j<ncX ; j++)
            lambdaqSigsqEps = lambdaqSigsqEps + trMat(j,j);

         if (dGeneral>0)
         {
            for (int j=0;j<dGeneral;j++)
            {
               muGammauFac = muqgammaU.col(j)%(1.0-muqgammaU.col(j));
               innVec = (muqgammaU.col(j)%sigsqqUtilde.col(j)
                        + muGammauFac%(muqUtilde.col(j)%muqUtilde.col(j)));
               lambdaqSigsqEps = lambdaqSigsqEps + 0.5*sum(innVec%wZmat.col(j)); 
            }
         }

         muqrecipSigsqEps = kappaqSigsqEps/lambdaqSigsqEps;

         /* Update the q(a_Eps) parameters: */

         lambdaqaEps = muqrecipSigsqEps + (1.0/(sEpsHYP*sEpsHYP));
         muqrecipaEps = kappaqaEps/lambdaqaEps; 
      }

      if (familyNum==2)
      {
         for (int i=0;i<n;i++)
            muqcAux(i) = omega20(i) + ySign(i)*zetad(ySign(i)*omega20(i));

         muqrecipSigsqEps = 1;  

         yT1adj = sum(muqcAux)  ; XTyAdj = X.t()*muqcAux ;
         if (dGeneral>0)
         {
            ZTyAdj = Z.t()*muqcAux;
         }
      }

      /* Compute and store the current value of the approximate marginal log-likelihood: */

      logMargLik(itNum-1) = approxLogML(muqBetaZero,sigsqqBetaZero,muqgammaBeta,AqrhoBeta,BqrhoBeta, 
                                        muqrecipSigsqBeta,muqbBeta,muqBetaTilde,SigmaqBetaTilde, 
                                        muqrecipaBeta,lambdaqSigsqBeta,sBetaHYP,lambdaqaBeta,muqgammaU, 
                                        AqrhoU,BqrhoU,muqrecipSigsqU,muqbU,muqUtilde,sigsqqUtilde,
                                        muqrecipaU,lambdaqSigsqU,sUHYP,lambdaqaU,muqrecipaEps,
                                        muqrecipSigsqEps,lambdaqaEps,lambdaqSigsqEps,sEpsHYP,n,ncZvec,
                                        familyNum,ySign,omega20);

      /* Check for convergence: */
   
      if (itNum>=maxIter) 
      { 
         converged = TRUE;
         stopType = 1;
      }
      
      if (itNum>2)
      {
         relErr = fabs((logMargLik(itNum-1)/logMargLik(itNum-2)) - 1.0); 
         if (relErr<toler)
         {
            converged = TRUE;
            stopType = 0;
            numIters = itNum;
         }
      }
   } 
   if (msgCode>0) {Rcout << "\n";};
  
   /* Return list of q-density parameters: */

   List MFVBresults = List::create(Named("muqBetaZero",muqBetaZero),
                                   Named("sigsqqBetaZero",sigsqqBetaZero),
                                   Named("muqBetaTilde",muqBetaTilde),
                                   Named("SigmaqBetaTilde",SigmaqBetaTilde),
                                   Named("muqgammaBeta",muqgammaBeta),
                                   Named("kappaqSigsqBeta",kappaqSigsqBeta),
                                   Named("lambdaqSigsqBeta",lambdaqSigsqBeta),
                                   Named("AqrhoBeta",AqrhoBeta),
                                   Named("BqrhoBeta",BqrhoBeta),
                                   Named("muqUtilde",muqUtilde),
                                   Named("sigsqqUtilde",sigsqqUtilde),
                                   Named("AqrhoU",AqrhoU),
                                   Named("BqrhoU",BqrhoU),
                                   Named("muqgammaU",muqgammaU),
                                   Named("kappaqSigsqEps",kappaqSigsqEps),
                                   Named("lambdaqSigsqEps",lambdaqSigsqEps),
                                   Named("logMargLik",logMargLik),
                                   Named("numIters",numIters),
                                   Named("relErr",relErr),
                                   Named("stopType",stopType));
                                               
   return MFVBresults;
} 

/************ End of gamselBayesMFVBinner ************/




