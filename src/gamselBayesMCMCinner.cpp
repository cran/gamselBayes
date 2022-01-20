/********** C++ function: gamselBayesMCMCinner **********/

/* Carries out Markov chain Monte Carlo for a gamsel-type model. */

/* Last changed: 06 DEC 2021 */

#include <RcppArmadillo.h>
#include "printPercMsgs.h"
#include "drawInvGaussVec.h"
#include "omitVecEnt.h"
#include "rTruncNormPos.h"

using namespace Rcpp;

// [[Rcpp::export]]

List gamselBayesMCMCinner(arma::vec y, arma::mat X, arma::mat Z,int familyNum,
                          arma::uvec ncZvec, int ncZmax, int dGeneral, arma::uvec ZsttInds, 
                          arma::uvec ZendInds, arma::vec XTy, arma::mat XTX, arma::vec ZTy,
                          arma::mat ZTX, arma::mat ZTZ,double sigmaBeta0HYP, double sepsHYP, 
                          double sbetaHYP, double suHYP, double AbetaHYP, double BbetaHYP,
                          double AuHYP, double BuHYP, int numMCMC, int msgCode) 
{
   /* Declare all non-input variables: */

   int n = y.n_elem;
   int ncX = X.n_cols;
   int ncZ = Z.n_cols;
   int percCnt; 
   arma::mat wZmat(ncZmax,dGeneral);
   arma::vec ySign(n);
   arma::vec beta0(numMCMC);
   arma::mat betaTilde(ncX,numMCMC); 
   arma::mat gammaBeta(ncX,numMCMC); 
   arma::mat bBeta(ncX,numMCMC); 
   arma::cube uTilde(ncZmax,dGeneral,numMCMC);
   arma::cube gammaU(ncZmax,dGeneral,numMCMC);
   arma::vec recipSigsqBeta(numMCMC);
   arma::vec recipaBeta(numMCMC);
   arma::vec recipSigsqEps(numMCMC); 
   arma::vec recipaEps(numMCMC);
   double meanBeta0; 
   double sdBeta0;
   arma::mat Omega(ncX,ncX);
   arma::mat UOmega(ncX,ncX);
   arma::mat VOmega(ncX,ncX);
   arma::vec dOmega(ncX);  
   arma::vec z1(ncX);
   arma::vec z2(ncZmax);
   double scaleVal;
   double shapeVal;
   double shapeValDash;
   arma::vec rhoBeta(numMCMC); 
   arma::mat rhoU(dGeneral,numMCMC);
   double probVal; 
   arma::mat bU(dGeneral,numMCMC);
   arma::mat recipSigsqU(dGeneral,numMCMC);
   arma::mat recipaU(dGeneral,numMCMC);
   arma::vec tmpSix(dGeneral);
   arma::vec cAux(n);
   double yT1adj;
   arma::vec XTyAdj(ncX);
   arma::vec ZTyAdj(ncZ); 
   arma::mat ZTXj(ncZmax,ncX);
   arma::mat ZTXjd(ncZmax,ncX);
   arma::mat ZTZjjd(ncZmax,ncZmax);
   arma::mat Zj(n,ncZmax);
   double omega1;
   double omega2;
   arma::vec omega3(ncX);
   double omega4;
   double omega5;
   arma::vec omega6(ncZmax);
   arma::vec omega7(ncZmax);
   arma::vec omega8(ncZmax);
   double omega9;
   arma::vec omega10(n);
   arma::vec firstTermBeta(ncX);
   arma::vec secondTermBeta(ncX);
   arma::vec invGaussParmBeta(ncX);
   arma::vec firstTermU(ncZmax);
   arma::vec secondTermUFac1; 
   arma::vec secondTermUFac2;
   arma::vec secondTermU(ncZmax);
   arma::vec invGaussParmU(dGeneral);
   double quadTerm;
   arma::vec betaCurr(ncX);
   arma::mat uCurr(ncZmax,dGeneral);
   arma::mat uTildeCurr(ncZmax,dGeneral);
   arma::mat gammaUcurr(ncZmax,dGeneral);
   arma::vec XTXej(ncX);
   arma::uvec Kminus1(dGeneral);

   /* Initialise chains: */
   
   beta0 = arma::ones(numMCMC);
   betaTilde = arma::zeros(ncX,numMCMC);
   bBeta = arma::ones(ncX,numMCMC);
   recipSigsqBeta = arma::ones(numMCMC);
   recipaBeta = arma::ones(numMCMC);
   gammaBeta = 0.5*arma::ones(ncX,numMCMC);
   rhoBeta = 0.5*arma::ones(numMCMC);
   recipSigsqEps = arma::ones(numMCMC);
   recipaEps = arma::ones(numMCMC);

   if (dGeneral>0)
   {
      bU = arma::ones(dGeneral,numMCMC);
      recipSigsqU = arma::ones(dGeneral,numMCMC);
      recipaU = arma::ones(dGeneral,numMCMC);
      rhoU = 0.5*arma::ones(dGeneral,numMCMC);
      
      for (int g=0; g<numMCMC; g++)
      {
         uTilde.slice(g) = arma::zeros(ncZmax,dGeneral);
         gammaU.slice(g) = 0.5*arma::ones(ncZmax,dGeneral);
      }
   }

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

   /* Perform Markov chain Monte Carlo sampling: */

   percCnt = 0; 
   if (msgCode>0)
   {
      Rcout << "   The percentage of Markov chain Monte Carlo sampling completed is:" <<"\n";
      Rcout << "   ";
   };

   for (int g = 1 ; g < numMCMC; g++)
   {

      /* Print percentage progess message: */

      percCnt = printPercMsgs(msgCode,numMCMC,g+1,percCnt); 

      /* Draw sample from the beta_0 full conditional distribution: */

      omega1 = yT1adj;
      omega2 = n*recipSigsqEps(g-1) + (1.0/(sigmaBeta0HYP*sigmaBeta0HYP));
      meanBeta0 = omega1*recipSigsqEps(g-1)/omega2;
      sdBeta0 = 1.0/sqrt(omega2);
      beta0(g) = R::rnorm(meanBeta0,sdBeta0);

      /* Draw sample from the betaTilde full conditional distribution: */

      Omega = recipSigsqEps(g-1)*((gammaBeta.col(g-1)*gammaBeta.col(g-1).t())%XTX);
      Omega = Omega + recipSigsqBeta(g-1)*diagmat(bBeta.col(g-1));
      omega3 = XTyAdj;
      if (dGeneral>0)
      {
         for (int j=0; j<dGeneral; j++)
         {
            ZTXj = ZTX.rows(ZsttInds(j),ZendInds(j));
            omega3 = omega3 - ZTXj.t()*((gammaU.slice(g-1).rows(0,Kminus1(j)).col(j))
                                        %(uTilde.slice(g-1).rows(0,Kminus1(j)).col(j)));
         }
      }     

      svd(UOmega,dOmega,VOmega,Omega);     

      z1 = rnorm(ncX); 
      firstTermBeta = (UOmega.t()*z1)/sqrt(dOmega);                         
      secondTermBeta  = recipSigsqEps(g-1)*(UOmega.t()*(gammaBeta.col(g-1)%omega3))/dOmega;   
      betaTilde.col(g) = UOmega*(firstTermBeta + secondTermBeta);                     

      /* Draw sample from the b_beta full conditional distribution: */
   
      invGaussParmBeta = 1.0/(sqrt(recipSigsqBeta(g-1))*abs(betaTilde.col(g)));
      bBeta.col(g) = drawInvGaussVec(invGaussParmBeta);  

      /* Draw sample from the 1/sigsq_beta full conditional distribution: */

      shapeVal = 0.5*(ncX + 1.0);    
      scaleVal = 1.0/(recipaBeta(g-1) + 0.5*sum(betaTilde.col(g)%bBeta.col(g)%betaTilde.col(g)));
      recipSigsqBeta(g) = R::rgamma(shapeVal,scaleVal);

      /* Draw sample from the 1/a_beta full conditional distribution: */

      shapeVal = 1.0;
      scaleVal = 1.0/(recipSigsqBeta(g) + (1.0/(sbetaHYP*sbetaHYP)));
      recipaBeta(g) = R::rgamma(shapeVal,scaleVal);

      /* Draw sample from the gamma_beta full conditional distribution: */
    
      betaCurr = gammaBeta.col(g-1)%betaTilde.col(g);
      if (dGeneral>0)
         uCurr = gammaU.slice(g-1)%uTilde.slice(g-1);

      for (int j = 0 ; j < ncX; j++)
      {  
         omega4 = XTyAdj(j); 
         XTXej = XTX.col(j);
         if (ncX>1)
            omega4 = omega4 - sum(omitVecEnt(XTXej,j)%omitVecEnt(betaCurr,j));

         if (dGeneral>0)
         {
            for (int jd = 0 ; jd < dGeneral; jd++)
            {  
               ZTXjd = ZTX.rows(ZsttInds(jd),ZendInds(jd));
               omega4 = omega4 - sum(ZTXjd.col(j)%uCurr.rows(0,Kminus1(jd)).col(jd));
            }
         } 
         quadTerm = betaTilde(j,g)*betaTilde(j,g)*XTX(j,j);
         omega5 = (log(rhoBeta(g-1)/(1.0-rhoBeta(g-1)))
                   - 0.5*recipSigsqEps(g-1)*(quadTerm - 2.0*betaTilde(j,g)*omega4));
         probVal = 1.0/(1.0 + exp(-omega5));
         gammaBeta(j,g) = R::rbinom(1,probVal);
      }  

      /* Draw sample from the rho_beta full conditional distribution: */

      shapeVal = AbetaHYP +  sum(gammaBeta.col(g));
      shapeValDash = BbetaHYP + ncX - sum(gammaBeta.col(g));
      rhoBeta(g) = R::rbeta(shapeVal,shapeValDash); 

      if (dGeneral>0)
      {
         /* Draw sample from the uTilde_j full conditional distributions: */

         betaCurr = gammaBeta.col(g)%betaTilde.col(g);
         uTildeCurr = uTilde.slice(g-1);

         for (int j = 0 ; j < dGeneral; j++)
         {  
            omega6 = ZTyAdj.rows(ZsttInds(j),ZendInds(j)) - ZTX.rows(ZsttInds(j),ZendInds(j))*betaCurr;
            for (int jd = 0 ; jd < dGeneral; jd++)                                                                   
            {                                                                                            
               ZTZjjd = ZTZ.rows(ZsttInds(j),ZendInds(j)).cols(ZsttInds(jd),ZendInds(jd));               
               omega6 = omega6 -  ZTZjjd*((gammaU.slice(g-1).rows(0,Kminus1(jd)).col(jd))
                                          %(uTildeCurr.rows(0,Kminus1(jd)).col(jd)));                            
            }                                                                                            
            omega6 = omega6 + (wZmat.rows(0,Kminus1(j)).col(j))%(gammaU.slice(g-1).rows(0,Kminus1(j)).col(j))
                                                               %(uTildeCurr.rows(0,Kminus1(j)).col(j));       
 
            omega7 = (recipSigsqEps(g-1)*(gammaU.slice(g-1).col(j).rows(0,Kminus1(j)))%(wZmat.col(j).rows(0,Kminus1(j)))
                                      + recipSigsqU(j,g-1)*bU(j,g-1));
            z2 = rnorm(ncZmax);
            firstTermU = z2.rows(0,Kminus1(j))/sqrt(omega7);
            secondTermUFac1 = recipSigsqEps(g-1)*gammaU.slice(g-1).rows(0,Kminus1(j)).col(j);
            secondTermUFac2 = omega6.rows(0,Kminus1(j))/omega7.rows(0,Kminus1(j));  
            secondTermU = secondTermUFac1%secondTermUFac2;
            uTildeCurr.rows(0,Kminus1(j)).col(j) =  firstTermU + secondTermU;
         }
         uTilde.slice(g) = uTildeCurr;

         /* Draw sample from the b_{uj} full conditional distributions: */
   
         for (int j = 0 ; j < dGeneral; j++)
            invGaussParmU(j) = 1.0/sqrt(recipSigsqU(j,g-1)*sum(uTilde.slice(g).col(j)%uTilde.slice(g).col(j)));  
         bU.col(g) = drawInvGaussVec(invGaussParmU);
        
         /* Draw samples from the sigsq_{uj} and a_{uj} full conditional distributions: */

         for (int j = 0 ; j < dGeneral; j++)
         {  
            shapeVal = 0.5*(ncZvec(j) + 1.0);
            scaleVal = 1.0/(recipaU(j,g-1) + 0.5*bU(j,g)*sum(uTilde.slice(g).col(j)%uTilde.slice(g).col(j)));
            recipSigsqU(j,g) = R::rgamma(shapeVal,scaleVal);  
            shapeVal = 1.0;
            scaleVal = 1.0/(recipSigsqU(j,g) + (1.0/(suHYP*suHYP)));
            recipaU(j,g) = R::rgamma(shapeVal,scaleVal);
         }

         /* Draw sample from the gamma_u full conditional distribution: */

         gammaUcurr = gammaU.slice(g-1);

         for (int j = 0 ; j < dGeneral; j++)
         {  
            omega8 = ZTyAdj.rows(ZsttInds(j),ZendInds(j)) - ZTX.rows(ZsttInds(j),ZendInds(j))*betaCurr;

            for (int jd = 0 ; jd < dGeneral; jd++)
            {
               ZTZjjd = ZTZ.rows(ZsttInds(j),ZendInds(j)).cols(ZsttInds(jd),ZendInds(jd));
               omega8 = omega8 - ZTZjjd*((gammaUcurr.rows(0,Kminus1(jd)).col(jd))
                                        %(uTilde.slice(g).rows(0,Kminus1(jd)).col(jd)));
            } 
            omega8 = omega8 + ((wZmat.rows(0,Kminus1(j)).col(j))
                                        %(gammaUcurr.rows(0,Kminus1(j)).col(j))
                                        %(uTilde.slice(g).rows(0,Kminus1(j)).col(j)));
  
            for (int k = 0; k < ncZvec(j); k++)
            {
               quadTerm = uTilde(k,j,g)*uTilde(k,j,g)*wZmat(k,j);
               omega9 = log(rhoU(j,g-1)/(1.0 - rhoU(j,g-1))) 
                       - 0.5*recipSigsqEps(g-1)*(quadTerm - 2.0*uTilde(k,j,g)*omega8(k));
               probVal = 1.0/(1.0 + exp(-omega9));
               gammaUcurr(k,j)  = R::rbinom(1,probVal);
            }
         }
         gammaU.slice(g) = gammaUcurr; 
            
         /* Draw sample from the rho_uj full conditional distributions: */

         for (int j = 0 ; j < dGeneral; j++)
         {  
            shapeVal = AuHYP + sum(gammaU.slice(g).rows(0,Kminus1(j)).col(j));
            shapeValDash = BuHYP + ncZvec(j) - sum(gammaU.slice(g).rows(0,Kminus1(j)).col(j));
            rhoU(j,g) = R::rbeta(shapeVal,shapeValDash); 
         }
      }

      omega10 = beta0(g)*arma::ones(n) + X*(gammaBeta.col(g)%betaTilde.col(g));

      if (dGeneral>0)
      {
         for (int j=0; j<dGeneral ; j++)
         {
            Zj = Z.cols(ZsttInds(j),ZendInds(j));
            omega10 = omega10 + Zj*((gammaU.slice(g).rows(0,Kminus1(j)).col(j))
                                   %(uTilde.slice(g).rows(0,Kminus1(j)).col(j)));
         }
      }

      if (familyNum == 1)
      {
         /* Draw sample from the 1/sigsq_epsilon full conditional distribution: */

         shapeVal = 0.5*(n + 1);
         scaleVal = 1.0/(recipaEps(g-1) + 0.5*sum((y-omega10)%(y-omega10)));
         recipSigsqEps(g) = R::rgamma(shapeVal,scaleVal);

         /* Draw sample from the 1/a_epsilon full conditional distribution: */

         shapeVal = 1.0;
         scaleVal = 1.0/(recipSigsqEps(g) + (1.0/(sepsHYP*sepsHYP)));
         recipaEps(g) = R::rgamma(shapeVal,scaleVal); 
      }

      if (familyNum == 2)
      {
         for (int i = 0 ; i < n; i++)
            cAux(i) = ySign(i)*rTruncNormPos(ySign(i)*omega10(i)); 

         yT1adj = sum(cAux)   ;  XTyAdj = X.t()*cAux ; 
         if (dGeneral>0)
         {
            ZTyAdj = Z.t()*cAux;
         }
      }
   }
   if (msgCode>0) {Rcout << "\n";};

   List MCMCsamples = List::create(
                           Named("beta0",beta0),
                           Named("betaTilde",betaTilde),
                           Named("gammaBeta",gammaBeta),
                           Named("recipSigsqBeta",recipSigsqBeta),
                           Named("rhoBeta",rhoBeta),
                           Named("uTilde",uTilde),
                           Named("gammaU",gammaU),
                           Named("rhoU",rhoU),
                           Named("recipSigsqU",recipSigsqU),
                           Named("recipSigsqEps",recipSigsqEps));

   return MCMCsamples;
}

/************ End of gamselBayesMCMCinner  ************/

