#include <cstdlib>
#include <cmath>
#include <ctime>
/*#include <math.h>
#include <stdlib.h>
#include <time.h>*/
#include <iostream>

using namespace std;

// return a uniformly distributed random number

double uniformRandom()
{
 return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

// return a normally distributed random number using Box-Muller Transform Sampling
double normalRandom()
{
 double u1=uniformRandom();
 double u2=uniformRandom();
 return cos(2.*M_PI*u2)*sqrt(-2.*log(u1));
}


// T-Bond volatility at time $t$, with tau = T-t
double BondSigmaVaciseck(double kappa, double sigma, double tau)
{
  double BondSigma = sigma * (1-exp(-kappa * tau))/kappa;
  return BondSigma;
}

// Forward Bond volatility with future time $T0$, $T$ at time $t$. tau0=T0-t, tau=T-t
double ForwardBondSigmaVaciseck(double kappa, double sigma, double tau0, double tau)
{
  double T0BondSigma = BondSigmaVaciseck(kappa, sigma, tau0);
  double TBondSigma  = BondSigmaVaciseck(kappa, sigma, tau);

  double FBondSigma = TBondSigma - T0BondSigma;
  return FBondSigma;
}

// Forward Bond variance with future time $T0$, $T$ at time $t$. tau0=T0-t, tau=T-t
double ForwardBondVariance(double kappa, double sigma, double tau0, double tau)
{
  // distance between time T and time T0
  double TTODistance = tau-tau0;

  // Integral G(kappa, T0, T)
  double TT0Integral = (1-exp(-kappa * TTODistance))/kappa;

  // Integral G(2*kappa, t, T0)
  double T0tIntegral = (1-exp(-2*kappa * tau0))/(2*kappa);

  // Forward bond price variance
  double FBondVar = pow(sigma,2) * pow(TT0Integral,2) * T0tIntegral;

  return FBondVar;
}

// Calculate T-Bond forward bond price with another future time T0 at time $t$
double ForwadBondPV(double r0, double tau, double tau0)
{
  double TBond = exp(-r0 * tau);
  double T0Bond =  exp(-r0 * tau0);

  double FBondPV = TBond/T0Bond;
  return FBondPV;
}

double randPath1(double F0, double FBondVar, double tau)
{
  // Forward bond process variance
  //double FBondVar = pow(FBondSigma, 2)*tau;
  // Forward bond process standard deviation
  double FBondSDev = sqrt(FBondVar);

  double ST = F0 * exp(-(FBondVar/2) - normalRandom() *  FBondSDev);
  return ST;
}

// Calculate the forward bond call option payoff at option maturity T0
double callT(double ST, double K)
{
    return max((ST-K), 0.0);
}

// FBondSigma is the volatility of the forward bond price process
// tau is the time to bond option maturity
// F0 is the forward bond price with future time $S$, $T$, $S<T$ at the calendar time $t$; B(t; T0,T)
// r0 is the quoted bond yield on the calendar time $t$
// Assuming the calendar time $t=0$
//double euroCall1(double F0, double FBondSigma, double tau, double r0, double K, int nrSim)
double euroCall(double r0, double kappa, double sigma, double tau0, double tau, double K, int nrSim)
{
  //Forward bond price at time $t=0$;
  double F0 = ForwadBondPV(r0, tau, tau0);
  // Forward bond variance at time $t=0$;
  double FBondVar = ForwardBondVariance(kappa, sigma, tau0, tau);

  //double test = pow(sigma, 2)/(2*pow(kappa,3))*pow((1-exp(-kappa * (tau))),2)*(1-exp(-2*kappa*tau0));

  double sumCalls = 0.0;
  for(int i=0; i<nrSim; i++)
  {
     //double path = randPath1( F0, test, tau);
     double path = randPath1( F0, FBondVar, tau);
     double call = callT(path, K);
     sumCalls+=call;
  }

  double mean = sumCalls/(double)nrSim;

  double T0Bond =  exp(-r0 * tau0);
  double PV = T0Bond * mean;
  //double PV = mean*T0Bond*exp(-r0*tau);

  //double PV = mean/pow((1+r), tau);
  return PV;
}


int main()
{
  double r0 = 0.03;
  double t = 0;
  double T = 2; //Bond maturity
  double T0 = 1; // Option maturity
  double tau = 2; //Time to bond maturity
  double tau0= 1; //Time to option maturity

  double sigma = 0.1;  // interest rate volatility
  double kappa = 0.1;
  double K = 0.8;

  //double F0 = ForwadBondPV(r0, tau, tau0); //Forward bond price at time $t=0$;
  // Forward Bond variance
  //double FBondVar = ForwardBondVariance(kappa, sigma, tau0, tau);


  //cout<<"Forward Bond standard deviation"<<sqrt(FBondVar)<<endl;

  //double test = pow(sigma, 2)/(2*pow(kappa,3))*pow((1-exp(-kappa * (tau))),2)*(1-exp(-2*kappa*tau0));
  //cout <<"test="<<test<<"sqrt="<<sqrt(test)<<endl;

  /*
  double tau = 10.0/12; // Time to Bond option maturity
  double F0 = 939.68; //Forward bond price at the calendar time $t=0$
  double K = 1008.33;
  double FBondSigma = 0.09; // Forward bond price volatility
  double r0 = 0.1;*/
  int nrSim = 1000000; // number of samples

  srand(time(0));

  cout << "European call on T-Bond price with option maturity T0 Using Monte Carlo Method is "
  << euroCall(r0, kappa, sigma, tau0, tau, K, nrSim) << endl;




  return 0;

}
