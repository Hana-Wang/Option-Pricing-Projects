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

// the T-bond price following a Geometric Brownian Motion,
// which is also a lognormal distribution
// Dynamics of the stock is dB(t, T) = mu(t) * B(t, T) dt + sigma(t,T) * B(t, T) dW(t)
// When using risk neutral pricing formula to calculate options on the bonds,
// mu becomes a risk-free interest rate, denoted by r.
// the forward bond price with future time point $T0$, $T$ following a Geometric Brownian Motion
// In the T0-forward measure, the dynamics of the forward bond with future time point $T0$, $T$
// is dB(t; T0, T) = [sigma(t,T)-sigma(t, T0)] * B(t; T0, T) dW_T0(t),
// where dW_T0(t) is a Brownian motion under the T0-forward measure.

// T-Bond volatility at time t is calculated by sigma(t, T) = - sigma * G(t, T), with G(t, T) = e^{-kappa * (T-t)}
// where sigma is the short rate volatility; in Vasicek model, sigma is a constant.
// in Hull-White model, sigma is sigma(t), depending on time parameter.
// sigma(t, T) = - sigma(t) * G(t, T), with G(t, T) = e^{-kappa * (T-t)}
// Calculate T-Bond volatility using short rate volatility(sigma), kappa, tau=T-t,

// T-Bond volatility at time $t$, with tau = T-t
double BondSigmaVaciseck(double kappa, double sigma, double tau)
{
  double BondSigma = sigma * (1-exp(-kappa * tau)/-kappa;
  return BondSigma;
}

// the forward bond price with future time point $T0$, $T$ following a Geometric Brownian Motion
// In the T0-forward measure, the dynamics of the forward bond with future time point $T0$, $T$
// is dB(t; T0, T) = [sigma(t,T)-sigma(t, T0)] * B(t; T0, T) dW_T0(t),
// where dW_T0(t) is a Brownian motion under the T0-forward measure.
// Calculate the volatility of the forward bond price with future time point $T0$, $T$
// using T-Bond volatility sigma(t, T) and T0-Bond volatility sigma(t, T0)
// or using the input arguments: short rate volatility(sigma), kappa, tau=T-t, tau0=T0-t

// Forward Bond volatility with future time $T0$, $T$ at time $t$. tau0=T0-t, tau=T-t
double ForwardBondSigmaVaciseck(double kappa, double sigma, double tau0, double tau)
{
  double T0BondSigma = BondSigmaVaciseck(kappa, sigma, tau0);
  double TBondSigma  = BondSigmaVaciseck(kappa, sigma, tau);

  double FBondSigma = TBondSigma - T0BondSigma;
  return FBondSigma;
}

double ForwardBondVariance(double kappa, double sigma, double tau0, double tau)
{
  // distance between time T and time T0
  double TTODistance = tau-tau0;

  // Integral G(kappa, T0, T)
  double TT0Integral = (1-exp(-kappa * TTODistance)/-kappa;

  // Integral G(2*kappa, t, T0)
  double T0tIntegral = (1-exp(-2*kappa * tau0)/-kappa;

  // Forward bond price variance
  double FBondVar = pow(sigma,2) * pow(TT0Integral,2) * T0tIntegral;

  return FBondVar;
}
// the forward bond price with future time point $T0$, $T$, $T0<T$ following a Geometric Brownian Motion
// In the T0-forward measure, the dynamics of the forward bond with future time point $T0$, $T$
// is dB(t; T0, T) = [sigma(t,T)-sigma(t, T0)] * B(t; T0, T) dW_T0(t),
// where dW_T0(t) is a Brownian motion under the T0-forward measure.
// The forward bond price with future time point $T0$, $T$ follows a lognormal distribution,
// with the variance Var(t, T0) = the integral of (sigma(t,T)-sigma(t, T0))^2 from time t to time S.
// after analytical calculation,
// in the Vaseicek model, Var(t, T0, T) = sigma^2/(2*(kappa^3)) * (1-exp(-kappa * (T-T0))) * (1-exp(-2* kappa * (T0-t)))
// in the Hull-White model, Var(t, T0, T) = ?????
// calculate the variance of log of a forward bond price with with future time point $T0$, $T$, $T0<T$
// using the input arguments: short rate volatility(sigma), kappa, tau=T-t, tau1=T0-t,
// notice that T-T0 = tau - tau1


// Calculate T-Bond forward bond price with another future time T0 at time 0
double ForwadBondPV(double r0, double tau, double tau0)
{
  double TBond = exp(-r0 * tau);
  double T0Bond =  exp(-r0 * tau0);

  double FBondPV = TBond/T0Bond;
  return FBondPV;
}

// calculate the forward bond price process
// B(S; S, T) = B(0; S,T) exp((r-(pow(sigma,2)/2))*tau + normalRandom() * sigma * sqrt(tau))

double randPath(double kappa, double sigma, double r0, double tau0, double tau)
{
  double F0 = ForwadBondPV(r0, tau, tau0);

  double FBondVar = ForwardBondVariance(kappa, sigma, tau0, tau);
  double FBondSDev = sqrt(FBondVar);

  double ST = S * exp(-(FBondVarV/2) - normalRandom() *  sqrt(FBondVarV));
  return ST;
}

// FBondSigma is the volatility of the forward bond price process
// tau is the time to bond option maturity
// r0 is the quoted bond yield on the calendar time $t$
// F0 is the forward bond price with future time $S$, $T$, $S<T$ at the calendar time $t$; B(t; T0,T)
double randPath1(double F0, double FBondSigma, double tau)
{
  // Forward bond process variance
  double FBondVar = pow(FBondSigma, 2)*tau;
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
double euroCall1(double F0, double FBondSigma, double tau, double r0, double K, int nrSim)
{
  double sumCalls = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPath1( F0, FBondSigma, tau);
     double call = callT(path, K);
     sumCalls+=call;
  }

  double mean = sumCalls/(double)nrSim;

  double T0Bond =  exp(-r0 * tau);
  double PV = T0Bond * mean;
  //double PV = mean*T0Bond*exp(-r0*tau);

  //double PV = mean/pow((1+r), tau);
  return PV;
}

double euroCall(double kappa, double sigma, double r0, double tau, double tau1, double K, int nrSim)
{
  double sumCalls = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPath( kappa, sigma, r0, tau, tau1);
     double call = callT(path, K);
     sumCalls+=call;
  }

  double mean = sumCalls/(double)nrSim;
  double T0Bond =  exp(-r0 * tau1);
  double PV = mean*T0Bond*exp(-r0*tau1);
  //double PV = mean/pow((1+r), tau);
  return PV;
}


int main()
{
/*  double r0 = 0.05;
  double T = 20; //Bond maturity
  double T0 = 10; // Option maturity
  double tau1= 10; //Time to option maturity
  double tau = 20; //Time to bond maturity
  double sigma = 0.01;
  double kappa = 0.86;
  double K=*/
  double tau = 10.0/12; // Time to Bond option maturity
  double F0 = 939.68; //Forward bond price at the calendar time $t=0$
  double K = 1008.33;
  double FBondSigma = 0.09; // Forward bond price volatility
  double r0 = 0.1;
  int nrSim = 1000000; // number of samples

  srand(time(0));

  cout << "European call on T-Bond price with option maturity T0 is "
  << euroCall1(F0, FBondSigma, tau, r0, K, nrSim) << endl;



  return 0;

}
