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
// the stock price following a Geometric Brownian Motion,
// which is also a lognormal distribution
// Dynamics of the stock is dS(t) = mu * S(t) dt + sigma * dW(t)
// When using risk neutral pricing formula to calculate options on the stocks,
// mu becomes a risk-free interest rate, denoted by r.
double randPath(double S, double sigma, double r, double tau)
{
  double ST = S *exp((r-(pow(sigma,2)/2))*tau + normalRandom() * sigma * sqrt(tau));
  return ST;
}

// the stock price following a Brownian Motion, which is also a normal distribution
// Dynamics of the stock is dS(t) = mu * dt + sigma * dW(t)
// When using risk neutral pricing formula to calculate options on the stocks,
// mu becomes a risk-free interest rate, denoted by r.
double randPathB(double S, double sigma, double r, double tau)
{
  double ST = S + r*tau + sigma * sqrt(tau) * normalRandom();
  return ST;
}


double callT(double ST, double K)
{
    return max((ST-K), 0.0);
}

double putT(double ST, double K)
{
   return max((K-ST), 0.0);
}

double digitalCallT(double ST, double K)
{
    double DT;
    if (ST > K){
          DT = 1;
    }
    else{
          DT = 0;
    }
    return DT;
}

double digitalPutT(double ST, double K)
{
    double DT;
    if (ST < K){
          DT = 1;
    }
    else{
          DT = 0;
    }
    return DT;
}

double euroCall(double S, double sigma, double r, double tau, double K, int nrSim)
{
  double sumCalls = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPath(S, sigma, r, tau);
     double call = callT(path, K);
     sumCalls+=call;
  }

  double mean = sumCalls/(double)nrSim;
  double PV = mean*exp(-r*tau);
  //double PV = mean/pow((1+r), tau);
  return PV;
}

double euroPut(double S, double sigma, double r, double tau, double K, int nrSim)
{

  double sumPuts = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPath(S, sigma, r, tau);
     double put = putT(path, K);
     sumPuts+=put;
   }

  double mean = sumPuts/(double)nrSim;
  double PV = mean/exp(r*tau);
  //double PV = mean/pow((1+r), tau);
  return PV;
}

double euroDigitalCall(double S, double sigma, double r, double tau, double K, int nrSim)
{
  double sumCalls = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPath(S, sigma, r, tau);
     double call = digitalCallT(path, K);
     sumCalls+=call;
  }

  double mean = sumCalls/(double)nrSim;
  //double mean = sumCalls/nrSim;
  double PV = mean*exp(-r*tau);
  //double PV = mean/pow((1+r), tau);
  return PV;
}

double euroDigitalPut(double S, double sigma, double r, double tau, double K, int nrSim)
{
  double sumCalls = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPath(S, sigma, r, tau);
     double call = digitalPutT(path, K);
     sumCalls+=call;
  }

  double mean = sumCalls/(double)nrSim;
  //double mean = sumCalls/nrSim;
  double PV = mean*exp(-r*tau);
  //double PV = mean/pow((1+r), tau);
  return PV;
}



// The options on stocks following Brownian Motion

double euroCallB(double S, double sigma, double r, double tau, double K, int nrSim)
{
  double sumCalls = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPathB(S, sigma, r, tau);
     double call = callT(path, K);
     sumCalls+=call;
  }

  double mean = sumCalls/(double)nrSim;
  //double mean = sumCalls/nrSim;
  double PV = mean*exp(-r*tau);
  //double PV = mean/pow((1+r), tau);
  return PV;
}

double euroPutB(double S, double sigma, double r, double tau, double K, int nrSim)
{

  double sumPuts = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPathB(S, sigma, r, tau);
     double put = putT(path, K);
     sumPuts+=put;
   }

  double mean = sumPuts/(double)nrSim;
  double PV = mean/exp(r*tau);
  //double PV = mean/pow((1+r), tau);
  return PV;
}



double euroDigitalCallB(double S, double sigma, double r, double tau, double K, int nrSim)
{
  double sumCalls = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPathB(S, sigma, r, tau);
     double call = digitalCallT(path, K);
     sumCalls+=call;
  }

  double mean = sumCalls/(double)nrSim;
  //double mean = sumCalls/nrSim;
  double PV = mean*exp(-r*tau);
  //double PV = mean/pow((1+r), tau);
  return PV;
}

double euroDigitalPutB(double S, double sigma, double r, double tau, double K, int nrSim)
{
  double sumCalls = 0.0;

  for(int i=0; i<nrSim; i++)
  {
     double path = randPathB(S, sigma, r, tau);
     double call = digitalPutT(path, K);
     sumCalls+=call;
  }

  double mean = sumCalls/(double)nrSim;
  //double mean = sumCalls/nrSim;
  double PV = mean*exp(-r*tau);
  //double PV = mean/pow((1+r), tau);
  return PV;
}


int main()
{
  double S = 100;
  double K = 90;
  double sigma = 0.02;
  double r = 0.05;
  double tau=0.17; //Time to maturity (say, 3months)
  int nrSim = 9;


  srand(time(0));

  cout << "European call on Black-Scholes Geometric Brownian Motion Stocks:" << euroCall(S, sigma, r, tau, K, nrSim) << endl;
  cout << "European put on Black-Scholes Geometric Brownian Motion Stocks:" << euroPut(S, sigma, r, tau, K, nrSim) << endl;
  cout << "European Digital call on Black-Scholes Geometric Brownian Motion Stocks:"
       << euroDigitalCall(S, sigma, r, tau, K, nrSim) << endl;
  cout << "European Digital call on Black-Scholes Geometric Brownian Motion Stocks:"
       << euroDigitalPut(S, sigma, r, tau, K, nrSim) << endl;



  cout << "European call on Brownian Motion Stocks:" << euroCallB(S, sigma, r, tau, K, nrSim) << endl;
  cout << "European call on Brownian Motion Stocks:" << euroPutB(S, sigma, r, tau, K, nrSim) << endl;
  cout << "European Digital call on Brownian Motion Stocks:"
       << euroDigitalCallB(S, sigma, r, tau, K, nrSim) << endl;
  cout << "European Digital call on on Brownian Motion Stocks:"
       << euroDigitalPutB(S, sigma, r, tau, K, nrSim) << endl;


  return 0;

}
