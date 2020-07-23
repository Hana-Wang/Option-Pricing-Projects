#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

double f(double x)
{
    double pi = M_PI;
    return exp(-pow(x,2) * 1.0/2.0)/sqrt(2.0 * M_PI);
}

// using composite simpson rule to find the definite integral
// with finite limit $a$ and $b$
// $n$ is the number of subintervals between lower limit $a$ and upper limmit $b$
double simpson(double a, double b, int n)
{
    // Construct a vector with $n+1$ elements with the value 0.0.
    vector<double> X(n+1, 0.0);
    vector<double> Y(n+1, 0.0);
    double h = (b-a)/double(n);

    for(int i=0; i <= n; i++)
    {
      X[i] = a + i * h;
      Y[i] = f(X[i]);
    }

    double sum = 0;
    for(int i=1; i <= n/2; i++)
    {
       sum += Y[2*i-2] + 4 * Y[2*i-1] + Y[2*i];
    }

    return sum*h/3;
}


// Cumulative Normal Distribution Function of N(0,1) by Simpson rule
double normalDistribution(double x)
{
  // -10.0 is a very small number for the distribution N(0,1)
  if (x < -10.0) return 0.0;
  if (x > 10.0) return 1.0;

  // number of steps;
  // number of subintervals you would choose to approximate integral
  int N = 200;

  //range of integration
  double a=0, b=x;
  double sum = simpson(a, b, N)+(1.0/2.0);
  return sum;
}


// the stock price following a Geometric Brownian Motion,
// which is also a lognormal distribution
// Dynamics of the stock is dS(t) = mu * S(t) dt + sigma * dW(t)
// When using risk neutral pricing formula to calculate options on the stocks,
// mu becomes a risk-free interest rate, denoted by r.

//Black-Scholes Call option price formula
// S(t)N(d_+)-KN(d_-)
// $tau=T-t$ is the time to option maturity
double euroCall(double S, double K, double sigma, double r, double t, double T)
{
   double d2 = (log(S/K)+(r-pow(sigma, 2)/2.0)*(T-t))/(sigma*sqrt(T-t));
   //cout << d2 <<endl;
   //double d1 = d2 + sigma*sqrt(T-t);
   double d1 = (log(S/K)+(r+pow(sigma, 2)/2.0)*(T-t))/(sigma * sqrt(T-t));
   //cout << d1 <<endl;

   double call = S*normalDistribution(d1)-K*exp(-r*(T-t))*normalDistribution(d2);

   return call;
}

int main()
{

  double S = 100;  //stock price at the calendar time $t$
  double K = 90;   // strike price of the European call option
  double sigma = 0.02;  // stock price volatility
  double r = 0.05;  //risk-free interest rate
  double t = 0;
  double T = 0.17;
  double tau=0.17; //Time to maturity

  cout << "European call on Black-Scholes Geometric Brownian Motion Stocks:"
     <<endl<< euroCall(S, K, sigma, r, t, T) << endl;


  return 0;

}
