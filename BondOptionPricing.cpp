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

//Bond European Call option price formula
// B(t,S)[B(t;S,T)N(d_+)-KN(d_-)] = B(t,T)N(d_+)-B(t,S)*K*N(d_-)
// $tau0=T0-t$ is the time to option maturity
double euroCall(double K, double r0, double sigma, double kappa, double t, double T, double T0)
{
   double tau = T-t;
   double tau0 = T0-t;
   // $T$-Bond price at time $t$
   double TBond = exp(-r0 * tau);
   //$T0$-Bond price at time $t$
   double T0Bond =  exp(-r0 * tau0);

   //Forward bond price at time $t$;
   double F0 = ForwadBondPV(r0, tau, tau0);

   // Forward bond variance at time $t$;
   double FBondVar = ForwardBondVariance(kappa, sigma, tau0, tau);
   // Forward bond standard deviation at time $t$
   double FBondSDev = sqrt(FBondVar);

   //double test = pow(sigma, 2)/(2*pow(kappa,3))*pow((1-exp(-kappa * (tau))),2)*(1-exp(-2*kappa*tau0));
   //double d2 = (log(F0/K)-(1.0/2.0)*test)/sqrt(test);
   //double d1 = (log(F0/K)+(1.0/2.0)*test)/sqrt(test);

   // calculate d_{-}=d2 in the Bond Call Option pricing formula
   double d2 = (log(F0/K)-(1.0/2.0)*FBondVar)/FBondSDev;

   // calculate d_{+}=d1 in the Bond Call Option pricing formula
   double d1 = (log(F0/K)+(1.0/2.0)*FBondVar)/FBondSDev;
   // or using d_{+}=d1=d2+FBondSDev, the relationship between $d1$ and $d2$
   // double d1 = d2 + FBondSDev;

   double call = TBond * normalDistribution(d1) - T0Bond * K * normalDistribution(d2);

   return call;
}

int main()
{
    double t = 0;   // calendar time $t$
    double T = 2; //Bond maturity; T-Bond
    double T0 = 1; // Option maturity
    double tau = T-t; //Time to bond maturity
    double tau0= T0-t; //Time to option maturity

    //r0 is the quoted bond yield on the calendar time $t$
    double r0 = 0.03;
    // interest rate process volatility $sigma$
    double sigma = 0.1;

    double kappa = 0.1;

    // strike price of T-Bond call Option at time $S$
    double K = 0.8;

    cout << "European call on T-Bond price with option maturity T0 is "
      << euroCall(K, r0, sigma, kappa, t, T, T0) << endl;

    return 0;
  }
