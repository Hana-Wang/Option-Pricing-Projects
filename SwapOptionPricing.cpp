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

double FBondVal(double F0, double FBonVar, double y)
{
  double FBondSDev = sqrt(FBondVar);
  double FBondVal = F0 * exp(y*FBondSDev - FBondVar/2);
  return FBondVal;
}


vector<double> FBondVarS(vector<double> &tau, double tau0, double sigma, double kappa)
{
  int n = tau.size();
  vector<double> FBondVars(n, 0.0);
  for(int i=0; i < n; i++)
  {
    // Forward bond variance at time $t$;
    FBondVars[i] = ForwardBondVariance(kappa, sigma, tau0, tau[i]);
  }

  return FBondVars;

}

// Calculate Forward bond present price vector at calendar time $t$; B(t;T0,Ti)
// assuming yield curve is flat, use r0
vector<double> FBondPVS(vector<double> &tau, double tau0, double r0)
{
  int n = tau.size();
  vector<double> FBondPVs(n, 0.0);
  for(int i=0; i < n; i++)
  {
    // Forward bond availabe value at calendar time $t$; B(t;T0,Ti)
    FBondPVs[i] = ForwardBondPV(r0, tau[i], tau0);
  }

  return FBondPVs;
}



// return forward bond price vector B(t;T0,Ti), for i=1,2,...,n
vector<double> FBondS(vector<double> &FBondPVs, vector<double> &FBondVars, double y)
{
    int n = FBondVars.size();
    vector<double> FBonds(n, 0.0);
    for(int i=0; i < n; i++)
      {
        //forward bond price B(T0;T0,Ti) at time T0 > t, $t$ is calendar time
        FBonds[i] = FBondVal(FBondPVs[i], FBonVars[i], y);
      }
    return FBonds;
}

// Calculate calculate Deltau[i]=T(i)-T(i-1)=tau(i)-tau(i-1)
// Here, tau(i)=T(i)-t; where $t$ is a calendar(quoted) time.
vector<double> DelTau(vector<double> &tau, double tau0)
{
  int n = tau.size();
  vector<double> Deltau(n, 0.0);

  // calculate Deltau[i]=T(i)-T(i-1)=tau(i)-tau(i-1)
  Deltau[0]=tau[0]-tau0;
  for(int i=1; i < n; i++)
  {
    Deltau[i] = tau[i]-tau[i-1];
  }

  return Deltau;
}

// using Taylor expantion approximation of exp(sqrt(FBondVar[i])*y-FBondVar[i]/2))
//approximatly equal to [1+sqrt(FBondVar[i])*y-FBondVar[i]/2]
// calculate initial guess y* using the derived formula
// KCoupR is fixed coupon rate to receive for the payer swap
double InitGuessY(vector<double> &FBondPVs, vector<double> &FBondVars, vector<double> &Deltau, double KCoupR)
{

  int n = FBondVars.size();

  // sum1 is the sum in the numerator of the calculation formula of y*
  double sum1 = 0;
  // sum2 is the sum in the denominator of the calculation formula of y*
  doulbe sum2 = 0;

  for(int i=0; i < n; i++)
  {
     sum1 += Deltau[i] * FBondPVs[i] * (1-FBondVars[i]/2);
     sum2 += Deltau[i] * FBondPVs[i] * sqrt(FBondVars[i]);
  }

  // calculate B(t;T0,TN)*(1-FBondVars(t;T0,TN)/2) in the y* formula numerator
  double bn = FBondPVs[n-1] * (1-FBondVars[n-1]/2);
  // calculate B(t;T0,TN) * sqrt(FBondVars(t;T0,TN)) in the y* formula denominator
  double dbn = FBondPVs[n-1] * sqrt(FBondVars[n-1]);

  //Calulate initial guess value y*
  doulbe y0 = (1-bn - KCoupR * sum1)/(KCoupR * sum2 + dbn);

  return y0;

}


//Newton-Raphson Method
//the function f(y)
double f(double y, vector<double> &Deltau, vector<double> &FBondPVs, vector<double> &FBondVars, double KCoupR)
{
  int n = Deltau.size();
  vector<double> FBonds = FBondS(FBondPVs, FBondVars, y);

  double Fsum = 0;
  for(int i=0; i < n; i++)
  {
     Fsum += Deltau[i] * FBonds[i];
  }

  doulbe a = FBonds[n-1] + KCoupR * Fsum - 1;

  return a;
}


//Newton-Raphson Method
//the derivative function fprime(y)
double fprime(double y, vector<double> &Deltau, vector<double> &FBondPVs, vector<double> &FBondVars, double KCoupR)
{
  int n = Deltau.size();
  vector<double> FBonds = FBondS(FBondPVs, FBondVars, y);

  double Fsum = 0;
  for(int i=0; i < n; i++)
  {
     Fsum += Deltau[i] * FBonds[i] * FBondVars[i];
  }

  doulbe b = FBondVars[i] * FBonds[n-1] + KCoupR * Fsum;

  return b;
}





//Newton-Raphson Method
//the function
double f(double y, double n, double kappa, double sigma)
{
  //consider a European swaption giving the holder have a right to
  // enter in $5$ years a swap into a $n$ years paying swap
  // consider the calendar time $t=0$
  double t = 0;
  double T0 = 5;
  double tau0 = T0 - t;

  // fixed coupon rate to receive for the payer swap
  double KCoupR = 0.062; // test number, assume the fixed coupon rate is 6.2%
  //vector<double> RI{ 0.03, 0.034, 0.037, 0.039, 0.04 };
  vector<double> tau(n, 0.0);
  vector<double> Deltau(n, 0.0);
  vector<double> Bonds(n, 0.0);
  vector<double> FBondPVs(n, 0.0);
  vector<double> FBonds(n, 0.0);
  vector<double> FBondVars(n, 0.0);


  for(int i=0; i < n; i++)
  {
    tau[i] = tau0 + i + 1;
    Bonds[i] = exp(-r0 * tau[i]);
    // Forward bond variance at time $t$;
    FBondVars[i] = ForwardBondVariance(kappa, sigma, tau0, tau[i]);
    FBondPVs[i] = ForwardBondPV(r0, tau[i], tau0);
    FBonds[i] = FBondVal(FBondPVs[i], FBonVars[i], y);
  }

  // calculate Deltau[i]=T(i)-T(i-1)
  Deltau[0]=tau[0]-tau0;
  for(int i=1; i < n; i++)
  {
    Deltau[i] = tau[i]-tau[i-1];
  }


  // sum1 is the sum in the numerator of the calculation formula of y*
  double sum1 = 0;
  // sum2 is the sum in the denominator of the calculation formula of y*
  doulbe sum2 = 0;

  for(int i=0; i < n; i++)
  {
     sum1 += Deltau[i] * Bonds[i] * (1-FBondVars[i]/2);
     sum2 += Deltau[i] * Bonds[i] * sqrt(FBondVars[i]);
  }

  // calculate B(t;T0,TN)*(1-FBondVars(t;T0,TN)) in the y* formula numerator
  double bn = Bonds[n-1] * (1-FBondVars[n]/2);
  double dbn = Bonds[n-1] * sqrt(FBondVars[i]);
  //Calulate initial guess value y*
  doulbe y0 = (1-bn - KCoupR * sum1)/(sum2 + dbn);

  double Fsum = 0;
  for(int i=0; i < n; i++)
  {
     Fsum += Deltau[i] * FBonds[i];
  }

  doulbe a = FBonds[n-1] + KCoupR * Fsum - 1 = 0;

  return a;
}
