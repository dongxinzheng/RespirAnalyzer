#include <Rcpp.h>
using namespace Rcpp;

int polyfit(const double* const dependentValues,
            const double* const independentValues,
            unsigned int        countOfElements,
            unsigned int        order,
            double*             coefficients);
List polyfitY(NumericVector x, NumericVector y, int n);
NumericVector polyfitCoef(NumericVector x, NumericVector y, int n);
// [[Rcpp::export]]
List cppMFDFA(const NumericVector tsx, const IntegerVector scale, const int m, const NumericVector q) 
{
  
  NumericVector Hq(q.size());
  NumericMatrix Fqi(scale.size(), q.size());
  NumericMatrix RegLine(scale.size(), q.size());
  
  
  NumericVector X = cumsum(tsx-mean(tsx));
  for(int i=0; i<scale.size(); i++)
  {
    int seg = X.size()/scale[i];
    NumericVector x(scale[i]);
    NumericVector y(scale[i]);
    NumericVector Y(scale[i]);
    NumericVector rmvi(seg);
    
    for(int vi=0; vi<seg; vi++)
    {
      x = Range(vi*scale[i], (vi+1)*scale[i]-1);
      y = X[x];
      Y = polyfitY(x,y,m)["Y"];
      rmvi[vi] = sqrt(mean(pow(Y-y, 2)));
    }
    
    for(int nq =0; nq<q.size(); nq++)
    {
      if(q[nq]==0)
        Fqi(i,nq) = exp(0.5*mean(log(pow(rmvi, 2))));
      else
        Fqi(i,nq) = pow(mean(pow(rmvi, q[nq])),1.0/q[nq]);
      
    }
    
  }
  
  double log2 = log(2);
  
  for(int nq=0; nq<q.size(); nq++)
  {
    List pf = polyfitY(log(scale),log((Fqi(_,nq))),1);
    Hq[nq] = as<NumericVector>(pf["coef"])[1];
    RegLine(_, nq) = as<NumericVector>(pf["Y"])/log2;
  }
  
  NumericVector tq(q.size());
  NumericVector hq(q.size()-1);
  NumericVector Dq(q.size()-1);
  tq =  Hq*q - 1;
  hq = diff(tq)/(q[2]-q[1]);
  Dq = head(q,q.size()-1)*hq-head(tq, tq.size()-1);

  return List::create(Named("Hq")=Hq, Named("tau_q")=tq, Named("hq")=hq, Named("Dq")=Dq,
                            Named("Fqi")=Fqi, Named("line")=RegLine);
}


//## intern functions: polyfit ####
NumericVector polyfitCoef(NumericVector x, NumericVector y, int n)
{
  NumericVector coef(n+1);
  polyfit(x.begin(), y.begin(), x.size(), n, coef.begin());
  return round(coef,4);
}

List polyfitY(NumericVector x, NumericVector y, int n)
{
  NumericVector coef(n+1);
  NumericVector predictedY(x.size());
  polyfit(x.begin(), y.begin(), x.size(), n, coef.begin());
  for(int i=0; i<x.size(); i++)
  {
    for(int j=0; j<coef.size(); j++)
    {
      predictedY[i] += coef[j]*pow(x[i], j);
    }
  }
  return List::create(Named("coef")=round(coef,4), Named("Y")=predictedY);
}

//----------------------------------------------------
//
// METHOD:  polyfit
//
// INPUTS:  dependentValues[0..(countOfElements-1)]
//          independentValues[0...(countOfElements-1)]
//          countOfElements
//          order - Order of the polynomial fitting
//
// OUTPUTS: coefficients[0..order] - indexed by term
//               (the (coef*x^3) is coefficients[3])
//
//----------------------------------------------------
int polyfit(const double* const dependentValues,
            const double* const independentValues,
            unsigned int        countOfElements,
            unsigned int        order,
            double*             coefficients)
{
  // Declarations...
  // ----------------------------------
  enum {maxOrder = 5};
  
  double B[maxOrder+1] = {0.0f};
  double P[((maxOrder+1) * 2)+1] = {0.0f};
  double A[(maxOrder + 1)*2*(maxOrder + 1)] = {0.0f};
  
  double x, y, powx;
  
  unsigned int ii, jj, kk;
  
  // Verify initial conditions....
  // ----------------------------------
  
  // This method requires that the countOfElements > 
  // (order+1) 
  if (countOfElements <= order)
    return -1;
  
  // This method has imposed an arbitrary bound of
  // order <= maxOrder.  Increase maxOrder if necessary.
  if (order > maxOrder)
    return -1;
  
  // Begin Code...
  // ----------------------------------
  
  // Identify the column vector
  for (ii = 0; ii < countOfElements; ii++)
  {
    x    = dependentValues[ii];
    y    = independentValues[ii];
    powx = 1;
    
    for (jj = 0; jj < (order + 1); jj++)
    {
      B[jj] = B[jj] + (y * powx);
      powx  = powx * x;
    }
  }
  
  // Initialize the PowX array
  P[0] = countOfElements;
  
  // Compute the sum of the Powers of X
  for (ii = 0; ii < countOfElements; ii++)
  {
    x    = dependentValues[ii];
    powx = dependentValues[ii];
    
    for (jj = 1; jj < ((2 * (order + 1)) + 1); jj++)
    {
      P[jj] = P[jj] + powx;
      powx  = powx * x;
    }
  }
  
  // Initialize the reduction matrix
  //
  for (ii = 0; ii < (order + 1); ii++)
  {
    for (jj = 0; jj < (order + 1); jj++)
    {
      A[(ii * (2 * (order + 1))) + jj] = P[ii+jj];
    }
    
    A[(ii*(2 * (order + 1))) + (ii + (order + 1))] = 1;
  }
  
  // Move the Identity matrix portion of the redux matrix
  // to the left side (find the inverse of the left side
  // of the redux matrix
  for (ii = 0; ii < (order + 1); ii++)
  {
    x = A[(ii * (2 * (order + 1))) + ii];
    if (x != 0)
    {
      for (kk = 0; kk < (2 * (order + 1)); kk++)
      {
        A[(ii * (2 * (order + 1))) + kk] = 
          A[(ii * (2 * (order + 1))) + kk] / x;
      }
      
      for (jj = 0; jj < (order + 1); jj++)
      {
        if ((jj - ii) != 0)
        {
          y = A[(jj * (2 * (order + 1))) + ii];
          for (kk = 0; kk < (2 * (order + 1)); kk++)
          {
            A[(jj * (2 * (order + 1))) + kk] = 
              A[(jj * (2 * (order + 1))) + kk] -
              y * A[(ii * (2 * (order + 1))) + kk];
          }
        }
      }
    }
    else
    {
      // Cannot work with singular matrices
      return -1;
    }
  }
  
  // Calculate and Identify the coefficients
  for (ii = 0; ii < (order + 1); ii++)
  {
    for (jj = 0; jj < (order + 1); jj++)
    {
      x = 0;
      for (kk = 0; kk < (order + 1); kk++)
      {
        x = x + (A[(ii * (2 * (order + 1))) + (kk + (order + 1))] *
          B[kk]);
      }
      coefficients[ii] = x;
    }
  }
  
  return 0;
}
