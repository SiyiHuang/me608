#include<stdlib.h>

double * tdmaSol(double * ap, double * ae, double * aw, double * b, double * X, int N)
{
  int i;
  double * AP[N];
  for(i = 0; i<N; i++) AP[i] = ap[i];
  double * AW[N-1];
  for(i = 0; i<N-1; i++) AW[i] = aw[i];
  double * AE[N-1];
  for(i = 0; i<N-1; i++) AE[i] = ae[i];
  double * B[N-1];
  for(i = 0; i<N; i++) B[i] = b[i];

  for(i = 1, i < N, i++)
    {
      double r = AW[i] / AP[i-1];
      AP[i] = AP [i] - r * AE[i-1];
      B[i] = B[i] - r * B[i-1];
    }
  X[N] = B[N] / AP[N];
  for(i = N-2; i > -1; i--)
    {
      X[i] = (B[i] - AW[i] * X[i+1]) / AP[i];
    }
  return X;
}


double * 
