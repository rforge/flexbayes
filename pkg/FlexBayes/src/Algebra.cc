#ifdef WIN32
#pragma warning(disable:4290)
#endif

#include <stdio.h>
#include <math.h>
#include "R.h"
#include "Algebra.h"


/* These comments were added by A. Murua (09/13/2002)

   It seems that this Cholesky decomposition takes a full
   matrix a of dimensions nxn as input.
   It returns the Cholesky decomposition in the lower triangular
   part of a, except for the diagonal which is returned in the
   vector p.

   Algorithm looks similar to that one in Num. Recipes in C, except
   that here several checks for numerical errors are made.
*/
void Algebra::choldc(double *a, long n, double *p, double mxff)  
throw(rtErr) 
{
    int i, j, k;
    int ind;
    double sum, minl, minl2, maxadd, minljj;
    const double EPS = 1.11022e-16;

    ind = 0;
    minl = pow( EPS, 0.25 ) * mxff;
    if ( mxff == 0.0 ) 
    {
      sum = 0.0;
      for ( i = 0; i < n; i++ ) 
      {
	if ( fabs(a[i*n + i]) > sum )
	{
	  sum = fabs(a[i*n + i]);
        }
      }
      mxff = sqrt(sum);
      minl2 = sqrt(EPS) * mxff;
    }//end if

    maxadd = 0.0;
    for (j = 0; j < n; j++) 
    {
      p[j] = a[j*n + j];
      for(i = j - 1; i >= 0; i--)
	p[j] -= a[j*n + i]*a[j*n + i];

      minljj = 0.0;
      for(i = j + 1; i < n; i++) 
      {
	a[i*n + j] = a[j*n + i];
	for(k = j - 1; k >= 0; k--) 
	  a[i*n + j] -= a[i*n + k]*a[j*n + k];
	if(fabs(a[i*n + j]) > minljj)
	  minljj = fabs(a[i*n + j]);
      }

      minljj = minljj/mxff;
      if(minljj < minl)
	minljj = minl;
      
      if(p[j] > minljj*minljj)
	p[j] = sqrt(p[j]);
      else {
	if(minljj < minl2)
	  minljj = minl2;
	if(maxadd < minljj*minljj - p[j])
	  maxadd = minljj*minljj - p[j];
	p[j] = minljj;
	ind++;
      }
     
      for(i = j + 1; i < n; i++) 
		a[i*n + j] /= p[j];

    }


    if ( ind ) 
    {
      //printf("algebra:choldc: ind = %d\n", ind ); fflush(stdout);

      //printf( "Algebra::choldc: Augmented Cholesky decomposition Error.\n\n" );
      //char the_error[] = "Algebra::choldc: Augmented Cholesky decomposition Error.";
      //rtErr runtime_error(the_error);
      //throw runtime_error;
      warning("Algebra::choldc: Augmented cholesky decomposition problem (not positive definite?): ind=%d", ind) ;
    }

}



bool Algebra::isCholdcOK( double *a, long n, double *p, double mxff )  
{
    bool is_chol_OK;
    int i, j, k;
    int ind;
    double sum, minl, minl2, maxadd, minljj;
    const double EPS = 1.11022e-16;

    ind = 0;
    minl = pow( EPS, 0.25 ) * mxff;
    if ( mxff == 0.0 ) 
    {
      sum = 0.0;
      for ( i = 0; i < n; i++ ) 
      {
	if ( fabs(a[i*n + i]) > sum )
	{
	  sum = fabs(a[i*n + i]);
        }
      }
      mxff = sqrt(sum);
      minl2 = sqrt(EPS) * mxff;
    }//end if

    maxadd = 0.0;
    for (j = 0; j < n; j++) 
    {
      p[j] = a[j*n + j];
      for(i = j - 1; i >= 0; i--)
	p[j] -= a[j*n + i]*a[j*n + i];

      minljj = 0.0;
      for(i = j + 1; i < n; i++) 
      {
	a[i*n + j] = a[j*n + i];
	for(k = j - 1; k >= 0; k--) 
	  a[i*n + j] -= a[i*n + k]*a[j*n + k];
	if(fabs(a[i*n + j]) > minljj)
	  minljj = fabs(a[i*n + j]);
      }

      minljj = minljj/mxff;
      if(minljj < minl)
	minljj = minl;
      
      if(p[j] > minljj*minljj)
	p[j] = sqrt(p[j]);
      else {
	if(minljj < minl2)
	  minljj = minl2;
	if(maxadd < minljj*minljj - p[j])
	  maxadd = minljj*minljj - p[j];
	p[j] = minljj;
	ind++;
      }
     
      for(i = j + 1; i < n; i++) 
		a[i*n + j] /= p[j];

    }


    if ( ind ) 
    {
      is_chol_OK = false;
      //printf("algebra:isCholOK: ind = %d\n", ind ); fflush(stdout);
    }
    else
    {
      is_chol_OK = true;
    }

    return ( is_chol_OK );

}//end



/*
   solves a linear equation: A x = b
   on input the lower triangle of a + p corresponds to
   the Cholesky decomposition of A. 
   b is the right-hand-side vector.
   The solution is returned in x.
*/

void Algebra::cholsl(double *a, long n, double *p, double *b, double *x) {

    double sum;
		long i = 0;
		long j = 0;
    
    for(i = 0; i < n; i++) {
      sum = b[i];
      for(j = i-1; j >= 0; j--)
        sum -= a[i*n + j]*x[j];
      x[i] = sum/p[i];
    }

    for(i = n-1; i >= 0; i--) {
      sum = x[i];
      for(j = i+1; j < n; j++)
        sum -= a[j*n + i]*x[j];
      x[i] = sum/p[i];
    }

}

void Algebra::ludcmp(double *a, int n, int *indx, double *d, double *vv)
  throw(rtErr) 
{

  int i, imax, j, k;
  const double TINY = 1.0e-20;
  double big, dum, sum, temp;

  *d = 1.0;
  for(i = 0; i < n; i++) {
    big = 0.0;
    for(j = 0; j < n; j++) {
      temp = fabs(a[i*n + j]);
      if(temp > big) 
	big = temp;
    }
    if(big == 0.0) {
			char the_error[] = "Singular matrix in LU decomposition.";
      rtErr runtime_error(the_error);
      throw runtime_error;
    }
    vv[i] = 1.0/big;
  }
  for(j = 0; j < n; j++) {
    for(i = 0; i < j; i++) {
      sum = a[i*n + j];
      for(k = 0; k < i; k++) {
	sum -= a[i*n + k]*a[k*n + j];
      }
      a[i*n + j] = sum;
    }
    big = 0.0;
    for(i = j; i < n; i++) {
      sum = a[i*n + j];
      for(k = 0; k < j; k++)
	sum -= a[i*n + k]*a[k*n + j];
      a[i*n + j] = sum;
      dum = vv[i]*fabs(sum);
      if(dum >= big) {
	big = dum;
	imax = i;
      }
    }
    if(j != imax) {
      for(k = 0; k < n; k++) {
	dum = a[imax*n + k];
	a[imax*n + k] = a[j*n + k];
	a[j*n + k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(a[j*n +j] == 0.0)
      a[j*n + j] = TINY;
    if(j != n-1) {
      dum = 1.0/a[j*n + j];
      for(i = j+1; i < n; i++)
	a[i*n + j] *= dum;
    }
  }

}

void Algebra::lubksb(double *a, int n, int *indx, double *b) {
  
  int i, ii = -1, ip, j;
  double sum;

  for(i = 0; i < n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if(ii != -1) {
      for(j = ii; j < i; j++)
	sum -= a[i*n + j]*b[j];
    } else if( sum != 0.0) {
      ii = i;
    }
    b[i] = sum;
  }
  for(i = n-1; i >= 0; i--) {
    sum = b[i];
    for(j = i+1; j < n; j++)
      sum -= a[i*n + j]*b[j];
    b[i] = sum/a[i*n + i];
  }

}

void Algebra::qrdcmp(double *a, int n, double *c, double *d, 
		     int *sing, double eps) {

  int i, j, k;
  double scale, sigma, sum, tau;

  *sing = 0;
  for(k = 0; k < n-1; k++) {
    scale = 0.0;
    for(i = k; i < n; i++) {
      if(fabs(a[i*n + k]) > scale)
	scale = fabs(a[i*n + k]);
    }
    if(scale <= eps){
      *sing = 1;
      c[k] = 0.0;
      d[k] = 0.0;
    } else {
      for(i = k; i < n; i++)
	a[i*n + k] /=scale;
      sum =0.0;
      for(i = k; i < n; i++)
	sum += a[i*n + k]*a[i*n + k];
      sigma = sqrt(sum);
      if(a[k*n + k] < 0.0)
	sigma = -sigma;
      a[k*n + k] += sigma;
      c[k] = sigma*a[k*n + k];
      d[k] = -scale*sigma;
      for(j = k+1; j < n; j++) {
	sum = 0.0;
	for(i = k; i < n; i++)
	  sum += a[i*n + k]*a[i*n + j];
	tau = sum/c[k];
	for(i = k; i < n; i++)
	  a[i*n + j] -= tau*a[i*n + k];
      }
    }
  }
  d[n-1] = a[(n-1)*n + n-1];
  if(fabs(d[n-1]) < eps) *sing = 1;

}

void Algebra::qrsolv(double *a, int n, double *c, double *d, double *b) {

  int i, j;
  double sum, tau;

  for(j = 0; j < n-1; j++) {
    sum = 0.0;
    for(i = j; i < n; i++)
      sum += a[i*n + j]*b[i];
    tau = sum/c[j];
    for(i = j; i < n; i++)
      b[i] -= tau*a[i*n + j];
  }
  rsolv(a, n, d, b);

}

void Algebra::rsolv(double *a, int n, double *d, double *b) {

  int i, j;
  double sum;

  b[n-1] /= d[n-1];
  for(i = n-2; i >= 0; i--) {
    sum = 0.0;
    for(j = i+1; j < n; j++)
      sum += a[i*n + j]*b[j];
    b[i] = (b[i] - sum)/d[i];
  }

}
