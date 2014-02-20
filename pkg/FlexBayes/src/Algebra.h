
#ifndef BAYES_ALGEBRA_H
#define BAYES_ALGEBRA_H

#include "rtErr.h"

class Algebra {

public:

  /* Cholesky decomposition  */
  static void choldc(double *a, int n, double *p, double mxff) 
	  throw(rtErr);
  static bool isCholdcOK( double *a, int n, double *p, double mxff );
  static void cholsl(double *a, int n, double *p, double *b, double *x);

  /* LU decomposition */
  static void ludcmp(double *a, int n, int *indx, double *d, double *vv)
    throw(rtErr);
  static void lubksb(double *a, int n, int *indx, double *b);

  /* QR decomposition */
  static void qrdcmp(double *a, int n, double *c, double *d, 
		             int *sing, double eps);
  static void qrsolv(double *a, int n, double *c, double *d, double *b);
  static void rsolv(double *a, int n, double *d, double *b);

};


#endif /* BAYES_ALGEBRA_H */
