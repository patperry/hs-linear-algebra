#ifndef _VECTOR_OPS_H
#define _VECTOR_OPS_H

#include <complex.h>

void dVectorScale (int n, double alpha, const double *x, double *z);
void dVectorShift (int n, double alpha, const double *x, double *z);
void dVectorNeg (int n, const double *x, double *z);
void dVectorAbs (int n, const double *x, double *z);
void dVectorSgn (int n, const double *x, double *z);
void dVectorInv (int n, const double *x, double *z);

void dVectorAxpby (int n, double alpha, const double *x, double beta, const double *y, double *z);
void dVectorAdd (int n, const double *x, const double *y, double *z);
void dVectorSub (int n, const double *x, const double *y, double *z);
void dVectorMul (int n, const double *x, const double *y, double *z);
void dVectorDiv (int n, const double *x, const double *y, double *z);

void dVectorExp (int n, const double *x, double *z);
void dVectorSqrt (int n, const double *x, double *z);
void dVectorLog (int n, const double *x, double *z);
void dVectorPow (int n, const double *x, const double *y, double *z);

void dVectorSin (int n, const double *x, double *z);
void dVectorCos (int n, const double *x, double *z);
void dVectorTan (int n, const double *x, double *z);
void dVectorASin (int n, const double *x, double *z);
void dVectorACos (int n, const double *x, double *z);
void dVectorATan (int n, const double *x, double *z);

void dVectorSinh (int n, const double *x, double *z);
void dVectorCosh (int n, const double *x, double *z);
void dVectorTanh (int n, const double *x, double *z);
void dVectorASinh (int n, const double *x, double *z);
void dVectorACosh (int n, const double *x, double *z);
void dVectorATanh (int n, const double *x, double *z);

void zVectorConj (int n, const double complex *x, double complex *z);
void zVectorScale (int n, const double complex *alpha, const double complex *x, double complex *z);
void zVectorShift (int n, const double complex *alpha, const double complex *x, double complex *z);
void zVectorNeg (int n, const double complex *x, double complex *z);
void zVectorAbs (int n, const double complex *x, double complex *z);
void zVectorSgn (int n, const double complex *x, double complex *z);
void zVectorInv (int n, const double complex *x, double complex *z);

void zVectorAxpby (int n, const double complex *alpha, const double complex *x, const double complex *beta, const double complex *y, double complex *z);
void zVectorAdd (int n, const double complex *x, const double complex *y, double complex *z);
void zVectorSub (int n, const double complex *x, const double complex *y, double complex *z);
void zVectorMul (int n, const double complex *x, const double complex *y, double complex *z);
void zVectorDiv (int n, const double complex *x, const double complex *y, double complex *z);

void zVectorExp (int n, const double complex *x, double complex *z);
void zVectorSqrt (int n, const double complex *x, double complex *z);
void zVectorLog (int n, const double complex *x, double complex *z);
void zVectorPow (int n, const double complex *x, const double complex *y, double complex *z);

void zVectorSin (int n, const double complex *x, double complex *z);
void zVectorCos (int n, const double complex *x, double complex *z);
void zVectorTan (int n, const double complex *x, double complex *z);
void zVectorASin (int n, const double complex *x, double complex *z);
void zVectorACos (int n, const double complex *x, double complex *z);
void zVectorATan (int n, const double complex *x, double complex *z);

void zVectorSinh (int n, const double complex *x, double complex *z);
void zVectorCosh (int n, const double complex *x, double complex *z);
void zVectorTanh (int n, const double complex *x, double complex *z);
void zVectorASinh (int n, const double complex *x, double complex *z);
void zVectorACosh (int n, const double complex *x, double complex *z);
void zVectorATanh (int n, const double complex *x, double complex *z);

#endif /* _VECTOR_OPS_H */
