#ifndef _VMATH_H
#define _VMATH_H

#include <complex.h>

void vdScale (int n, double alpha, const double *x, double *z);
void vdShift (int n, double alpha, const double *x, double *z);
void vdNeg (int n, const double *x, double *z);
void vdAbs (int n, const double *x, double *z);
void vdSgn (int n, const double *x, double *z);
void vdInv (int n, const double *x, double *z);

void vdAxpby (int n, double alpha, const double *x, double beta, const double *y, double *z);
void vdAdd (int n, const double *x, const double *y, double *z);
void vdSub (int n, const double *x, const double *y, double *z);
void vdMul (int n, const double *x, const double *y, double *z);
void vdDiv (int n, const double *x, const double *y, double *z);

void vdExp (int n, const double *x, double *z);
void vdSqrt (int n, const double *x, double *z);
void vdLog (int n, const double *x, double *z);
void vdPow (int n, const double *x, const double *y, double *z);

void vdSin (int n, const double *x, double *z);
void vdCos (int n, const double *x, double *z);
void vdTan (int n, const double *x, double *z);
void vdASin (int n, const double *x, double *z);
void vdACos (int n, const double *x, double *z);
void vdATan (int n, const double *x, double *z);

void vdSinh (int n, const double *x, double *z);
void vdCosh (int n, const double *x, double *z);
void vdTanh (int n, const double *x, double *z);
void vdASinh (int n, const double *x, double *z);
void vdACosh (int n, const double *x, double *z);
void vdATanh (int n, const double *x, double *z);

void vzConj (int n, const double complex *x, double complex *z);
void vzScale (int n, const double complex *alpha, const double complex *x, double complex *z);
void vzShift (int n, const double complex *alpha, const double complex *x, double complex *z);
void vzNeg (int n, const double complex *x, double complex *z);
void vzAbs (int n, const double complex *x, double complex *z);
void vzSgn (int n, const double complex *x, double complex *z);
void vzInv (int n, const double complex *x, double complex *z);

void vzAxpby (int n, const double complex *alpha, const double complex *x, const double complex *beta, const double complex *y, double complex *z);
void vzAdd (int n, const double complex *x, const double complex *y, double complex *z);
void vzSub (int n, const double complex *x, const double complex *y, double complex *z);
void vzMul (int n, const double complex *x, const double complex *y, double complex *z);
void vzDiv (int n, const double complex *x, const double complex *y, double complex *z);

void vzExp (int n, const double complex *x, double complex *z);
void vzSqrt (int n, const double complex *x, double complex *z);
void vzLog (int n, const double complex *x, double complex *z);
void vzPow (int n, const double complex *x, const double complex *y, double complex *z);

void vzSin (int n, const double complex *x, double complex *z);
void vzCos (int n, const double complex *x, double complex *z);
void vzTan (int n, const double complex *x, double complex *z);
void vzASin (int n, const double complex *x, double complex *z);
void vzACos (int n, const double complex *x, double complex *z);
void vzATan (int n, const double complex *x, double complex *z);

void vzSinh (int n, const double complex *x, double complex *z);
void vzCosh (int n, const double complex *x, double complex *z);
void vzTanh (int n, const double complex *x, double complex *z);
void vzASinh (int n, const double complex *x, double complex *z);
void vzACosh (int n, const double complex *x, double complex *z);
void vzATanh (int n, const double complex *x, double complex *z);

#endif /* _VMATH_H */
