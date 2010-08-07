
#include "config.h"
#include <math.h>
#include <string.h>
#include "vmath.h"

#define la_int int

extern void
F77_FUNC(daxpy) (const la_int *n, const double *da, const double *dx,
	         const la_int *incx, double *dy, const la_int *incy);



static void
vdClear (int n, double *z)
{
        memset(z, 0, n * sizeof(double));
}

static void
vdCopy (int n, const double *x, double *z)
{
        if (x != z) {
                memcpy(z, x, n * sizeof(double));
        }
}

void
vdScale (int n, double alpha, const double *x, double *z)
{
        if (alpha == 1) {
                vdCopy(n, x, z);
        } else if (alpha == 0) {
                vdClear(n, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i];
                }
        }
        
        /* Note: Using dscal sometimes gives differences in the least
         * significant bit of the answer */
}

void
vdShift (int n, double alpha, const double *x, double *z)
{
        if (alpha == 0) {
                vdCopy(n, x, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha + x[i];
                }
        }
}

void
vdNeg (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = -x[i];
        }
}

void
vdAbs (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = fabs(x[i]);
        }
}

static double
dsgn (double x)
{
        if (x == 0 || isnan(x)) {
                return x;
        } else {
                return copysign(1, x);
        }
}

void
vdSgn (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = dsgn(x[i]);
        }
}

void
vdInv (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = 1.0 / x[i];
        }
}

static void
vdAxpy (int n, double alpha, const double *x, const double *y, double *z)
{
        la_int one = 1;
        double done = 1.0;
        
        if (alpha == 0) {
                vdCopy(n, y, z);
        } else if (y == z) {
                F77_FUNC(daxpy) (&n, &alpha, x, &one, z, &one);
        } else if (alpha == 1 && x == z) {
                F77_FUNC(daxpy) (&n, &done, y, &one, z, &one);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i] + y[i];
                }
        }
}

void
vdAxpby (int n, double alpha, const double *x, double beta, const double *y, double *z)
{
        if (alpha == 0) {
                vdScale(n, beta, y, z);
        } else if (alpha == 1) {
                vdAxpy(n, beta, y, x, z);
        } else if (beta == 0) {
                vdScale(n, alpha, x, z);
        } else if (beta == 1) {
                vdAxpy(n, alpha, x, y, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i] + beta * y[i];
                }
        }
}

void
vdAdd (int n, const double *x, const double *y, double *z)
{
        vdAxpy(n, 1, x, y, z);
}

void vdSub (int n, const double *x, const double *y, double *z)
{
        vdAxpy(n, -1, y, x, z);
}

void
vdMul (int n, const double *x, const double *y, double *z)
{
        /* Note: using dtbmv sometimes gives different answers in the lower
         * bits of precision */
        
        /*if (y == z) {
                blas_dtbmv(BlasUpper, BlasNoTrans, BlasNonUnit, n, 0, x, 1,
                           z, 1);
        } else if (x == z) {
                blas_dtbmv(BlasUpper, BlasNoTrans, BlasNonUnit, n, 0, y, 1,
                           z, 1);
        } else { */
                
        int i;
        for (i = 0; i < n; i++) {
                z[i] = x[i] * y[i];
        }
}

void
vdDiv (int n, const double *x, const double *y, double *z)
{
        /* Note: using dtbsv gives slightly different answers
         */
        /*if (y == z) {
                blas_dtbsv(BlasUpper, BlasNoTrans, BlasNonUnit, n, 0, x, 1,
                           z, 1);
        } else */
        
        if (x == z) {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = z[i] / y[i];
                }
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = x[i] / y[i];
                }
        }
}

void
vdExp (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = exp(x[i]);
        }
}

void
vdSqrt (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = sqrt(x[i]);
        }
}

void
vdLog (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = log(x[i]);
        }        
}

void
vdPow (int n, const double *x, const double *y, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = pow(x[i], y[i]);
        }        
}

void
vdSin (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = sin(x[i]);
        }                
}

void
vdCos (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cos(x[i]);
        }                        
}

void
vdTan (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = tan(x[i]);
        }                                
}

void
vdASin (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = asin(x[i]);
        }                                        
}

void
vdACos (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = acos(x[i]);
        }                                                
}

void
vdATan (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = atan(x[i]);
        }
}

void
vdSinh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = sinh(x[i]);
        }                
}

void
vdCosh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cosh(x[i]);
        }                        
}

void
vdTanh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = tanh(x[i]);
        }                                
}

void
vdASinh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = asinh(x[i]);
        }                                        
}

void
vdACosh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = acosh(x[i]);
        }                                                
}

void
vdATanh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = atanh(x[i]);
        }
}

