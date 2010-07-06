
#include <complex.h>
#include <string.h>
#include "BLAS.h"
#include "vectorOps.h"

static double complex
complex_mul (double complex x, double complex y)
{
	double xr = creal(x), xi = cimag(x);
	double yr = creal(y), yi = cimag(y);
	double zr = xr * yr - xi * yi;
	double zi = xr * yi + xi * yr;
        double complex z = zr + I * zi;	
	return z;
}


static void
zVectorClear (int n, double complex *z)
{
        memset(z, 0, n * 2 * sizeof(double));
}

static void
zVectorCopy (int n, const double complex *x, double complex *z)
{
        if (x != z) {
                memcpy(z, x, n * 2 * sizeof(double));
        }
}

void
zVectorConj (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = conj(x[i]);
        }
}

void
zVectorScale (int n, const double complex *palpha, const double complex *x, double complex *z)
{
        double complex alpha = *palpha;
        
        if (alpha == 1) {
                zVectorCopy(n, x, z);
        } else if (alpha == 0) {
                zVectorClear(n, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = complex_mul(alpha, x[i]);
                }
        }
        
        /* Note: Using zscal sometimes gives differences in the least
         * significant bit of the answer */
}

void
zVectorShift (int n, const double complex *palpha, const double complex *x, double complex *z)
{
        double complex alpha = *palpha;
        
        if (alpha == 0) {
                zVectorCopy(n, x, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha + x[i];
                }
        }
}

void
zVectorNeg (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = -x[i];
        }
}

void
zVectorAbs (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cabs(x[i]) + 0*I;
        }
}

static void
zsgn (const double complex x, double complex *z)
{
        if (x == 0) {
                *z = x;
        } else {
                double r = cabs(x);
                *z = creal(x)/r + cimag(x)/r*I;
        }
}

void
zVectorSgn (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                zsgn(x[i], &(z[i]));
        }
}

void
zVectorInv (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = 1.0 / x[i];
        }
}

static void
zVectorAxpy (int n, const double complex *palpha, const double complex *x, const double complex *y, double complex *z)
{
        double complex alpha = *palpha;
        
        if (alpha == 0) {
                zVectorCopy(n, y, z);
        } else if (y == z) {
                blas_zaxpy(n, palpha, x, 1, z, 1);
        } else if (alpha == 1 && x == z) {
                blas_zaxpy(n, palpha, y, 1, z, 1);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i] + y[i];
                }
        }
}

void
zVectorAxpby (int n, const double complex *palpha, const double complex *x, const double complex *pbeta, const double complex *y, double complex *z)
{
        double complex alpha = *palpha;
        double complex beta = *pbeta;        
                
        if (alpha == 0) {
                zVectorScale(n, pbeta, y, z);
        } else if (alpha == 1) {
                zVectorAxpy(n, pbeta, y, x, z);
        } else if (beta == 0) {
                zVectorScale(n, palpha, x, z);
        } else if (beta == 1) {
                zVectorAxpy(n, palpha, x, y, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i] + beta * y[i];
                }
        }
}

void
zVectorAdd (int n, const double complex *x, const double complex *y, double complex *z)
{
        dVectorAdd(2 * n, (const double *)x, (const double *)y, (double *)z);
}

void zVectorSub (int n, const double complex *x, const double complex *y, double complex *z)
{
        dVectorSub(2 * n, (const double *)x, (const double *)y, (double *)z);
}

void
zVectorMul (int n, const double complex *x, const double complex *y, double complex *z)
{
        /* Note: using ztbmv sometimes gives different answers in the lower
         * bits of precision */
        /*
        if (y == z) {
                blas_ztbmv(BlasUpper, BlasNoTrans, BlasNonUnit, n, 0, x, 1,
                           z, 1);
        } else if (x == z) {
                blas_ztbmv(BlasUpper, BlasNoTrans, BlasNonUnit, n, 0, y, 1,
                           z, 1);
        } */
        int i;
        for (i = 0; i < n; i++) {
                z[i] = complex_mul(x[i], y[i]);
        }
}

void
zVectorDiv (int n, const double complex *x, const double complex *y, double complex *z)
{
        /* Note: using dtbsv gives slightly different answers
         */
        /*
        if (y == z) {
                blas_ztbsv(BlasUpper, BlasNoTrans, BlasNonUnit, n, 0, x, 1,
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
zVectorExp (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cexp(x[i]);
        }
}

void
zVectorSqrt (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = csqrt(x[i]);
        }
}

void
zVectorLog (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = clog(x[i]);
        }        
}

void
zVectorPow (int n, const double complex *x, const double complex *y, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cpow(x[i], y[i]);
        }        
}

void
zVectorSin (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = csin(x[i]);
        }                
}

void
zVectorCos (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = ccos(x[i]);
        }                        
}

void
zVectorTan (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = ctan(x[i]);
        }                                
}

void
zVectorASin (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = casin(x[i]);
        }                                        
}

void
zVectorACos (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cacos(x[i]);
        }                                                
}

void
zVectorATan (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = catan(x[i]);
        }
}

void
zVectorSinh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = csinh(x[i]);
        }                
}

void
zVectorCosh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = ccosh(x[i]);
        }                        
}

void
zVectorTanh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = ctanh(x[i]);
        }                                
}

void
zVectorASinh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = casinh(x[i]);
        }                                        
}

void
zVectorACosh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cacosh(x[i]);
        }                                                
}

void
zVectorATanh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = catanh(x[i]);
        }
}

