
#include "config.h"
#include <complex.h>
#include <string.h>
#include "vmath.h"

#define la_int int

extern void
F77_FUNC(zaxpy) (const la_int *n, const void *za, const void *zx,
	         const la_int *incx, void *zy, const la_int *incy);


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
vzClear (int n, double complex *z)
{
        memset(z, 0, n * 2 * sizeof(double));
}

static void
vzCopy (int n, const double complex *x, double complex *z)
{
        if (x != z) {
                memcpy(z, x, n * 2 * sizeof(double));
        }
}

void
vzConj (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = conj(x[i]);
        }
}

void
vzScale (int n, const double complex *palpha, const double complex *x, double complex *z)
{
        double complex alpha = *palpha;
        
        if (alpha == 1) {
                vzCopy(n, x, z);
        } else if (alpha == 0) {
                vzClear(n, z);
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
vzShift (int n, const double complex *palpha, const double complex *x, double complex *z)
{
        double complex alpha = *palpha;
        
        if (alpha == 0) {
                vzCopy(n, x, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha + x[i];
                }
        }
}

void
vzNeg (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = -x[i];
        }
}

void
vzAbs (int n, const double complex *x, double complex *z)
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
vzSgn (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                zsgn(x[i], &(z[i]));
        }
}

void
vzInv (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = 1.0 / x[i];
        }
}

static void
vzAxpy (int n, const double complex *palpha, const double complex *x, const double complex *y, double complex *z)
{
        double complex alpha = *palpha;
        la_int one = 1;
        
        if (alpha == 0) {
                vzCopy(n, y, z);
        } else if (y == z) {
                F77_FUNC(zaxpy) (&n, palpha, x, &one, z, &one);
        } else if (alpha == 1 && x == z) {
                F77_FUNC(zaxpy) (&n, palpha, y, &one, z, &one);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i] + y[i];
                }
        }
}

void
vzAxpby (int n, const double complex *palpha, const double complex *x, const double complex *pbeta, const double complex *y, double complex *z)
{
        double complex alpha = *palpha;
        double complex beta = *pbeta;        
                
        if (alpha == 0) {
                vzScale(n, pbeta, y, z);
        } else if (alpha == 1) {
                vzAxpy(n, pbeta, y, x, z);
        } else if (beta == 0) {
                vzScale(n, palpha, x, z);
        } else if (beta == 1) {
                vzAxpy(n, palpha, x, y, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i] + beta * y[i];
                }
        }
}

void
vzAdd (int n, const double complex *x, const double complex *y, double complex *z)
{
        vdAdd(2 * n, (const double *)x, (const double *)y, (double *)z);
}

void vzSub (int n, const double complex *x, const double complex *y, double complex *z)
{
        vdSub(2 * n, (const double *)x, (const double *)y, (double *)z);
}

void
vzMul (int n, const double complex *x, const double complex *y, double complex *z)
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
vzDiv (int n, const double complex *x, const double complex *y, double complex *z)
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
vzExp (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cexp(x[i]);
        }
}

void
vzSqrt (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = csqrt(x[i]);
        }
}

void
vzLog (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = clog(x[i]);
        }        
}

void
vzPow (int n, const double complex *x, const double complex *y, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cpow(x[i], y[i]);
        }        
}

void
vzSin (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = csin(x[i]);
        }                
}

void
vzCos (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = ccos(x[i]);
        }                        
}

void
vzTan (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = ctan(x[i]);
        }                                
}

void
vzASin (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = casin(x[i]);
        }                                        
}

void
vzACos (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cacos(x[i]);
        }                                                
}

void
vzATan (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = catan(x[i]);
        }
}

void
vzSinh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = csinh(x[i]);
        }                
}

void
vzCosh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = ccosh(x[i]);
        }                        
}

void
vzTanh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = ctanh(x[i]);
        }                                
}

void
vzASinh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = casinh(x[i]);
        }                                        
}

void
vzACosh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cacosh(x[i]);
        }                                                
}

void
vzATanh (int n, const double complex *x, double complex *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = catanh(x[i]);
        }
}

