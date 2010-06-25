
#include <math.h>
#include <string.h>
#include "BLAS.h"
#include "vectorOps.h"

static void
dVectorClear (int n, double *z)
{
        memset(z, 0, n * sizeof(double));
}

static void
dVectorCopy (int n, const double *x, double *z)
{
        if (x != z) {
                memcpy(z, x, n * sizeof(double));
        }
}

void
dVectorScale (int n, double alpha, const double *x, double *z)
{
        if (alpha == 1) {
                dVectorCopy(n, x, z);
        } else if (alpha == 0) {
                dVectorClear(n, z);
        } else if (x == z) {
                blas_dscal(n, alpha, z, 1);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i];
                }
        }
}

void
dVectorShift (int n, double alpha, const double *x, double *z)
{
        if (alpha == 0) {
                dVectorCopy(n, x, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha + x[i];
                }
        }
}

void
dVectorNeg (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = -x[i];
        }
}

void
dVectorAbs (int n, const double *x, double *z)
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
dVectorSgn (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = dsgn(x[i]);
        }
}

void
dVectorInv (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = 1.0 / x[i];
        }
}

static void
dVectorAxpy (int n, double alpha, const double *x, const double *y, double *z)
{
        if (alpha == 0) {
                dVectorCopy(n, y, z);
        } else if (y == z) {
                blas_daxpy(n, alpha, x, 1, z, 1);
        } else if (alpha == 1 && x == z) {
                blas_daxpy(n, 1, y, 1, z, 1);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i] + y[i];
                }
        }
}

void
dVectorAxpby (int n, double alpha, const double *x, double beta, const double *y, double *z)
{
        if (alpha == 0) {
                dVectorScale(n, beta, y, z);
        } else if (alpha == 1) {
                dVectorAxpy(n, beta, y, x, z);
        } else if (beta == 0) {
                dVectorScale(n, alpha, x, z);
        } else if (beta == 1) {
                dVectorAxpy(n, alpha, x, y, z);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = alpha * x[i] + beta * y[i];
                }
        }
}

void
dVectorAdd (int n, const double *x, const double *y, double *z)
{
        dVectorAxpy(n, 1, x, y, z);
}

void dVectorSub (int n, const double *x, const double *y, double *z)
{
        dVectorAxpy(n, -1, y, x, z);
}

void
dVectorMul (int n, const double *x, const double *y, double *z)
{
        if (y == z) {
                blas_dtbmv(BlasUpper, BlasNoTrans, BlasNonUnit, n, 0, x, 1,
                           z, 1);
        } else if (x == z) {
                blas_dtbmv(BlasUpper, BlasNoTrans, BlasNonUnit, n, 0, y, 1,
                           z, 1);
        } else {
                int i;
                for (i = 0; i < n; i++) {
                        z[i] = x[i] * y[i];
                }
        }
}

void
dVectorDiv (int n, const double *x, const double *y, double *z)
{
        if (y == z) {
                blas_dtbsv(BlasUpper, BlasNoTrans, BlasNonUnit, n, 0, x, 1,
                           z, 1);
        } else if (x == z) {
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
dVectorExp (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = exp(x[i]);
        }
}

void
dVectorSqrt (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = sqrt(x[i]);
        }
}

void
dVectorLog (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = log(x[i]);
        }        
}

void
dVectorPow (int n, const double *x, const double *y, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = pow(x[i], y[i]);
        }        
}

void
dVectorSin (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = sin(x[i]);
        }                
}

void
dVectorCos (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cos(x[i]);
        }                        
}

void
dVectorTan (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = tan(x[i]);
        }                                
}

void
dVectorASin (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = asin(x[i]);
        }                                        
}

void
dVectorACos (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = acos(x[i]);
        }                                                
}

void
dVectorATan (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = atan(x[i]);
        }
}

void
dVectorSinh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = sinh(x[i]);
        }                
}

void
dVectorCosh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = cosh(x[i]);
        }                        
}

void
dVectorTanh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = tanh(x[i]);
        }                                
}

void
dVectorASinh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = asinh(x[i]);
        }                                        
}

void
dVectorACosh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = acosh(x[i]);
        }                                                
}

void
dVectorATanh (int n, const double *x, double *z)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = atanh(x[i]);
        }
}

