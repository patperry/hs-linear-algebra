
#include "BLAS.h"
#include "config.h"

static char *BLAS_TRANS_CODES[] = { "N", "T", "C" };
#define TRANS(x) BLAS_TRANS_CODES[(int) (x) - (int) BlasNoTrans]

static char *BLAS_UPLO_CODES[] = { "U", "L" };
#define UPLO(x) BLAS_UPLO_CODES[(int) (x) - (int) BlasUpper]

static char *BLAS_DIAG_CODES[] = { "N", "U" };
#define DIAG(x) BLAS_DIAG_CODES[(int) (x) - (int) BlasNonUnit]

static char *BLAS_SIDE_CODES[] = { "L", "R" };
#define SIDE(x) BLAS_SIDE_CODES[(int) (x) - (int) BlasLeft]


extern double
F77_FUNC(ddot) (const int *N, const double *X, const int *incX,
	            const double *Y, const int *incY);
extern double
F77_FUNC(dnrm2) (const int *N, const double *X, const int *incX);
extern double
F77_FUNC(dasum) (const int *N, const double *X, const int *incX);

extern BLAS_INDEX 
F77_FUNC(idamax) (const int *N, const double *X, const int *incX);
extern void 
F77_FUNC(dswap) (const int *N, double *X, const int *incX,
                 double *Y, const int *incY);
extern void 
F77_FUNC(dcopy) (const int *N, const double *X, const int *incX,
                 double *Y, const int *incY);
extern void 
F77_FUNC(daxpy) (const int *N, const double *alpha, const double *X,
                 const int *incX, double *Y, const int *incY);

extern void
F77_FUNC(drotg) (double *a, double *b, double *c, double *s);
extern void
F77_FUNC(drotmg) (double *d1, double *d2, double *b1, const double *b2, double *P);
extern void
F77_FUNC(drot) (const int *N, double *X, const int *incX,
                double *Y, const int *incY, const double *c, const double *s);
extern void 
F77_FUNC(drotm) (const int *N, double *X, const int *incX,
                 double *Y, const int *incY, const double *P);
extern void 
F77_FUNC(dscal) (const int *N, const double *alpha, double *X, const int *incX);


extern void 
F77_FUNC(dgemv) (const char *TransA, const int *M, const int *N,
                 const double *alpha, const double *A, const int *lda,
                 const double *X, const int *incX, const double *beta,
                 double *Y, const int *incY);

extern void
F77_FUNC(dgbmv)(const char *trans, const int *m, const int *n,
		const int *kl,const int *ku,
		const double *alpha, const double *a, const int *lda,
		const double *x, const int *incx,
		const double *beta, double *y, const int *incy);

extern void
F77_FUNC(dsbmv)(const char *uplo, const int *n, const int *k,
		const double *alpha, const double *a, const int *lda,
		const double *x, const int *incx,
		const double *beta, double *y, const int *incy);
extern void
F77_FUNC(dspmv)(const char *uplo, const int *n,
		const double *alpha, const double *ap,
		const double *x, const int *incx,
		const double *beta, double *y, const int *incy);
extern void
F77_FUNC(dsymv)(const char *uplo, const int *n, const double *alpha,
		const double *a, const int *lda,
		const double *x, const int *incx,
		const double *beta, double *y, const int *incy);
extern void
F77_FUNC(dtbmv)(const char *uplo, const char *trans,
		const char *diag, const int *n, const int *k,
		const double *a, const int *lda,
		double *x, const int *incx);
extern void
F77_FUNC(dtpmv)(const char *uplo, const char *trans, const char *diag,
		const int *n, const double *ap,
		double *x, const int *incx);
extern void
F77_FUNC(dtrmv)(const char *uplo, const char *trans, const char *diag,
		const int *n, const double *a, const int *lda,
		double *x, const int *incx);
extern void
F77_FUNC(dtbsv)(const char *uplo, const char *trans,
		const char *diag, const int *n, const int *k,
		const double *a, const int *lda,
		double *x, const int *incx);
extern void
F77_FUNC(dtpsv)(const char *uplo, const char *trans,
		const char *diag, const int *n,
		const double *ap, double *x, const int *incx);
extern void
F77_FUNC(dtrsv)(const char *uplo, const char *trans,
		const char *diag, const int *n,
		const double *a, const int *lda,
		double *x, const int *incx);
extern void
F77_FUNC(dger)(const int *m, const int *n, const double *alpha,
	       const double *x, const int *incx,
	       const double *y, const int *incy,
	       double *a, const int *lda);
extern void
F77_FUNC(dsyr)(const char *uplo, const int *n, const double *alpha,
	       const double *x, const int *incx,
	       double *a, const int *lda);
extern void
F77_FUNC(dspr)(const char *uplo, const int *n, const double *alpha,
	       const double *x, const int *incx, double *ap);
extern void
F77_FUNC(dsyr2)(const char *uplo, const int *n, const double *alpha,
		const double *x, const int *incx,
		const double *y, const int *incy,
		double *a, const int *lda);
extern void
F77_FUNC(dspr2)(const char *uplo, const int *n, const double *alpha,
		const double *x, const int *incx,
		const double *y, const int *incy, double *ap);
extern void
F77_FUNC(dgemm)(const char *transa, const char *transb, const int *m,
		const int *n, const int *k, const double *alpha,
		const double *a, const int *lda,
		const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);
extern void
F77_FUNC(dtrsm)(const char *side, const char *uplo,
		const char *transa, const char *diag,
		const int *m, const int *n, const double *alpha,
		const double *a, const int *lda,
		double *b, const int *ldb);
extern void
F77_FUNC(dtrmm)(const char *side, const char *uplo, const char *transa,
		const char *diag, const int *m, const int *n,
		const double *alpha, const double *a, const int *lda,
		double *b, const int *ldb);
extern void
F77_FUNC(dsymm)(const char *side, const char *uplo, const int *m,
		const int *n, const double *alpha,
		const double *a, const int *lda,
		const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);
extern void
F77_FUNC(dsyrk)(const char *uplo, const char *trans,
		const int *n, const int *k,
		const double *alpha, const double *a, const int *lda,
		const double *beta, double *c, const int *ldc);
extern void
F77_FUNC(dsyr2k)(const char *uplo, const char *trans,
		 const int *n, const int *k,
		 const double *alpha, const double *a, const int *lda,
		 const double *b, const int *ldb,
		 const double *beta, double *c, const int *ldc);


double 
blas_ddot (const int N, const double *X, const int incX,
           const double *Y, const int incY)
{
    return F77_FUNC(ddot) (&N, X, &incX, Y, &incY);
}


double 
blas_dnrm2 (const int N, const double *X, const int incX)
{
    return F77_FUNC(dnrm2) (&N, X, &incX);
}

double 
blas_dasum (const int N, const double *X, const int incX)
{
    return F77_FUNC(dasum) (&N, X, &incX);    
}

BLAS_INDEX 
blas_idamax (const int N, const double *X, const int incX)
{
    return (F77_FUNC(idamax) (&N, X, &incX) - 1);
}

void 
blas_dswap (const int N, double *X, const int incX,
            double *Y, const int incY)
{
    F77_FUNC(dswap) (&N, X, &incX, Y, &incY);
}

void 
blas_dcopy (const int N, const double *X, const int incX,
            double *Y, const int incY)
{
    F77_FUNC(dcopy) (&N, X, &incX, Y, &incY);
}

void 
blas_daxpy (const int N, const double alpha, const double *X,
            const int incX, double *Y, const int incY)
{
    F77_FUNC(daxpy) (&N, &alpha, X, &incX, Y, &incY);
}

void 
blas_drotg (double *a, double *b, double *c, double *s)
{
    F77_FUNC(drotg) (a, b, c, s);
}

void 
blas_drotmg (double *d1, double *d2, double *b1, const double b2, double *P)
{
    F77_FUNC(drotmg) (d1, d2, b1, &b2, P);
}

void 
blas_drot (const int N, double *X, const int incX,
           double *Y, const int incY, const double c, const double s)
{
    F77_FUNC(drot) (&N, X, &incX, Y, &incY, &c, &s);
}

void 
blas_drotm (const int N, double *X, const int incX,
            double *Y, const int incY, const double *P)
{
    F77_FUNC(drotm) (&N, X, &incX, Y, &incY, P);
}

void 
blas_dscal (const int N, const double alpha, double *X, const int incX)
{
    F77_FUNC(dscal) (&N, &alpha, X, &incX);
}

void 
blas_dgemv (const enum BLAS_TRANSPOSE TransA, const int M, const int N,
            const double alpha, const double *A, const int lda,
            const double *X, const int incX, const double beta,
            double *Y, const int incY)
{
    F77_FUNC(dgemv) (TRANS(TransA), &M, &N, &alpha, A, &lda, X, &incX,
                     &beta, Y, &incY);
}

void 
blas_dgbmv (const enum BLAS_TRANSPOSE TransA, const int M, const int N,
            const int KL, const int KU, const double alpha,
            const double *A, const int lda, const double *X,
            const int incX, const double beta, double *Y, const int incY)
{
    F77_FUNC(dgbmv) (TRANS(TransA), &M, &N, &KL, &KU, &alpha, A, &lda,
                     X, &incX, &beta, Y, &incY);
}
 
void 
blas_dtrmv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const double *A, const int lda,
            double *X, const int incX)
{
    F77_FUNC(dtrmv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, A, &lda,
                     X, &incX);
}

void 
blas_dtbmv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const int K, const double *A, const int lda,
            double *X, const int incX)
{
    F77_FUNC(dtbmv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, &K, A, &lda,
                     X, &incX);
}

void 
blas_dtpmv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const double *Ap, double *X, const int incX)
{
    F77_FUNC(dtpmv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, Ap,
                     X, &incX);
}

void 
blas_dtrsv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const double *A, const int lda, double *X,
            const int incX)
{
    F77_FUNC(dtrsv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, A, &lda,
                     X, &incX);
}

void 
blas_dtbsv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const int K, const double *A, const int lda,
            double *X, const int incX)

{
    F77_FUNC(dtbsv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, &K, A, &lda,
                     X, &incX);
}

void 
blas_dtpsv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const double *Ap, double *X, const int incX)
{
    F77_FUNC(dtpsv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, Ap,
                     X, &incX);
}


void 
blas_dsymv (const enum BLAS_UPLO Uplo,
            const int N, const double alpha, const double *A,
            const int lda, const double *X, const int incX,
            const double beta, double *Y, const int incY)
{
    F77_FUNC(dsymv) (UPLO(Uplo), &N, &alpha, A, &lda, X, &incX,
                     &beta, Y, &incY);    
}

void 
blas_dsbmv (const enum BLAS_UPLO Uplo,
            const int N, const int K, const double alpha, const double *A,
            const int lda, const double *X, const int incX,
            const double beta, double *Y, const int incY)
{
    F77_FUNC(dsbmv) (UPLO(Uplo), &N, &K, &alpha, A, &lda, X, &incX,
                     &beta, Y, &incY);    
}

void
blas_dspmv (const enum BLAS_UPLO Uplo,
            const int N, const double alpha, const double *Ap,
            const double *X, const int incX,
            const double beta, double *Y, const int incY)
{
    F77_FUNC(dspmv) (UPLO(Uplo), &N, &alpha, Ap, X, &incX, &beta, Y, &incY);
}

void
blas_dger (const int M, const int N,
           const double alpha, const double *X, const int incX,
           const double *Y, const int incY, double *A, const int lda)
{
    F77_FUNC(dger) (&M, &N, &alpha, X, &incX, Y, &incY, A, &lda);        
}

void
blas_dsyr (const enum BLAS_UPLO Uplo,
           const int N, const double alpha, const double *X,
           const int incX, double *A, const int lda)
{
    F77_FUNC(dsyr) (UPLO(Uplo), &N, &alpha, X, &incX, A, &lda);    
}

void
blas_dspr (const enum BLAS_UPLO Uplo,
           const int N, const double alpha, const double *X,
           const int incX, double *Ap)
{
    F77_FUNC(dspr) (UPLO(Uplo), &N, &alpha, X, &incX, Ap);    
}

void 
blas_dsyr2 (const enum BLAS_UPLO Uplo,
            const int N, const double alpha, const double *X,
            const int incX, const double *Y, const int incY, double *A,
            const int lda)
{
    F77_FUNC(dsyr2) (UPLO(Uplo), &N, &alpha, X, &incX, Y, &incY, A, &lda);    
}

void
blas_dspr2 (const enum BLAS_UPLO Uplo,
            const int N, const double alpha, const double *X,
            const int incX, const double *Y, const int incY, double *A)
{
    F77_FUNC(dspr2) (UPLO(Uplo), &N, &alpha, X, &incX, Y, &incY, A);
}

void 
blas_dgemm (const enum BLAS_TRANSPOSE TransA,
            const enum BLAS_TRANSPOSE TransB, const int M, const int N,
            const int K, const double alpha, const double *A,
            const int lda, const double *B, const int ldb,
            const double beta, double *C, const int ldc)
{
    F77_FUNC(dgemm) (TRANS(TransA), TRANS(TransB), &M, &N, &K,
                     &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

void 
blas_dsymm (const enum BLAS_SIDE Side,
            const enum BLAS_UPLO Uplo, const int M, const int N,
            const double alpha, const double *A, const int lda,
            const double *B, const int ldb, const double beta,
            double *C, const int ldc)
{
    F77_FUNC(dsymm) (SIDE(Side), UPLO(Uplo), &M, &N,
                     &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

void 
blas_dsyrk (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE Trans, const int N, const int K,
            const double alpha, const double *A, const int lda,
            const double beta, double *C, const int ldc)
{
    F77_FUNC(dsyrk) (UPLO(Uplo), TRANS(Trans), &N, &K,
                     &alpha, A, &lda, &beta, C, &ldc);
    
}

void
blas_dsyr2k (const enum BLAS_UPLO Uplo,
             const enum BLAS_TRANSPOSE Trans, const int N, const int K,
             const double alpha, const double *A, const int lda,
             const double *B, const int ldb, const double beta,
             double *C, const int ldc)
{
    F77_FUNC(dsyr2k) (UPLO(Uplo), TRANS(Trans), &N, &K,
                      &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

void
blas_dtrmm (const enum BLAS_SIDE Side,
            const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
            const enum BLAS_DIAG Diag, const int M, const int N,
            const double alpha, const double *A, const int lda,
            double *B, const int ldb)
{
    F77_FUNC(dtrmm) (SIDE(Side), UPLO(Uplo), TRANS(TransA), DIAG(Diag), &M, &N,
                     &alpha, A, &lda, B, &ldb);
}

void
blas_dtrsm (const enum BLAS_SIDE Side,
            const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
            const enum BLAS_DIAG Diag, const int M, const int N,
            const double alpha, const double *A, const int lda,
            double *B, const int ldb)
{
    F77_FUNC(dtrsm) (SIDE(Side), UPLO(Uplo), TRANS(TransA), DIAG(Diag), &M, &N,
                     &alpha, A, &lda, B, &ldb);
}
