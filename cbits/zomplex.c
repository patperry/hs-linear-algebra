
#include "BLAS.h"
#include "config.h"
#include <stdio.h>

static char *BLAS_TRANS_CODES[] = { "N", "T", "C" };
#define TRANS(x) BLAS_TRANS_CODES[(int) (x) - (int) BlasNoTrans]

static char *BLAS_UPLO_CODES[] = { "U", "L" };
#define UPLO(x) BLAS_UPLO_CODES[(int) (x) - (int) BlasUpper]

static char *BLAS_DIAG_CODES[] = { "N", "U" };
#define DIAG(x) BLAS_DIAG_CODES[(int) (x) - (int) BlasNonUnit]

static char *BLAS_SIDE_CODES[] = { "L", "R" };
#define SIDE(x) BLAS_SIDE_CODES[(int) (x) - (int) BlasLeft]


extern double
F77_FUNC(dzasum)(const int *n, const void *zx, const int *incx);
extern double
F77_FUNC(dznrm2)(const int *n, const void *x, const int *incx);
extern int
F77_FUNC(izamax)(const int *n, const void *zx, const int *incx);
extern void
F77_FUNC(zaxpy)(const int *n, const void *za, const void *zx,
	    const int *incx, void *zy, const int *incy);
extern void 
F77_FUNC(zcopy)(const int *n, const void *zx, const int *incx,
	    void *zy, const int *incy);
extern const void 
F77_FUNC(zdotc)(const void * ret_val, const int *n,
	    const void *zx, const int *incx, const void *zy, const int *incy);
extern const void 
F77_FUNC(zdotu)(const void * ret_val, const int *n,
	    const void *zx, const int *incx, const void *zy, const int *incy);
extern void 
F77_FUNC(zdrot)(const int *n, void *zx, const int *incx, void *zy,
	const int *incy, const double *c, const double *s);
extern void 
F77_FUNC(zdscal)(const int *n, const double *da, void *zx, const int *incx);
extern void 
F77_FUNC(zgbmv)(const char *trans, const int *m, const int *n, const int *kl,
	    const int *ku, const void *alpha, const void *a, const int *lda,
	    const void *x, const int *incx, const void *beta, void *y,
	    const int *incy);
extern void
F77_FUNC(zgemm)(const char *transa, const char *transb, const int *m,
	    const int *n, const int *k, const const void *alpha,
	    const const void *a, const int *lda,
	    const const void *b, const int *ldb,
	    const const void *beta, const void *c, const int *ldc);
extern void 
F77_FUNC(zgemv)(const char *trans, const int *m, const int *n, const void *alpha,
	    const void *a, const int *lda, const void *x, const int *incx,
	    const void *beta, void *y, const int *incy);
extern void 
F77_FUNC(zgerc)(const int *m, const int *n, const void *alpha, const void *x,
	    const int *incx, const void *y, const int *incy, void *a, const int *lda);
extern void 
F77_FUNC(zgeru)(const int *m, const int *n, const void *alpha, const void *x,
	    const int *incx, const void *y, const int *incy, void *a, const int *lda);
extern void 
F77_FUNC(zhbmv)(const char *uplo, const int *n, const int *k, const void *alpha,
	    const void *a, const int *lda, const void *x, const int *incx,
	    const void *beta, const void *y, const int *incy);
extern void 
F77_FUNC(zhemm)(const char *side, const char *uplo, const int *m, const int *n,
	    const void *alpha, const void *a, const int *lda, const void *b,
	    const int *ldb, const void *beta, const void *c, const int *ldc);
extern void 
F77_FUNC(zhemv)(const char *uplo, const int *n, const void *alpha, const void *a,
	    const int *lda, const void *x, const int *incx, const void *beta,
	    const void *y, const int *incy);
extern void 
F77_FUNC(zher)(const char *uplo, const int *n, const double *alpha, const void *x,
	   const int *incx, void *a, const int *lda);
extern void 
F77_FUNC(zher2)(const char *uplo, const int *n, const void *alpha, const void *x,
	    const int *incx, const void *y, const int *incy, void *a, const int *lda);
extern void 
F77_FUNC(zher2k)(const char *uplo, const char *trans, const int *n, const int *k,
	     const void *alpha, const void *a, const int *lda, const void *b,
	     const int *ldb, const double *beta, void *c, const int *ldc);
extern void 
F77_FUNC(zherk)(const char *uplo, const char *trans, const int *n, const int *k,
	    const double *alpha, const void *a, const int *lda, const double *beta,
	    void *c, const int *ldc);
extern void 
F77_FUNC(zhpmv)(const char *uplo, const int *n, const void *alpha, const void *ap,
	    const void *x, const int *incx, const void * beta, void *y,
	    const int *incy);
extern void 
F77_FUNC(zhpr)(const char *uplo, const int *n, const double *alpha,
	   const void *x, const int *incx, void *ap);
extern void 
F77_FUNC(zhpr2)(const char *uplo, const int *n, const void *alpha, const void *x,
	    const int *incx, const void *y, const int *incy, void *ap);
extern void 
F77_FUNC(zrotg)(const void *ca, const void *cb, void *c, const void *s);
extern void 
F77_FUNC(zscal)(const int *n, const void *za, const void *zx, const int *incx);
extern void 
F77_FUNC(zswap)(const int *n, const void *zx, const int *incx, const void *zy, const int *incy);
extern void 
F77_FUNC(zsymm)(const char *side, const char *uplo, const int *m, const int *n,
	    const void *alpha, const void *a, const int *lda, const void *b,
	    const int *ldb, const void *beta, const void *c, const int *ldc);
extern void 
F77_FUNC(zsyr2k)(const char *uplo, const char *trans, const int *n, const int *k,
	     const void *alpha, const void *a, const int *lda, const void *b,
	     const int *ldb, const void *beta, void *c, const int *ldc);
extern void 
F77_FUNC(zsyrk)(const char *uplo, const char *trans, const int *n, const int *k,
	    const void *alpha, const void *a, const int *lda,
	    const void *beta, void *c, const int *ldc);
extern void 
F77_FUNC(ztbmv)(const char *uplo, const char *trans, const char *diag, const int *n, const int *k,
	    const void *a, const int *lda, void *x, const int *incx);
extern void 
F77_FUNC(ztbsv)(const char *uplo, const char *trans, const char *diag, const int *n, const int *k,
	    const void *a, const int *lda, void *x, const int *incx);
extern void 
F77_FUNC(ztpmv)(const char *uplo, const char *trans, const char *diag, const int *n,
	    const void *ap, void *x, const int *incx);
extern void 
F77_FUNC(ztpsv)(const char *uplo, const char *trans, const char *diag, const int *n,
	    const void *ap, void *x, const int *incx);
extern void 
F77_FUNC(ztrmm)(const char *side, const char *uplo, const char *transa, const char *diag,
	    const int *m, const int *n, const void *alpha, const void *a,
	    const int *lda, void *b, const int *ldb);
extern void 
F77_FUNC(ztrmv)(const char *uplo, const char *trans, const char *diag, const int *n,
	    const void *a, const int *lda, void *x, const int *incx);
extern void 
F77_FUNC(ztrsm)(const char *side, const char *uplo, const char *transa, const char *diag,
	    const int *m, const int *n, const void *alpha, const void *a,
	    const int *lda, void *b, const int *ldb);
extern void 
F77_FUNC(ztrsv)(const char *uplo, const char *trans, const char *diag, const int *n,
	    const void *a, const int *lda, void *x, const int *incx);


void 
blas_zdotu_sub (const int N, const void *X, const int incX,
                const void *Y, const int incY, void *dotu)
{
    return F77_FUNC(zdotu) (dotu, &N, X, &incX, Y, &incY);
}

void 
blas_zdotc_sub (const int N, const void *X, const int incX,
                const void *Y, const int incY, void *dotc)
{
    return F77_FUNC(zdotc) (dotc, &N, X, &incX, Y, &incY);
}


double 
blas_dznrm2 (const int N, const void *X, const int incX)
{
    return F77_FUNC(dznrm2) (&N, X, &incX);
}

double 
blas_dzasum (const int N, const void *X, const int incX)
{
    return F77_FUNC(dzasum) (&N, X, &incX);    
}

BLAS_INDEX 
blas_izamax (const int N, const void *X, const int incX)
{
    return (F77_FUNC(izamax) (&N, X, &incX) - 1);
}

void 
blas_zswap (const int N, void *X, const int incX,
            void *Y, const int incY)
{
    F77_FUNC(zswap) (&N, X, &incX, Y, &incY);
}

void 
blas_zcopy (const int N, const void *X, const int incX,
            void *Y, const int incY)
{
    F77_FUNC(zcopy) (&N, X, &incX, Y, &incY);
}

void 
blas_zaxpy (const int N, const void *alpha, const void *X,
            const int incX, void *Y, const int incY)
{
    F77_FUNC(zaxpy) (&N, alpha, X, &incX, Y, &incY);
}

void 
blas_zdscal (const int N, const double alpha, void *X, const int incX)
{
    F77_FUNC(zdscal) (&N, &alpha, X, &incX);
}

void 
blas_zscal (const int N, const void *alpha, void *X, const int incX)
{
    F77_FUNC(zscal) (&N, alpha, X, &incX);
}

void 
blas_zrotg (void *a, void *b, void *c, void *s)
{
    F77_FUNC(zrotg) (a, b, c, s);
}

void 
blas_zdrot (const int N, void *X, const int incX,
           void *Y, const int incY, const double c, const double s)
{
    F77_FUNC(zdrot) (&N, X, &incX, Y, &incY, &c, &s);
}


void 
blas_zgemv (const enum BLAS_TRANSPOSE TransA, const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX, const void *beta,
            void *Y, const int incY)
{
    F77_FUNC(zgemv) (TRANS(TransA), &M, &N, alpha, A, &lda, X, &incX,
                     beta, Y, &incY);
}

void 
blas_zgbmv (const enum BLAS_TRANSPOSE TransA, const int M, const int N,
            const int KL, const int KU, const void *alpha,
            const void *A, const int lda, const void *X,
            const int incX, const void *beta, void *Y, const int incY)
{
    F77_FUNC(zgbmv) (TRANS(TransA), &M, &N, &KL, &KU, alpha, A, &lda,
                     X, &incX, beta, Y, &incY);
}
 
void 
blas_ztrmv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const void *A, const int lda,
            void *X, const int incX)
{
    F77_FUNC(ztrmv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, A, &lda,
                     X, &incX);
}

void 
blas_ztbmv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const int K, const void *A, const int lda,
            void *X, const int incX)
{
    F77_FUNC(ztbmv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, &K, A, &lda,
                     X, &incX);
}

void 
blas_ztpmv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const void *Ap, void *X, const int incX)
{
    F77_FUNC(ztpmv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, Ap,
                     X, &incX);
}

void 
blas_ztrsv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const void *A, const int lda, void *X,
            const int incX)
{
    F77_FUNC(ztrsv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, A, &lda,
                     X, &incX);
}

void 
blas_ztbsv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const int K, const void *A, const int lda,
            void *X, const int incX)

{
    F77_FUNC(ztbsv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, &K, A, &lda,
                     X, &incX);
}

void 
blas_ztpsv (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
            const int N, const void *Ap, void *X, const int incX)
{
    F77_FUNC(ztpsv) (UPLO(Uplo), TRANS(TransA), DIAG(Diag), &N, Ap,
                     X, &incX);
}


void 
blas_zhemv (const enum BLAS_UPLO Uplo,
            const int N, const void *alpha, const void *A,
            const int lda, const void *X, const int incX,
            const void *beta, void *Y, const int incY)
{
    F77_FUNC(zhemv) (UPLO(Uplo), &N, alpha, A, &lda, X, &incX,
                     beta, Y, &incY);    
}

void 
blas_zhbmv (const enum BLAS_UPLO Uplo,
            const int N, const int K, const void *alpha, const void *A,
            const int lda, const void *X, const int incX,
            const void *beta, void *Y, const int incY)
{
    F77_FUNC(zhbmv) (UPLO(Uplo), &N, &K, alpha, A, &lda, X, &incX,
                     beta, Y, &incY);    
}

void
blas_zhpmv (const enum BLAS_UPLO Uplo,
            const int N, const void *alpha, const void *Ap,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY)
{
    F77_FUNC(zhpmv) (UPLO(Uplo), &N, alpha, Ap, X, &incX, beta, Y, &incY);
}

void
blas_zgeru (const int M, const int N,
           const void *alpha, const void *X, const int incX,
           const void *Y, const int incY, void *A, const int lda)
{
    F77_FUNC(zgeru) (&M, &N, alpha, X, &incX, Y, &incY, A, &lda);        
}

void
blas_zgerc (const int M, const int N,
           const void *alpha, const void *X, const int incX,
           const void *Y, const int incY, void *A, const int lda)
{
    F77_FUNC(zgerc) (&M, &N, alpha, X, &incX, Y, &incY, A, &lda);        
}

void
blas_zher (const enum BLAS_UPLO Uplo,
           const int N, const double alpha, const void *X,
           const int incX, void *A, const int lda)
{
    F77_FUNC(zher) (UPLO(Uplo), &N, &alpha, X, &incX, A, &lda);    
}

void
blas_zhpr (const enum BLAS_UPLO Uplo,
           const int N, const double alpha, const void *X,
           const int incX, void *Ap)
{
    F77_FUNC(zhpr) (UPLO(Uplo), &N, &alpha, X, &incX, Ap);    
}

void 
blas_zher2 (const enum BLAS_UPLO Uplo,
            const int N, const void * alpha, const void *X,
            const int incX, const void *Y, const int incY, void *A,
            const int lda)
{
    F77_FUNC(zher2) (UPLO(Uplo), &N, alpha, X, &incX, Y, &incY, A, &lda);    
}

void
blas_zhpr2 (const enum BLAS_UPLO Uplo,
            const int N, const void * alpha, const void *X,
            const int incX, const void *Y, const int incY, void *A)
{
    F77_FUNC(zhpr2) (UPLO(Uplo), &N, alpha, X, &incX, Y, &incY, A);
}

void 
blas_zgemm (const enum BLAS_TRANSPOSE TransA,
            const enum BLAS_TRANSPOSE TransB, const int M, const int N,
            const int K, const void * alpha, const void *A,
            const int lda, const void *B, const int ldb,
            const void * beta, void *C, const int ldc)
{
    F77_FUNC(zgemm) (TRANS(TransA), TRANS(TransB), &M, &N, &K,
                     alpha, A, &lda, B, &ldb, beta, C, &ldc);
}

void 
blas_zsymm (const enum BLAS_SIDE Side,
            const enum BLAS_UPLO Uplo, const int M, const int N,
            const void * alpha, const void *A, const int lda,
            const void *B, const int ldb, const void * beta,
            void *C, const int ldc)
{
    F77_FUNC(zsymm) (SIDE(Side), UPLO(Uplo), &M, &N,
                     alpha, A, &lda, B, &ldb, beta, C, &ldc);
}

void 
blas_zsyrk (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE Trans, const int N, const int K,
            const void * alpha, const void *A, const int lda,
            const void * beta, void *C, const int ldc)
{
    F77_FUNC(zsyrk) (UPLO(Uplo), TRANS(Trans), &N, &K,
                     alpha, A, &lda, beta, C, &ldc);
    
}

void
blas_zsyr2k (const enum BLAS_UPLO Uplo,
             const enum BLAS_TRANSPOSE Trans, const int N, const int K,
             const void * alpha, const void *A, const int lda,
             const void *B, const int ldb, const void * beta,
             void *C, const int ldc)
{
    F77_FUNC(zsyr2k) (UPLO(Uplo), TRANS(Trans), &N, &K,
                      alpha, A, &lda, B, &ldb, beta, C, &ldc);
}

void 
blas_zhemm (const enum BLAS_SIDE Side,
            const enum BLAS_UPLO Uplo, const int M, const int N,
            const void * alpha, const void *A, const int lda,
            const void *B, const int ldb, const void * beta,
            void *C, const int ldc)
{
    F77_FUNC(zhemm) (SIDE(Side), UPLO(Uplo), &M, &N,
                     alpha, A, &lda, B, &ldb, beta, C, &ldc);
}

void 
blas_zherk (const enum BLAS_UPLO Uplo,
            const enum BLAS_TRANSPOSE Trans, const int N, const int K,
            const double alpha, const void *A, const int lda,
            const double beta, void *C, const int ldc)
{
    F77_FUNC(zherk) (UPLO(Uplo), TRANS(Trans), &N, &K,
                     &alpha, A, &lda, &beta, C, &ldc);
    
}

void
blas_zher2k (const enum BLAS_UPLO Uplo,
             const enum BLAS_TRANSPOSE Trans, const int N, const int K,
             const void * alpha, const void *A, const int lda,
             const void *B, const int ldb, const double beta,
             void *C, const int ldc)
{
    F77_FUNC(zher2k) (UPLO(Uplo), TRANS(Trans), &N, &K,
                      alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

void
blas_ztrmm (const enum BLAS_SIDE Side,
            const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
            const enum BLAS_DIAG Diag, const int M, const int N,
            const void * alpha, const void *A, const int lda,
            void *B, const int ldb)
{
    F77_FUNC(ztrmm) (SIDE(Side), UPLO(Uplo), TRANS(TransA), DIAG(Diag), &M, &N,
                     alpha, A, &lda, B, &ldb);
}

void
blas_ztrsm (const enum BLAS_SIDE Side,
            const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
            const enum BLAS_DIAG Diag, const int M, const int N,
            const void * alpha, const void *A, const int lda,
            void *B, const int ldb)
{
    F77_FUNC(ztrsm) (SIDE(Side), UPLO(Uplo), TRANS(TransA), DIAG(Diag), &M, &N,
                     alpha, A, &lda, B, &ldb);
}
