#ifndef BLAS_H
#define BLAS_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef BLAS_ENUM_DEFINED_H
   #define BLAS_ENUM_DEFINED_H
   enum BLAS_TRANSPOSE {BlasNoTrans=111, BlasTrans=112, BlasConjTrans=113};
   enum BLAS_UPLO  {BlasUpper=121, BlasLower=122};
   enum BLAS_DIAG  {BlasNonUnit=131, BlasUnit=132};
   enum BLAS_SIDE  {BlasLeft=141, BlasRight=142};
#endif

#define BLAS_INDEX int

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */
float  blas_sdsdot(const int N, const float alpha, const float *X,
                    const int incX, const float *Y, const int incY);
double blas_dsdot(const int N, const float *X, const int incX, const float *Y,
                   const int incY);
float  blas_sdot(const int N, const float  *X, const int incX,
                  const float  *Y, const int incY);
double blas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);
/*
 * Functions having prefixes Z and C only
 */
void   blas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu);
void   blas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc);

void   blas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu);
void   blas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc);


/*
 * Functions having prefixes S D SC DZ
 */
float  blas_snrm2(const int N, const float *X, const int incX);
float  blas_sasum(const int N, const float *X, const int incX);

double blas_dnrm2(const int N, const double *X, const int incX);
double blas_dasum(const int N, const double *X, const int incX);

float  blas_scnrm2(const int N, const void *X, const int incX);
float  blas_scasum(const int N, const void *X, const int incX);

double blas_dznrm2(const int N, const void *X, const int incX);
double blas_dzasum(const int N, const void *X, const int incX);


/*
 * Functions having standard 4 prefixes (S D C Z)
 */
BLAS_INDEX blas_isamax(const int N, const float  *X, const int incX);
BLAS_INDEX blas_idamax(const int N, const double *X, const int incX);
BLAS_INDEX blas_icamax(const int N, const void   *X, const int incX);
BLAS_INDEX blas_izamax(const int N, const void   *X, const int incX);

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void blas_sswap(const int N, float *X, const int incX,
                 float *Y, const int incY);
void blas_scopy(const int N, const float *X, const int incX,
                 float *Y, const int incY);
void blas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY);

void blas_dswap(const int N, double *X, const int incX,
                 double *Y, const int incY);
void blas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY);
void blas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);
                       
void blas_cswap(const int N, void *X, const int incX,
                 void *Y, const int incY);
void blas_ccopy(const int N, const void *X, const int incX,
                 void *Y, const int incY);
void blas_caxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);

void blas_zswap(const int N, void *X, const int incX,
                 void *Y, const int incY);
void blas_zcopy(const int N, const void *X, const int incX,
                 void *Y, const int incY);
void blas_zaxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);

/*
 * Routines with S and D prefix only
 */
void blas_srotg(float *a, float *b, float *c, float *s);
void blas_srotmg(float *d1, float *d2, float *b1, const float b2, float *P);
void blas_srot(const int N, float *X, const int incX,
                float *Y, const int incY, const float c, const float s);
void blas_srotm(const int N, float *X, const int incX,
                float *Y, const int incY, const float *P);

void blas_drotg(double *a, double *b, double *c, double *s);
void blas_drotmg(double *d1, double *d2, double *b1, const double b2, double *P);
void blas_drot(const int N, double *X, const int incX,
                double *Y, const int incY, const double c, const double s);
void blas_drotm(const int N, double *X, const int incX,
                double *Y, const int incY, const double *P);


/*
 * Routines with S D C Z CS and ZD prefixes
 */
void blas_sscal(const int N, const float alpha, float *X, const int incX);
void blas_dscal(const int N, const double alpha, double *X, const int incX);
void blas_cscal(const int N, const void *alpha, void *X, const int incX);
void blas_zscal(const int N, const void *alpha, void *X, const int incX);
void blas_csscal(const int N, const float alpha, void *X, const int incX);
void blas_zdscal(const int N, const double alpha, void *X, const int incX);

/*
 * Extra reference routines provided by ATLAS, but not mandated by the standard
 */
void blas_crotg(void *a, void *b, void *c, void *s);
void blas_zrotg(void *a, void *b, void *c, void *s);
void blas_csrot(const int N, void *X, const int incX, void *Y, const int incY,
                 const float c, const float s);
void blas_zdrot(const int N, void *X, const int incX, void *Y, const int incY,
                 const double c, const double s);

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void blas_sgemv(
                 const enum BLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY);
void blas_sgbmv(
                 const enum BLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const float alpha,
                 const float *A, const int lda, const float *X,
                 const int incX, const float beta, float *Y, const int incY);
void blas_strmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const float *A, const int lda,
                 float *X, const int incX);
void blas_stbmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const int K, const float *A, const int lda,
                 float *X, const int incX);
void blas_stpmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const float *Ap, float *X, const int incX);
void blas_strsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const float *A, const int lda, float *X,
                 const int incX);
void blas_stbsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const int K, const float *A, const int lda,
                 float *X, const int incX);
void blas_stpsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const float *Ap, float *X, const int incX);

void blas_dgemv(
                 const enum BLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
void blas_dgbmv(
                 const enum BLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const double alpha,
                 const double *A, const int lda, const double *X,
                 const int incX, const double beta, double *Y, const int incY);
void blas_dtrmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const double *A, const int lda,
                 double *X, const int incX);
void blas_dtbmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const int K, const double *A, const int lda,
                 double *X, const int incX);
void blas_dtpmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX);
void blas_dtrsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const double *A, const int lda, double *X,
                 const int incX);
void blas_dtbsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const int K, const double *A, const int lda,
                 double *X, const int incX);
void blas_dtpsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX);

void blas_cgemv(
                 const enum BLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY);
void blas_cgbmv(
                 const enum BLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const void *alpha,
                 const void *A, const int lda, const void *X,
                 const int incX, const void *beta, void *Y, const int incY);
void blas_ctrmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const void *A, const int lda,
                 void *X, const int incX);
void blas_ctbmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void blas_ctpmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);
void blas_ctrsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const void *A, const int lda, void *X,
                 const int incX);
void blas_ctbsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void blas_ctpsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);

void blas_zgemv(
                 const enum BLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY);
void blas_zgbmv(
                 const enum BLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const void *alpha,
                 const void *A, const int lda, const void *X,
                 const int incX, const void *beta, void *Y, const int incY);
void blas_ztrmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const void *A, const int lda,
                 void *X, const int incX);
void blas_ztbmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void blas_ztpmv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);
void blas_ztrsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const void *A, const int lda, void *X,
                 const int incX);
void blas_ztbsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void blas_ztpsv(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE TransA, const enum BLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);


/*
 * Routines with S and D prefixes only
 */
void blas_ssymv(const enum BLAS_UPLO Uplo,
                 const int N, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY);
void blas_ssbmv(const enum BLAS_UPLO Uplo,
                 const int N, const int K, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY);
void blas_sspmv(const enum BLAS_UPLO Uplo,
                 const int N, const float alpha, const float *Ap,
                 const float *X, const int incX,
                 const float beta, float *Y, const int incY);
void blas_sger(const int M, const int N,
                const float alpha, const float *X, const int incX,
                const float *Y, const int incY, float *A, const int lda);
void blas_ssyr(const enum BLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, float *A, const int lda);
void blas_sspr(const enum BLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, float *Ap);
void blas_ssyr2(const enum BLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, const float *Y, const int incY, float *A,
                const int lda);
void blas_sspr2(const enum BLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, const float *Y, const int incY, float *A);

void blas_dsymv(const enum BLAS_UPLO Uplo,
                 const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
void blas_dsbmv(const enum BLAS_UPLO Uplo,
                 const int N, const int K, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
void blas_dspmv(const enum BLAS_UPLO Uplo,
                 const int N, const double alpha, const double *Ap,
                 const double *X, const int incX,
                 const double beta, double *Y, const int incY);
void blas_dger(const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda);
void blas_dsyr(const enum BLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *A, const int lda);
void blas_dspr(const enum BLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *Ap);
void blas_dsyr2(const enum BLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, const double *Y, const int incY, double *A,
                const int lda);
void blas_dspr2(const enum BLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, const double *Y, const int incY, double *A);


/*
 * Routines with C and Z prefixes only
 */
void blas_chemv(const enum BLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void blas_chbmv(const enum BLAS_UPLO Uplo,
                 const int N, const int K, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void blas_chpmv(const enum BLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *Ap,
                 const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void blas_cgeru(const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void blas_cgerc(const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void blas_cher(const enum BLAS_UPLO Uplo,
                const int N, const float alpha, const void *X, const int incX,
                void *A, const int lda);
void blas_chpr(const enum BLAS_UPLO Uplo,
                const int N, const float alpha, const void *X,
                const int incX, void *A);
void blas_cher2(const enum BLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *A, const int lda);
void blas_chpr2(const enum BLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *Ap);

void blas_zhemv(const enum BLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void blas_zhbmv(const enum BLAS_UPLO Uplo,
                 const int N, const int K, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void blas_zhpmv(const enum BLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *Ap,
                 const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void blas_zgeru(const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void blas_zgerc(const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void blas_zher(const enum BLAS_UPLO Uplo,
                const int N, const double alpha, const void *X, const int incX,
                void *A, const int lda);
void blas_zhpr(const enum BLAS_UPLO Uplo,
                const int N, const double alpha, const void *X,
                const int incX, void *A);
void blas_zher2(const enum BLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *A, const int lda);
void blas_zhpr2(const enum BLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void blas_sgemm(const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc);
void blas_ssymm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *B, const int ldb, const float beta,
                 float *C, const int ldc);
void blas_ssyrk(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const float *A, const int lda,
                 const float beta, float *C, const int ldc);
void blas_ssyr2k(const enum BLAS_UPLO Uplo,
                  const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                  const float alpha, const float *A, const int lda,
                  const float *B, const int ldb, const float beta,
                  float *C, const int ldc);
void blas_strmm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);
void blas_strsm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);

void blas_dgemm(const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
void blas_dsymm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double beta,
                 double *C, const int ldc);
void blas_dsyrk(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);
void blas_dsyr2k(const enum BLAS_UPLO Uplo,
                  const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                  const double alpha, const double *A, const int lda,
                  const double *B, const int ldb, const double beta,
                  double *C, const int ldc);
void blas_dtrmm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);
void blas_dtrsm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);

void blas_cgemm(const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);
void blas_csymm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void blas_csyrk(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc);
void blas_csyr2k(const enum BLAS_UPLO Uplo,
                  const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc);
void blas_ctrmm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);
void blas_ctrsm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);

void blas_zgemm(const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);
void blas_zsymm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void blas_zsyrk(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc);
void blas_zsyr2k(const enum BLAS_UPLO Uplo,
                  const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc);
void blas_ztrmm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);
void blas_ztrsm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const enum BLAS_TRANSPOSE TransA,
                 const enum BLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);


/*
 * Routines with prefixes C and Z only
 */
void blas_chemm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void blas_cherk(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const void *A, const int lda,
                 const float beta, void *C, const int ldc);
void blas_cher2k(const enum BLAS_UPLO Uplo,
                  const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const float beta,
                  void *C, const int ldc);
void blas_zhemm(const enum BLAS_SIDE Side,
                 const enum BLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void blas_zherk(const enum BLAS_UPLO Uplo,
                 const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const void *A, const int lda,
                 const double beta, void *C, const int ldc);
void blas_zher2k(const enum BLAS_UPLO Uplo,
                  const enum BLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const double beta,
                  void *C, const int ldc);

#ifdef __cplusplus
}
#endif

#endif
