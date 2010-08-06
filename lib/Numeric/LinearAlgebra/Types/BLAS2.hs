{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Types.BLAS2
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Matrix-Vector operations.
--

module Numeric.LinearAlgebra.Types.BLAS2
    where
     
import Data.Complex 
import Foreign( Ptr, with )

import Numeric.LinearAlgebra.Types.Enums
import Numeric.LinearAlgebra.Types.CEnums
import Numeric.LinearAlgebra.Types.BLAS1
import Numeric.LinearAlgebra.Types.Double 
import Numeric.LinearAlgebra.Types.Zomplex
   
-- | Types with matrix-vector operations.
class (BLAS1 a) => BLAS2 a where
    gemv :: Trans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    gbmv :: Trans -> Int -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmv :: Uplo -> Trans -> Diag -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbmv :: Uplo -> Trans -> Diag -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsv :: Uplo -> Trans -> Diag -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbsv :: Uplo -> Trans -> Diag -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    hemv :: Uplo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hbmv :: Uplo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    gerc :: Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    geru :: Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her  :: Uplo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her2 :: Uplo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

instance BLAS2 Double where
    gemv t = dgemv (cblasTrans t)
    gbmv t = dgbmv (cblasTrans t)
    trmv u t d = dtrmv (cblasUplo u) (cblasTrans t) (cblasDiag d)
    tbmv u t d = dtbmv (cblasUplo u) (cblasTrans t) (cblasDiag d)
    trsv u t d = dtrsv (cblasUplo u) (cblasTrans t) (cblasDiag d)
    tbsv u t d = dtbsv (cblasUplo u) (cblasTrans t) (cblasDiag d)
    hemv u = dsymv (cblasUplo u)
    hbmv u = dsbmv (cblasUplo u)
    gerc = dger
    geru = dger    
    her  u = dsyr  (cblasUplo u)
    her2 u = dsyr2 (cblasUplo u)
    
instance BLAS2 (Complex Double) where
    gemv transA m n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgemv (cblasTrans transA) m n pAlpha pA ldA pX incX pBeta pY incY
    
    gbmv transA m n kl ku alpha pA ldA pX incX beta pY incY =
         with alpha $ \pAlpha -> with beta $ \pBeta ->
             zgbmv (cblasTrans transA) m n kl ku pAlpha pA ldA pX incX pBeta pY incY

    trmv u t d n pA ldA pX incX =
        ztrmv (cblasUplo u) (cblasTrans t) (cblasDiag d) n pA ldA pX incX
            
    tbmv u t d n k pA ldA pX incX =
        ztbmv (cblasUplo u) (cblasTrans t) (cblasDiag d) n k pA ldA pX incX

    trsv u t d n pA ldA pX incX =
        ztrsv (cblasUplo u) (cblasTrans t) (cblasDiag d) n pA ldA pX incX

    tbsv u t d n k pA ldA pX incX =
        ztbsv (cblasUplo u) (cblasTrans t) (cblasDiag d) n k pA ldA pX incX
    
    hemv uplo n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhemv (cblasUplo uplo) n pAlpha pA ldA pX incX pBeta pY incY
    
    hbmv uplo n k alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhbmv (cblasUplo uplo) n k pAlpha pA ldA pX incX pBeta pY incY

    gerc m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha ->
            zgerc m n pAlpha pX incX pY incY pA ldA

    geru m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha ->
            zgeru m n pAlpha pX incX pY incY pA ldA
            
    her uplo n alpha pX incX pA ldA = 
        with alpha $ \pAlpha -> 
            zher (cblasUplo uplo) n pAlpha pX incX pA ldA
    
    her2 uplo n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha ->
            zher2 (cblasUplo uplo) n pAlpha pX incX pY incY pA ldA
