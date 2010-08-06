{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Elem.BLAS2
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Matrix-Vector operations.
--

module Numeric.LinearAlgebra.Elem.BLAS2
    where
     
import Data.Complex 
import Foreign( Ptr, with )

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.CTypes
import Numeric.LinearAlgebra.Elem.BLAS1
import Numeric.LinearAlgebra.Elem.Double 
import Numeric.LinearAlgebra.Elem.Zomplex
   
-- | Types with matrix-vector operations.
class (BLAS1 a) => BLAS2 a where
    gemv :: Trans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    gbmv :: Trans -> Int -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmv :: UpLo -> Trans -> Diag -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbmv :: UpLo -> Trans -> Diag -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsv :: UpLo -> Trans -> Diag -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbsv :: UpLo -> Trans -> Diag -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    hemv :: UpLo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hbmv :: UpLo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    gerc :: Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    geru :: Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her  :: UpLo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her2 :: UpLo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

instance BLAS2 Double where
    gemv t = dgemv (cblasTrans t)
    gbmv t = dgbmv (cblasTrans t)
    trmv u t d = dtrmv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    tbmv u t d = dtbmv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    trsv u t d = dtrsv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    tbsv u t d = dtbsv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    hemv u = dsymv (cblasUpLo u)
    hbmv u = dsbmv (cblasUpLo u)
    gerc = dger
    geru = dger    
    her  u = dsyr  (cblasUpLo u)
    her2 u = dsyr2 (cblasUpLo u)
    
instance BLAS2 (Complex Double) where
    gemv transA m n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgemv (cblasTrans transA) m n pAlpha pA ldA pX incX pBeta pY incY
    
    gbmv transA m n kl ku alpha pA ldA pX incX beta pY incY =
         with alpha $ \pAlpha -> with beta $ \pBeta ->
             zgbmv (cblasTrans transA) m n kl ku pAlpha pA ldA pX incX pBeta pY incY

    trmv u t d n pA ldA pX incX =
        ztrmv (cblasUpLo u) (cblasTrans t) (cblasDiag d) n pA ldA pX incX
            
    tbmv u t d n k pA ldA pX incX =
        ztbmv (cblasUpLo u) (cblasTrans t) (cblasDiag d) n k pA ldA pX incX

    trsv u t d n pA ldA pX incX =
        ztrsv (cblasUpLo u) (cblasTrans t) (cblasDiag d) n pA ldA pX incX

    tbsv u t d n k pA ldA pX incX =
        ztbsv (cblasUpLo u) (cblasTrans t) (cblasDiag d) n k pA ldA pX incX
    
    hemv uplo n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhemv (cblasUpLo uplo) n pAlpha pA ldA pX incX pBeta pY incY
    
    hbmv uplo n k alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhbmv (cblasUpLo uplo) n k pAlpha pA ldA pX incX pBeta pY incY

    gerc m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha ->
            zgerc m n pAlpha pX incX pY incY pA ldA

    geru m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha ->
            zgeru m n pAlpha pX incX pY incY pA ldA
            
    her uplo n alpha pX incX pA ldA = 
        with alpha $ \pAlpha -> 
            zher (cblasUpLo uplo) n pAlpha pX incX pA ldA
    
    her2 uplo n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha ->
            zher2 (cblasUpLo uplo) n pAlpha pX incX pY incY pA ldA
