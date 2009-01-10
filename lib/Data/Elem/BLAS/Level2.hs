{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Elem.BLAS.Level2
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Matrix-Vector operations.
--

module Data.Elem.BLAS.Level2
    where
     
import Data.Complex 
import Foreign ( Ptr, with )   

import BLAS.Types
import BLAS.CTypes
import Data.Elem.BLAS.Level1
import Data.Elem.BLAS.Double 
import Data.Elem.BLAS.Zomplex
   
-- | Types with matrix-vector operations.
class (BLAS1 a) => BLAS2 a where
    gemv :: TransEnum -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    gbmv :: TransEnum -> Int -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmv :: UpLoEnum -> TransEnum -> DiagEnum -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbmv :: UpLoEnum -> TransEnum -> DiagEnum -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsv :: UpLoEnum -> TransEnum -> DiagEnum -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbsv :: UpLoEnum -> TransEnum -> DiagEnum -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    hemv :: UpLoEnum -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hbmv :: UpLoEnum -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    geru :: Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    gerc :: Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her  :: UpLoEnum -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her2 :: UpLoEnum -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

instance BLAS2 Double where
    gemv t = dgemv (cblasTrans t)
    gbmv t = dgbmv (cblasTrans t)
    trmv u t d = dtrmv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    tbmv u t d = dtbmv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    trsv u t d = dtrsv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    tbsv u t d = dtbsv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    hemv u = dsymv (cblasUpLo u)
    hbmv u = dsbmv (cblasUpLo u)
    geru = dger
    gerc = dger
    her  u = dsyr  (cblasUpLo u)
    her2 u = dsyr2 (cblasUpLo u)
    
instance BLAS2 (Complex Double) where
    gemv transA m n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgemv (cblasTrans transA) m n pAlpha pA ldA pX incX pBeta pY incY
    
    gbmv transA m n kl ku alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgbmv (cblasTrans transA) m n kl ku pAlpha pA ldA pX incX pBeta pY incY

    trmv u t d = ztrmv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    tbmv u t d = ztbmv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    trsv u t d = ztrsv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    tbsv u t d = ztbsv (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    
    hemv uplo n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhemv (cblasUpLo uplo) n pAlpha pA ldA pX incX pBeta pY incY
    
    hbmv uplo n k alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhbmv (cblasUpLo uplo) n k pAlpha pA ldA pX incX pBeta pY incY

    geru m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha -> 
            zgeru m n pAlpha pX incX pY incY pA ldA

    gerc m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha -> 
            zgerc m n pAlpha pX incX pY incY pA ldA

    her uplo n alpha pX incX pA ldA = 
        with alpha $ \pAlpha -> 
            zher (cblasUpLo uplo) n pAlpha pX incX pA ldA
    
    her2 uplo n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha ->
            zher2 (cblasUpLo uplo) n pAlpha pX incX pY incY pA ldA
