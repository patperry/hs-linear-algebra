{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Elem.BLAS.Level2
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Elem.BLAS.Level2
    where
     
import Data.Complex 
import Foreign ( Ptr, with )   

import Data.Elem.BLAS.Types
import Data.Elem.BLAS.Level1
import Data.Elem.BLAS.Double 
import Data.Elem.BLAS.Zomplex
        
class (BLAS1 a) => BLAS2 a where
    gemv :: CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    gbmv :: CBLASTrans -> Int -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmv :: CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbmv :: CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsv :: CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbsv :: CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    hemv :: CBLASUpLo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hbmv :: CBLASUpLo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    geru :: Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    gerc :: Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her  :: CBLASUpLo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her2 :: CBLASUpLo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

instance BLAS2 Double where
    gemv = dgemv
    gbmv = dgbmv
    trmv = dtrmv
    tbmv = dtbmv
    trsv = dtrsv
    tbsv = dtbsv
    hemv = dsymv
    hbmv = dsbmv
    geru = dger
    gerc = dger
    her  = dsyr
    her2 = dsyr2
    
instance BLAS2 (Complex Double) where
    gemv transA m n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgemv transA m n pAlpha pA ldA pX incX pBeta pY incY
    
    gbmv transA m n kl ku alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgbmv transA m n kl ku pAlpha pA ldA pX incX pBeta pY incY

    trmv = ztrmv
    tbmv = ztbmv
    trsv = ztrsv
    tbsv = ztbsv 
    
    hemv uplo n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhemv uplo n pAlpha pA ldA pX incX pBeta pY incY
    
    hbmv uplo n k alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhbmv uplo n k pAlpha pA ldA pX incX pBeta pY incY

    geru m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha -> zgeru m n pAlpha pX incX pY incY pA ldA

    gerc m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha -> zgerc m n pAlpha pX incX pY incY pA ldA

    her uplo n alpha pX incX pA ldA = 
        with alpha $ \pAlpha -> zher uplo n pAlpha pX incX pA ldA
    
    her2 uplo n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha -> zher2 uplo n pAlpha pX incX pY incY pA ldA
