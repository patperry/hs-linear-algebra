{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.C.Level2
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.C.Level2
    where
     
import Data.Complex 
import Foreign ( Ptr, with )   

import BLAS.C.Types
import BLAS.C.Level1
import BLAS.C.Double 
import BLAS.C.Zomplex
        
class (BLAS1 a) => BLAS2 a where
    gemv :: CBLASOrder -> CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    gbmv :: CBLASOrder -> CBLASTrans -> Int -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbmv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    tbsv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    hemv :: CBLASOrder -> CBLASUpLo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hbmv :: CBLASOrder -> CBLASUpLo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    geru :: CBLASOrder -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    gerc :: CBLASOrder -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her  :: CBLASOrder -> CBLASUpLo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    her2 :: CBLASOrder -> CBLASUpLo -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

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
    gemv order transA m n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgemv order transA m n pAlpha pA ldA pX incX pBeta pY incY
    
    gbmv order transA m n kl ku alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgbmv order transA m n kl ku pAlpha pA ldA pX incX pBeta pY incY

    trmv = ztrmv
    tbmv = ztbmv
    trsv = ztrsv
    tbsv = ztbsv 
    
    hemv order uplo n alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhemv order uplo n pAlpha pA ldA pX incX pBeta pY incY
    
    hbmv order uplo n k alpha pA ldA pX incX beta pY incY =
        with alpha $ \pAlpha -> with beta $ \pBeta -> 
            zhbmv order uplo n k pAlpha pA ldA pX incX pBeta pY incY

    geru order m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha -> zgeru order m n pAlpha pX incX pY incY pA ldA

    gerc order m n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha -> zgerc order m n pAlpha pX incX pY incY pA ldA

    her order uplo n alpha pX incX pA ldA = 
        with alpha $ \pAlpha -> zher order uplo n pAlpha pX incX pA ldA
    
    her2 order uplo n alpha pX incX pY incY pA ldA = 
        with alpha $ \pAlpha -> zher2 order uplo n pAlpha pX incX pY incY pA ldA
