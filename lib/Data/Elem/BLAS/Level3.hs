{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Elem.BLAS.Level3
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Elem.BLAS.Level3
    where
     
import Data.Complex 
import Foreign ( Ptr, with )   

import BLAS.C.Types
import Data.Elem.BLAS.Level2
import Data.Elem.BLAS.Double  
import Data.Elem.BLAS.Zomplex 
        
        
class (BLAS2 a) => BLAS3 a where
    gemm  :: CBLASTrans -> CBLASTrans -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    symm  :: CBLASSide -> CBLASUpLo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hemm  :: CBLASSide -> CBLASUpLo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmm  :: CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsm  :: CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    syrk  :: CBLASUpLo -> CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    syr2k :: CBLASUpLo -> CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    herk  :: CBLASUpLo -> CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    her2k :: CBLASUpLo -> CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    
    
instance BLAS3 Double where
    gemm  = dgemm
    symm  = dsymm
    hemm  = dsymm
    trmm  = dtrmm
    trsm  = dtrsm
    syrk  = dsyrk
    syr2k = dsyr2k
    herk  = dsyrk
    her2k = dsyr2k
    
    
instance BLAS3 (Complex Double) where
    gemm transA transB m n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgemm transA transB m n k pAlpha pA ldA pB ldB pBeta pC ldC
    
    symm side uplo m n alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsymm side uplo m n pAlpha pA ldA pB ldB pBeta pC ldC

    hemm side uplo m n alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zhemm side uplo m n pAlpha pA ldA pB ldB pBeta pC ldC
    
    trmm side uplo transA diag m n alpha pA ldA pB ldB =
        with alpha $ \pAlpha -> 
            ztrmm side uplo transA diag m n pAlpha pA ldA pB ldB
            
    trsm side uplo transA diag m n alpha pA ldA pB ldB =
        with alpha $ \pAlpha -> 
            ztrsm side uplo transA diag m n pAlpha pA ldA pB ldB
            
    syrk uplo transA n k alpha pA ldA beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsyrk uplo transA n k pAlpha pA ldA pBeta pC ldC
            
    syr2k uplo transA n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsyr2k uplo transA n k pAlpha pA ldA pB ldB pBeta pC ldC

    herk uplo transA n k alpha pA ldA beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zherk uplo transA n k pAlpha pA ldA pBeta pC ldC
            
    her2k uplo transA n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zher2k uplo transA n k pAlpha pA ldA pB ldB pBeta pC ldC
