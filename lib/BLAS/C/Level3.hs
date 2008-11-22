{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.C.Level3
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.C.Level3
    where
     
import Data.Complex 
import Foreign ( Ptr, with )   

import BLAS.C.Types
import BLAS.C.Level2
import BLAS.C.Double  
import BLAS.C.Zomplex 
        
        
class (BLAS2 a) => BLAS3 a where
    gemm  :: CBLASOrder -> CBLASTrans -> CBLASTrans -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    symm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hemm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    syrk  :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    syr2k :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    herk  :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    her2k :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    
    
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
    gemm order transA transB m n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgemm order transA transB m n k pAlpha pA ldA pB ldB pBeta pC ldC
    
    symm order side uplo m n alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsymm order side uplo m n pAlpha pA ldA pB ldB pBeta pC ldC

    hemm order side uplo m n alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zhemm order side uplo m n pAlpha pA ldA pB ldB pBeta pC ldC
    
    trmm order side uplo transA diag m n alpha pA ldA pB ldB =
        with alpha $ \pAlpha -> 
            ztrmm order side uplo transA diag m n pAlpha pA ldA pB ldB
            
    trsm order side uplo transA diag m n alpha pA ldA pB ldB =
        with alpha $ \pAlpha -> 
            ztrsm order side uplo transA diag m n pAlpha pA ldA pB ldB
            
    syrk order uplo transA n k alpha pA ldA beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsyrk order uplo transA n k pAlpha pA ldA pBeta pC ldC
            
    syr2k order uplo transA n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsyr2k order uplo transA n k pAlpha pA ldA pB ldB pBeta pC ldC

    herk order uplo transA n k alpha pA ldA beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zherk order uplo transA n k pAlpha pA ldA pBeta pC ldC
            
    her2k order uplo transA n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zher2k order uplo transA n k pAlpha pA ldA pB ldB pBeta pC ldC
