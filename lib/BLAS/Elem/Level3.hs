{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Elem.Level3
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Matrix-Matrix operations.
--

module BLAS.Elem.Level3
    where
     
import Data.Complex 
import Foreign ( Ptr, with )   

import BLAS.Types
import BLAS.CTypes
import BLAS.Elem.Level2
import BLAS.Elem.Double  
import BLAS.Elem.Zomplex 
        
-- | Types with matrix-matrix operations.        
class (BLAS2 a) => BLAS3 a where
    gemm  :: TransEnum -> TransEnum -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    symm  :: SideEnum -> UpLoEnum -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hemm  :: SideEnum -> UpLoEnum -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmm  :: SideEnum -> UpLoEnum -> TransEnum -> DiagEnum -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsm  :: SideEnum -> UpLoEnum -> TransEnum -> DiagEnum -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    syrk  :: UpLoEnum -> TransEnum -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    syr2k :: UpLoEnum -> TransEnum -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    herk  :: UpLoEnum -> TransEnum -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    her2k :: UpLoEnum -> TransEnum -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    
    
instance BLAS3 Double where
    gemm ta tb = dgemm (cblasTrans ta) (cblasTrans tb)
    symm  s u = dsymm (cblasSide s) (cblasUpLo u) 
    hemm  s u = dsymm (cblasSide s) (cblasUpLo u) 
    trmm  s u t d = dtrmm (cblasSide s) (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    trsm  s u t d = dtrsm (cblasSide s) (cblasUpLo u) (cblasTrans t) (cblasDiag d)
    syrk  u t = dsyrk  (cblasUpLo u) (cblasTrans t)
    syr2k u t = dsyr2k (cblasUpLo u) (cblasTrans t)
    herk  u t = dsyrk  (cblasUpLo u) (cblasTrans t)
    her2k u t = dsyr2k (cblasUpLo u) (cblasTrans t)
    
    
instance BLAS3 (Complex Double) where
    gemm transA transB m n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgemm (cblasTrans transA) (cblasTrans transB) m n k pAlpha pA ldA pB ldB pBeta pC ldC
    
    symm side uplo m n alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsymm (cblasSide side) (cblasUpLo uplo) m n pAlpha pA ldA pB ldB pBeta pC ldC

    hemm side uplo m n alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zhemm (cblasSide side) (cblasUpLo uplo) m n pAlpha pA ldA pB ldB pBeta pC ldC
    
    trmm side uplo transA diag m n alpha pA ldA pB ldB =
        with alpha $ \pAlpha -> 
            ztrmm (cblasSide side) (cblasUpLo uplo) (cblasTrans transA) (cblasDiag diag) m n pAlpha pA ldA pB ldB
            
    trsm side uplo transA diag m n alpha pA ldA pB ldB =
        with alpha $ \pAlpha -> 
            ztrsm (cblasSide side) (cblasUpLo uplo) (cblasTrans transA) (cblasDiag diag) m n pAlpha pA ldA pB ldB
            
    syrk uplo transA n k alpha pA ldA beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsyrk (cblasUpLo uplo) (cblasTrans transA) n k pAlpha pA ldA pBeta pC ldC
            
    syr2k uplo transA n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsyr2k (cblasUpLo uplo) (cblasTrans transA) n k pAlpha pA ldA pB ldB pBeta pC ldC

    herk uplo transA n k alpha pA ldA beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zherk (cblasUpLo uplo) (cblasTrans transA) n k pAlpha pA ldA pBeta pC ldC
            
    her2k uplo transA n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zher2k (cblasUpLo uplo) (cblasTrans transA) n k pAlpha pA ldA pB ldB pBeta pC ldC
