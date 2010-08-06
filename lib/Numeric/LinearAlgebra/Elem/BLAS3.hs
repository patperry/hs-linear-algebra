{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Elem.BLAS3
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Matrix-Matrix operations.
--

module Numeric.LinearAlgebra.Elem.BLAS3
    where
     
import Data.Complex 
import Foreign ( Ptr, with )   

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.CTypes
import Numeric.LinearAlgebra.Elem.BLAS2
import Numeric.LinearAlgebra.Elem.Double  
import Numeric.LinearAlgebra.Elem.Zomplex 
        
-- | Types with matrix-matrix operations.        
class (BLAS2 a) => BLAS3 a where
    gemm  :: Trans -> Trans -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    symm  :: Side -> Uplo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hemm  :: Side -> Uplo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmm  :: Side -> Uplo -> Trans -> Diag -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsm  :: Side -> Uplo -> Trans -> Diag -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    syrk  :: Uplo -> Trans -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    syr2k :: Uplo -> Trans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    herk  :: Uplo -> Trans -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    her2k :: Uplo -> Trans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    
    
instance BLAS3 Double where
    gemm ta tb = dgemm (cblasTrans ta) (cblasTrans tb)
    symm  s u = dsymm (cblasSide s) (cblasUplo u) 
    hemm  s u = dsymm (cblasSide s) (cblasUplo u) 
    trmm  s u t d = dtrmm (cblasSide s) (cblasUplo u) (cblasTrans t) (cblasDiag d)
    trsm  s u t d = dtrsm (cblasSide s) (cblasUplo u) (cblasTrans t) (cblasDiag d)
    syrk  u t = dsyrk  (cblasUplo u) (cblasTrans t)
    syr2k u t = dsyr2k (cblasUplo u) (cblasTrans t)
    herk  u t = dsyrk  (cblasUplo u) (cblasTrans t)
    her2k u t = dsyr2k (cblasUplo u) (cblasTrans t)
    
    
instance BLAS3 (Complex Double) where
    gemm transA transB m n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zgemm (cblasTrans transA) (cblasTrans transB) m n k pAlpha pA ldA pB ldB pBeta pC ldC
    
    symm side uplo m n alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsymm (cblasSide side) (cblasUplo uplo) m n pAlpha pA ldA pB ldB pBeta pC ldC

    hemm side uplo m n alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zhemm (cblasSide side) (cblasUplo uplo) m n pAlpha pA ldA pB ldB pBeta pC ldC
    
    trmm side uplo transA diag m n alpha pA ldA pB ldB =
        with alpha $ \pAlpha -> 
            ztrmm (cblasSide side) (cblasUplo uplo) (cblasTrans transA) (cblasDiag diag) m n pAlpha pA ldA pB ldB
            
    trsm side uplo transA diag m n alpha pA ldA pB ldB =
        with alpha $ \pAlpha -> 
            ztrsm (cblasSide side) (cblasUplo uplo) (cblasTrans transA) (cblasDiag diag) m n pAlpha pA ldA pB ldB
            
    syrk uplo transA n k alpha pA ldA beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsyrk (cblasUplo uplo) (cblasTrans transA) n k pAlpha pA ldA pBeta pC ldC
            
    syr2k uplo transA n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zsyr2k (cblasUplo uplo) (cblasTrans transA) n k pAlpha pA ldA pB ldB pBeta pC ldC

    herk uplo transA n k alpha pA ldA beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zherk (cblasUplo uplo) (cblasTrans transA) n k pAlpha pA ldA pBeta pC ldC
            
    her2k uplo transA n k alpha pA ldA pB ldB beta pC ldC =
        with alpha $ \pAlpha -> with beta $ \pBeta ->
            zher2k (cblasUplo uplo) (cblasTrans transA) n k pAlpha pA ldA pB ldB pBeta pC ldC
