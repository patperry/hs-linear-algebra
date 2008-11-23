{-# LANGUAGE ForeignFunctionInterface #-}
{-# CFILES cbits/zomplex.c #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.C.Zomplex
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.C.Zomplex
    where
        
import Data.Complex ( Complex )
import Foreign.Ptr  ( Ptr )
import BLAS.C.Types

---------------------------- Level 1 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_zdotu_sub"
    zdotu_sub   :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "BLAS.h blas_zdotc_sub"
    zdotc_sub   :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> IO ()


foreign import ccall unsafe "BLAS.h blas_dznrm2"
    znrm2  :: Int -> Ptr (Complex Double) -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_dzasum"
    zasum  :: Int -> Ptr (Complex Double) -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_izamax"
    izamax :: Int -> Ptr (Complex Double) -> Int -> IO Int

foreign import ccall unsafe "BLAS.h blas_zscal"
    zscal  :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zswap"
    zswap  :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zcopy"
    zcopy  :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zaxpy"
    zaxpy  :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zrotg"
    zrotg  :: Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "BLAS.h blas_zdrot"
    zdrot :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Double -> Double -> IO ()


---------------------------- Level 2 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_zgemv"
    zgemv :: CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zgbmv"
    zgbmv ::  CBLASTrans -> Int -> Int -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztrmv"
    ztrmv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztbmv"
    ztbmv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()
                 
foreign import ccall unsafe "BLAS.h blas_ztrsv"
    ztrsv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztbsv"
    ztbsv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_zhemv"
    zhemv ::  CBLASUpLo -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zhbmv"
    zhbmv ::  CBLASUpLo -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_zgeru"
    zgeru  ::  Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zgerc"
    zgerc  ::  Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()
        
foreign import ccall unsafe "BLAS.h blas_zher"
    zher  ::  CBLASUpLo -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zher2"
    zher2 ::  CBLASUpLo -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()


---------------------------- Level 3 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_zgemm"
    zgemm  ::  CBLASTrans -> CBLASTrans -> Int -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zsymm"
    zsymm  ::  CBLASSide -> CBLASUpLo -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zhemm"
    zhemm  ::  CBLASSide -> CBLASUpLo -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztrmm"
    ztrmm  ::  CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztrsm"
    ztrsm  ::  CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zsyrk"
    zsyrk  ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
           
foreign import ccall unsafe "BLAS.h blas_zsyr2k"           
    zsyr2k ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zherk"
    zherk  ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
           
foreign import ccall unsafe "BLAS.h blas_zher2k"           
    zher2k ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
    