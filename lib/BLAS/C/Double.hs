{-# LANGUAGE ForeignFunctionInterface #-}
{-# CFILES cbits/double.c #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.C.Double
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.C.Double 
    where
        
import Foreign.Ptr ( Ptr )
import BLAS.C.Types

---------------------------- Level 1 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_ddot"
    ddot :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_dnrm2"
    dnrm2  :: Int -> Ptr Double -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_dasum"
    dasum  :: Int -> Ptr Double -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_idamax"
    idamax :: Int -> Ptr Double -> Int -> IO Int

foreign import ccall unsafe "BLAS.h blas_dscal"
    dscal  :: Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dswap"
    dswap  :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dcopy"
    dcopy  :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_daxpy"
    daxpy  :: Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_drotg"
    drotg  :: Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "BLAS.h blas_drot"
    drot :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Double -> IO ()

foreign import ccall unsafe "BLAS.h blas_drotmg"
    drotmg :: Ptr Double -> Ptr Double -> Ptr Double -> Double -> Ptr Double -> IO ()

foreign import ccall unsafe "BLAS.h blas_drotm"
    drotm :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> Ptr Double -> IO ()


---------------------------- Level 2 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_dgemv"
    dgemv :: CBLASTrans -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dgbmv"
    dgbmv ::  CBLASTrans -> Int -> Int -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtrmv"
    dtrmv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtbmv"
    dtbmv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()
                 
foreign import ccall unsafe "BLAS.h blas_dtrsv"
    dtrsv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtbsv"
    dtbsv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_dsymv"
    dsymv ::  CBLASUpLo -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_dsbmv"
    dsbmv ::  CBLASUpLo -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_dger"
    dger  ::  Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()
        
foreign import ccall unsafe "BLAS.h blas_dsyr"
    dsyr  ::  CBLASUpLo -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dsyr2"
    dsyr2 ::  CBLASUpLo -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()


---------------------------- Level 3 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_dgemm"
    dgemm  ::  CBLASTrans -> CBLASTrans -> Int -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dsymm"
    dsymm  ::  CBLASSide -> CBLASUpLo -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtrmm"
    dtrmm  ::  CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtrsm"
    dtrsm  ::  CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dsyrk"
    dsyrk  ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Double -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
           
foreign import ccall unsafe "BLAS.h blas_dsyr2k"           
    dsyr2k ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
    