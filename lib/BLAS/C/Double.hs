{-# LANGUAGE ForeignFunctionInterface #-}
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

foreign import ccall unsafe "cblas.h cblas_ddot"
    ddot   :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO Double

foreign import ccall unsafe "cblas.h cblas_dnrm2"
    dnrm2  :: Int -> Ptr Double -> Int -> IO Double

foreign import ccall unsafe "cblas.h cblas_dasum"
    dasum  :: Int -> Ptr Double -> Int -> IO Double

foreign import ccall unsafe "cblas.h cblas_idamax"
    idamax :: Int -> Ptr Double -> Int -> IO Int

foreign import ccall unsafe "cblas.h cblas_dscal"
    dscal  :: Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dswap"
    dswap  :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dcopy"
    dcopy  :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_daxpy"
    daxpy  :: Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_drotg"
    drotg  :: Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "cblas.h cblas_drot"
    drot :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Double -> IO ()

foreign import ccall unsafe "cblas.h cblas_drotmg"
    drotmg :: Ptr Double -> Ptr Double -> Ptr Double -> Double -> Ptr Double -> IO ()

foreign import ccall unsafe "cblas.h cblas_drotm"
    drotm :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> Ptr Double -> IO ()


---------------------------- Level 2 Routines -------------------------------

foreign import ccall unsafe "cblas.h cblas_dgemv"
    dgemv :: CBLASOrder -> CBLASTrans -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dgbmv"
    dgbmv :: CBLASOrder -> CBLASTrans -> Int -> Int -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dtrmv"
    dtrmv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dtbmv"
    dtbmv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()
                 
foreign import ccall unsafe "cblas.h cblas_dtrsv"
    dtrsv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dtbsv"
    dtbsv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()
    
foreign import ccall unsafe "cblas.h cblas_dsymv"
    dsymv :: CBLASOrder -> CBLASUpLo -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
    
foreign import ccall unsafe "cblas.h cblas_dsbmv"
    dsbmv :: CBLASOrder -> CBLASUpLo -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
    
foreign import ccall unsafe "cblas.h cblas_dger"
    dger  :: CBLASOrder -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()
        
foreign import ccall unsafe "cblas.h cblas_dsyr"
    dsyr  :: CBLASOrder -> CBLASUpLo -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dsyr2"
    dsyr2 :: CBLASOrder -> CBLASUpLo -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()


---------------------------- Level 3 Routines -------------------------------

foreign import ccall unsafe "cblas.h cblas_dgemm"
    dgemm  :: CBLASOrder -> CBLASTrans -> CBLASTrans -> Int -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dsymm"
    dsymm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dtrmm"
    dtrmm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dtrsm"
    dtrsm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_dsyrk"
    dsyrk  :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> Double -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
           
foreign import ccall unsafe "cblas.h cblas_dsyr2k"           
    dsyr2k :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
    