{-# LANGUAGE ForeignFunctionInterface #-}
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

foreign import ccall unsafe "cblas.h cblas_zdotu_sub"
    zdotu_sub   :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "cblas.h cblas_zdotc_sub"
    zdotc_sub   :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> IO ()


foreign import ccall unsafe "cblas.h cblas_dznrm2"
    znrm2  :: Int -> Ptr (Complex Double) -> Int -> IO Double

foreign import ccall unsafe "cblas.h cblas_dzasum"
    zasum  :: Int -> Ptr (Complex Double) -> Int -> IO Double

foreign import ccall unsafe "cblas.h cblas_izamax"
    izamax :: Int -> Ptr (Complex Double) -> Int -> IO Int

foreign import ccall unsafe "cblas.h cblas_zscal"
    zscal  :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zswap"
    zswap  :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zcopy"
    zcopy  :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zaxpy"
    zaxpy  :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zrotg"
    zrotg  :: Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "cblas.h cblas_zdrot"
    zdrot :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Double -> Double -> IO ()


---------------------------- Level 2 Routines -------------------------------

foreign import ccall unsafe "cblas.h cblas_zgemv"
    zgemv :: CBLASOrder -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zgbmv"
    zgbmv :: CBLASOrder -> CBLASTrans -> Int -> Int -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_ztrmv"
    ztrmv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_ztbmv"
    ztbmv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()
                 
foreign import ccall unsafe "cblas.h cblas_ztrsv"
    ztrsv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_ztbsv"
    ztbsv :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()
    
foreign import ccall unsafe "cblas.h cblas_zhemv"
    zhemv :: CBLASOrder -> CBLASUpLo -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zhbmv"
    zhbmv :: CBLASOrder -> CBLASUpLo -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
    
foreign import ccall unsafe "cblas.h cblas_zgeru"
    zgeru  :: CBLASOrder -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zgerc"
    zgerc  :: CBLASOrder -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()
        
foreign import ccall unsafe "cblas.h cblas_zher"
    zher  :: CBLASOrder -> CBLASUpLo -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zher2"
    zher2 :: CBLASOrder -> CBLASUpLo -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()


---------------------------- Level 3 Routines -------------------------------

foreign import ccall unsafe "cblas.h cblas_zgemm"
    zgemm  :: CBLASOrder -> CBLASTrans -> CBLASTrans -> Int -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zsymm"
    zsymm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zhemm"
    zhemm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_ztrmm"
    ztrmm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_ztrsm"
    ztrsm  :: CBLASOrder -> CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zsyrk"
    zsyrk  :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
           
foreign import ccall unsafe "cblas.h cblas_zsyr2k"           
    zsyr2k :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "cblas.h cblas_zherk"
    zherk  :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
           
foreign import ccall unsafe "cblas.h cblas_zher2k"           
    zher2k :: CBLASOrder -> CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
    