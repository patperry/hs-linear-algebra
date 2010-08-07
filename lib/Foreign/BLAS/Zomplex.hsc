{-# LANGUAGE ForeignFunctionInterface #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.BLAS.Zomplex
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Foreign.BLAS.Zomplex 
    where

import Data.Complex( Complex )
import Foreign
import Foreign.BLAS.Types

#include "config.h"
#include "f77_func-hsc.h"
#define la_int int

---------------------------- Level 1 Routines -------------------------------

foreign import ccall unsafe #f77_func zdotu
    zdotu :: Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zdotc
    zdotc :: Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()


foreign import ccall unsafe #f77_func dznrm2
    znrm2  :: Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO Double

foreign import ccall unsafe #f77_func dzasum
    zasum  :: Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO Double

foreign import ccall unsafe #f77_func izamax
    izamax :: Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO (#type la_int)

foreign import ccall unsafe #f77_func zscal
    zscal  :: Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zswap
    zswap  :: Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zcopy
    zcopy  :: Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zaxpy
    zaxpy  :: Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zrotg
    zrotg  :: Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe #f77_func zdrot
    zdrot :: Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr Double -> Ptr Double -> IO ()


---------------------------- Level 2 Routines -------------------------------

foreign import ccall unsafe #f77_func zgemv
    zgemv :: BLASTrans -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zgbmv
    zgbmv ::  BLASTrans -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func ztrmv
    ztrmv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func ztbmv
    ztbmv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()
                 
foreign import ccall unsafe #f77_func ztrsv
    ztrsv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func ztbsv
    ztbsv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()
    
foreign import ccall unsafe #f77_func zhemv
    zhemv ::  BLASUplo -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zhbmv
    zhbmv ::  BLASUplo -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()
    
foreign import ccall unsafe #f77_func zgeru
    zgeru  ::  Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zgerc
    zgerc  ::  Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()
        
foreign import ccall unsafe #f77_func zher
    zher  ::  BLASUplo -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zher2
    zher2 ::  BLASUplo -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()


---------------------------- Level 3 Routines -------------------------------

foreign import ccall unsafe #f77_func zgemm
    zgemm  ::  BLASTrans -> BLASTrans -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zsymm
    zsymm  ::  BLASSide -> BLASUplo -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zhemm
    zhemm  ::  BLASSide -> BLASUplo -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func ztrmm
    ztrmm  ::  BLASSide -> BLASUplo -> BLASTrans -> BLASDiag -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func ztrsm
    ztrsm  ::  BLASSide -> BLASUplo -> BLASTrans -> BLASDiag -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zsyrk
    zsyrk  ::  BLASUplo -> BLASTrans -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()
           
foreign import ccall unsafe #f77_func zsyr2k           
    zsyr2k ::  BLASUplo -> BLASTrans -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()

foreign import ccall unsafe #f77_func zherk
    zherk  ::  BLASUplo -> BLASTrans -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()
           
foreign import ccall unsafe #f77_func zher2k           
    zher2k ::  BLASUplo -> BLASTrans -> Ptr (#type la_int) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (#type la_int) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (#type la_int) -> IO ()
