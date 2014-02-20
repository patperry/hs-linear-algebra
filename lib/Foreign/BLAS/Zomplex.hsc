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
import Foreign.C.Types

#include "config.h"
#include "f77_func-hsc.h"
#define la_int int

---------------------------- Level 1 Routines -------------------------------

foreign import ccall unsafe #f77_func zdotu
    zdotu :: Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zdotc
    zdotc :: Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()


foreign import ccall unsafe #f77_func dznrm2
    znrm2  :: Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO Double

foreign import ccall unsafe #f77_func dzasum
    zasum  :: Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO Double

foreign import ccall unsafe #f77_func izamax
    izamax_hidden :: Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO CInt 
izamax :: Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO LAInt    
izamax a b c = do 
                res <- izamax_hidden a b c 
                return $! LAInt res 


foreign import ccall unsafe #f77_func zscal
    zscal  :: Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zswap
    zswap  :: Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zcopy
    zcopy  :: Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zaxpy
    zaxpy  :: Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zrotg
    zrotg  :: Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe #f77_func zdrot
    zdrot :: Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr Double -> Ptr Double -> IO ()


---------------------------- Level 2 Routines -------------------------------

foreign import ccall unsafe #f77_func zgemv
    zgemv :: BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zgbmv
    zgbmv ::  BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func ztrmv
    ztrmv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func ztpmv
    ztpmv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func ztpsv
    ztpsv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func ztbmv
    ztbmv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()
                 
foreign import ccall unsafe #f77_func ztrsv
    ztrsv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func ztbsv
    ztbsv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()
    
foreign import ccall unsafe #f77_func zhemv
    zhemv ::  BLASUplo -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zhbmv
    zhbmv ::  BLASUplo -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()
    
foreign import ccall unsafe #f77_func zgeru
    zgeru  ::  Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zgerc
    zgerc  ::  Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()
        
foreign import ccall unsafe #f77_func zher
    zher  ::  BLASUplo -> Ptr LAInt -> Ptr Double -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zher2
    zher2 ::  BLASUplo -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zhpmv
    zhpmv ::  BLASUplo -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zhpr
    zhpr  ::  BLASUplo -> Ptr LAInt -> Ptr Double -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe #f77_func zhpr2
    zhpr2 ::  BLASUplo -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> IO ()


---------------------------- Level 3 Routines -------------------------------

foreign import ccall unsafe #f77_func zgemm
    zgemm  ::  BLASTrans -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zsymm
    zsymm  ::  BLASSide -> BLASUplo -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zhemm
    zhemm  ::  BLASSide -> BLASUplo -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func ztrmm
    ztrmm  ::  BLASSide -> BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func ztrsm
    ztrsm  ::  BLASSide -> BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zsyrk
    zsyrk  ::  BLASUplo -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()
           
foreign import ccall unsafe #f77_func zsyr2k           
    zsyr2k ::  BLASUplo -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zherk
    zherk  ::  BLASUplo -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()
           
foreign import ccall unsafe #f77_func zher2k           
    zher2k ::  BLASUplo -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt -> IO ()
