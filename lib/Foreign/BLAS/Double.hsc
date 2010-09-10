{-# LANGUAGE ForeignFunctionInterface #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.BLAS.Double
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Foreign.BLAS.Double 
    where
    
import Foreign
import Foreign.BLAS.Types

#include "config.h"
#include "f77_func-hsc.h"


---------------------------- Level 1 Routines -------------------------------

foreign import ccall unsafe #f77_func ddot
    ddot :: Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO Double

foreign import ccall unsafe #f77_func dnrm2
    dnrm2  :: Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO Double

foreign import ccall unsafe #f77_func dasum
    dasum  :: Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO Double

foreign import ccall unsafe #f77_func idamax
    idamax :: Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO LAInt

foreign import ccall unsafe #f77_func dscal
    dscal  :: Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dswap
    dswap  :: Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dcopy
    dcopy  :: Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func daxpy
    daxpy  :: Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func drotg
    drotg  :: Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe #f77_func drot
    drot :: Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe #f77_func drotmg
    drotmg :: Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe #f77_func drotm
    drotm :: Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> IO ()


---------------------------- Level 2 Routines -------------------------------

foreign import ccall unsafe #f77_func dgemv
    dgemv :: BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dgbmv
    dgbmv ::  BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dtrmv
    dtrmv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dtbmv
    dtbmv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()
                 
foreign import ccall unsafe #f77_func dtrsv
    dtrsv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dtbsv
    dtbsv ::  BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()
    
foreign import ccall unsafe #f77_func dsymv
    dsymv ::  BLASUplo -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()
    
foreign import ccall unsafe #f77_func dsbmv
    dsbmv ::  BLASUplo -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()
    
foreign import ccall unsafe #f77_func dger
    dger  ::  Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()
        
foreign import ccall unsafe #f77_func dsyr
    dsyr  ::  BLASUplo -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dsyr2
    dsyr2 ::  BLASUplo -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dspmv
    dspmv ::  BLASUplo -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dspr
    dspr  ::  BLASUplo -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> IO ()

foreign import ccall unsafe #f77_func dspr2
    dspr2 ::  BLASUplo -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> IO ()


---------------------------- Level 3 Routines -------------------------------

foreign import ccall unsafe #f77_func dgemm
    dgemm  ::  BLASTrans -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dsymm
    dsymm  ::  BLASSide -> BLASUplo -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dtrmm
    dtrmm  ::  BLASSide -> BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dtrsm
    dtrsm  ::  BLASSide -> BLASUplo -> BLASTrans -> BLASDiag -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dsyrk
    dsyrk  ::  BLASUplo -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()
           
foreign import ccall unsafe #f77_func dsyr2k           
    dsyr2k ::  BLASUplo -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> IO ()
