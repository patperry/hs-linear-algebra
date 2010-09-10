{-# LANGUAGE  ForeignFunctionInterface #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.LAPACK.Zomplex
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Foreign.LAPACK.Zomplex
    where

import Data.Complex( Complex )
import Foreign( Ptr )
import Foreign.BLAS.Types

#include "f77_func-hsc.h"


foreign import ccall unsafe #f77_func zgeqrf
    zgeqrf :: Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double)
           -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double)
           -> Ptr LAInt -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zgelqf
    zgelqf :: Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double)
           -> Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double)
           -> Ptr LAInt -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func zlarfg
    zlarfg :: Ptr LAInt -> Ptr (Complex Double) -> Ptr (Complex Double)
           -> Ptr LAInt -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe #f77_func zpotrf
    zpotrf :: BLASUplo -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt
           -> Ptr LAInt -> IO ()
           
foreign import ccall unsafe #f77_func zpotrs
    zpotrs :: BLASUplo -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double)
           -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt -> Ptr LAInt
           -> IO ()

foreign import ccall unsafe #f77_func zpptrf
    zpptrf :: BLASUplo -> Ptr LAInt -> Ptr (Complex Double)
           -> Ptr LAInt -> IO ()
           
foreign import ccall unsafe #f77_func zpptrs
    zpptrs :: BLASUplo -> Ptr LAInt -> Ptr LAInt -> Ptr (Complex Double)
           -> Ptr (Complex Double) -> Ptr LAInt -> Ptr LAInt
           -> IO ()

foreign import ccall unsafe #f77_func zunmqr
    zunmqr :: BLASSide -> BLASTrans -> Ptr LAInt -> Ptr LAInt
           -> Ptr LAInt -> Ptr (Complex Double) -> Ptr LAInt
           -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr LAInt
           -> Ptr (Complex Double) -> Ptr LAInt -> Ptr LAInt
           -> IO ()

foreign import ccall unsafe #f77_func zunmlq
    zunmlq :: BLASSide -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr LAInt
           -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double)
           -> Ptr (Complex Double) -> Ptr LAInt -> Ptr (Complex Double)
           -> Ptr LAInt -> Ptr LAInt -> IO ()

