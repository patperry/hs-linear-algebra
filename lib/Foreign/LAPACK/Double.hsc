{-# LANGUAGE ForeignFunctionInterface #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.LAPACK.Double
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Foreign.LAPACK.Double
    where

import Foreign( Ptr )
import Foreign.BLAS.Types

#include "f77_func-hsc.h"


foreign import ccall unsafe #f77_func dgeqrf
    dgeqrf :: Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double
           -> Ptr Double -> Ptr LAInt -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dgelqf
    dgelqf :: Ptr LAInt -> Ptr LAInt -> Ptr Double -> Ptr LAInt -> Ptr Double
           -> Ptr Double -> Ptr LAInt -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dormqr
    dormqr :: BLASSide -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr LAInt
           -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt
           -> Ptr Double -> Ptr LAInt -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dormlq
    dormlq :: BLASSide -> BLASTrans -> Ptr LAInt -> Ptr LAInt -> Ptr LAInt
           -> Ptr Double -> Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt
           -> Ptr Double -> Ptr LAInt -> Ptr LAInt -> IO ()

foreign import ccall unsafe #f77_func dlarfg
    dlarfg :: Ptr LAInt -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr Double
           -> IO ()
