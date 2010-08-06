{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Types.BLAS1
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector operations.
--

module Numeric.LinearAlgebra.Types.BLAS1
    where
     
import Foreign( Storable, Ptr, peek, with )
import Foreign.Storable.Complex()
import Data.Complex( Complex(..) )

import Numeric.LinearAlgebra.Types.Double  
import Numeric.LinearAlgebra.Types.Zomplex
        
-- | Types with vector-vector operations.
class (Storable a, Fractional a) => BLAS1 a where
    copy  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()    
    swap  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    dotc  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO a
    dotu  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO a    
    nrm2  :: Int -> Ptr a -> Int -> IO Double
    asum  :: Int -> Ptr a -> Int -> IO Double
    iamax :: Int -> Ptr a -> Int -> IO Int
    
    scal  :: Int -> a -> Ptr a -> Int -> IO () 

    axpy  :: Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    rotg  :: Ptr a -> Ptr a -> Ptr a -> Ptr a -> IO ()
    rot   :: Int -> Ptr a -> Int -> Ptr a -> Int -> Double -> Double -> IO ()


instance BLAS1 Double where
    copy = dcopy
    {-# INLINE copy #-}
    swap = dswap
    {-# INLINE swap #-}
    dotc = ddot
    {-# INLINE dotc #-}
    dotu = ddot
    {-# INLINE dotu #-}
    nrm2    = dnrm2
    {-# INLINE nrm2 #-}
    asum    = dasum
    {-# INLINE asum #-}
    iamax   = idamax
    {-# INLINE iamax #-}
    axpy    = daxpy
    {-# INLINE axpy #-}
    scal    = dscal
    {-# INLINE scal #-}
    rotg    = drotg
    {-# INLINE rotg #-}
    rot     = drot
    {-# INLINE rot #-}


instance BLAS1 (Complex Double) where
    dotc n pX incX pY incY =
        with 0 $ \pDotc -> do
            zdotc_sub n pX incX pY incY pDotc
            peek pDotc
    {-# INLINE dotc #-}
    dotu n pX incX pY incY =
        with 0 $ \pDotu -> do
            zdotu_sub n pX incX pY incY pDotu
            peek pDotu
    {-# INLINE dotu #-}

    copy = zcopy
    {-# INLINE copy #-}
    swap = zswap
    {-# INLINE swap #-}
    nrm2  = znrm2
    {-# INLINE nrm2 #-}
    asum  = zasum
    {-# INLINE asum #-}
    iamax = izamax
    {-# INLINE iamax #-}
    
    axpy n alpha pX incX pY incY = 
        with alpha $ \pAlpha ->
            zaxpy n pAlpha pX incX pY incY
    {-# INLINE axpy #-}
    
    scal n alpha pX incX =
        with alpha $ \pAlpha ->
            zscal n pAlpha pX incX
    {-# INLINE scal #-}
            
    rotg  = zrotg
    {-# INLINE rotg #-}
    rot = zdrot
    {-# INLINE rot #-}
