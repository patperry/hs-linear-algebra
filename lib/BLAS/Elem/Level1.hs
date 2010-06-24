{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Elem.Level1
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector operations.
--

module BLAS.Elem.Level1
    where
     
import Foreign( Storable, Ptr, peek, with )
import Foreign.Storable.Complex()
import Data.Complex( Complex(..) )

import BLAS.Elem.Double  
import BLAS.Elem.Zomplex
        
-- | Types with vector-vector operations.
class (Storable a, Fractional a) => BLAS1 a where
    copy  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()    
    swap  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    dot   :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO a
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
    dot     = ddot
    {-# INLINE dot #-}
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
    dot n pX incX pY incY =
        with 0 $ \pDotc -> do
            zdotc_sub n pX incX pY incY pDotc
            peek pDotc
    {-# INLINE dot #-}

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
