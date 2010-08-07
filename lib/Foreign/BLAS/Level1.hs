{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.BLAS.Level1
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector operations.
--

module Foreign.BLAS.Level1 (
    BLAS1(..),
    ) where
     
import Foreign( Storable, Ptr, peek, with )
import Foreign.Storable.Complex()
import Data.Complex( Complex(..) )

import Foreign.BLAS.Double  
import Foreign.BLAS.Zomplex
        
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


withEnum :: (Enum a, Storable a) => Int -> (Ptr a -> IO b) -> IO b
withEnum = with . toEnum
{-# INLINE withEnum #-}

instance BLAS1 Double where
    copy n px incx py incy =
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
            dcopy pn px pincx py pincy
    {-# INLINE copy #-}

    swap n px incx py incy = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
            dswap pn px pincx py pincy
    {-# INLINE swap #-}

    dotc n px incx py incy = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
            ddot pn px pincx py pincy
    {-# INLINE dotc #-}

    dotu = dotc
    {-# INLINE dotu #-}
    
    nrm2 n px incx = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
            dnrm2 pn px pincx
    {-# INLINE nrm2 #-}

    asum n px incx = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
            dasum pn px pincx
    {-# INLINE asum #-}

    iamax n px incx = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx -> do
            i <- idamax pn px pincx
            return $! fromEnum (i - 1)
    {-# INLINE iamax #-}

    axpy n alpha px incx py incy = 
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
            daxpy pn palpha px pincx py pincy
    {-# INLINE axpy #-}
    
    scal n alpha px incx = 
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum incx $ \pincx ->
            dscal pn palpha px pincx
    {-# INLINE scal #-}

    rotg = drotg
    {-# INLINE rotg #-}
    
    rot n px incx py incy c s = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
        with c $ \pc ->
        with s $ \ps ->
            drot pn px pincx py pincy pc ps
    {-# INLINE rot #-}


instance BLAS1 (Complex Double) where
    copy n px incx py incy =
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
            zcopy pn px pincx py pincy
    {-# INLINE copy #-}

    swap n px incx py incy = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
            zswap pn px pincx py pincy
    {-# INLINE swap #-}
    
    dotc n px incx py incy =
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
        with 0 $ \pdotc -> do
            zdotc pdotc pn px pincx py pincy
            peek pdotc
    {-# INLINE dotc #-}

    dotu n px incx py incy =
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
        with 0 $ \pdotu -> do
            zdotu pdotu pn px pincx py pincy
            peek pdotu
    {-# INLINE dotu #-}

    nrm2 n px incx = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
            znrm2 pn px pincx
    {-# INLINE nrm2 #-}

    asum n px incx = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
            zasum pn px pincx
    {-# INLINE asum #-}

    iamax n px incx = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx -> do
            i <- izamax pn px pincx
            return $! fromEnum (i - 1)
    {-# INLINE iamax #-}

    axpy n alpha px incx py incy = 
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
            zaxpy pn palpha px pincx py pincy
    {-# INLINE axpy #-}
    
    scal n alpha px incx = 
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum incx $ \pincx ->
            zscal pn palpha px pincx
    {-# INLINE scal #-}

    rotg = zrotg
    {-# INLINE rotg #-}
    
    rot n px incx py incy c s = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        withEnum incy $ \pincy ->
        with c $ \pc ->
        with s $ \ps ->
            zdrot pn px pincx py pincy pc ps
    {-# INLINE rot #-}
