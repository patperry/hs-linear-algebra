{-# LANGUAGE FlexibleInstances #-}
{-# OPTIONS_GHC -fno-excess-precision #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.C.Level1
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.C.Level1
    where
     
import Prelude hiding ( div )

import Foreign ( Ptr, Storable, advancePtr, castPtr, peek, poke, with )
import Foreign.Storable.Complex ()
import Data.Complex

import BLAS.Elem.Base
import BLAS.C.Types
import BLAS.C.Double  
import BLAS.C.Zomplex
        
class (Elem a) => BLAS1 a where
    dotu  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO a
    dotc  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO a
    nrm2  :: Int -> Ptr a -> Int -> IO Double
    asum  :: Int -> Ptr a -> Int -> IO Double
    iamax :: Int -> Ptr a -> Int -> IO Int
    
    scal  :: Int -> a -> Ptr a -> Int -> IO () 

    swap  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    copy  :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    axpy  :: Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    rotg  :: Ptr a -> Ptr a -> Ptr a -> Ptr a -> IO ()
    rot   :: Int -> Ptr a -> Int -> Ptr a -> Int -> Double -> Double -> IO ()

    -- | conjugate all elements of a vector
    conj  :: Int -> Ptr a -> Int -> IO ()

    -- | Replaces @y@ with @alpha (conj x) + y@
    acxpy :: Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    -- | Replaces @y@ with @x*y@.
    mul :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    -- | Replaces @y@ with @conj(x)*y@.
    cmul :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    -- | Replaces @y@ with @y/x@.
    div :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    -- | Replaces @y@ with @y/conj(x)@.
    cdiv :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()



instance BLAS1 Double where
    dotu  = ddot
    dotc  = ddot
    nrm2  = dnrm2
    asum  = dasum
    iamax = idamax
    swap  = dswap
    copy  = dcopy
    axpy  = daxpy
    scal  = dscal
    rotg  = drotg
    rot   = drot
    conj _ _ _ = return ()
    acxpy = daxpy
    mul n = dtbmv colMajor upper noTrans nonUnit n 0
    cmul  = mul
    div n = dtbsv colMajor upper noTrans nonUnit n 0
    cdiv = div

instance BLAS1 (Complex Double) where
    dotu n pX incX pY incY =
        with 0 $ \pDotu -> do
            zdotu_sub n pX incX pY incY pDotu
            peek pDotu

    dotc n pX incX pY incY =
        with 0 $ \pDotc -> do
            zdotc_sub n pX incX pY incY pDotc
            peek pDotc

    nrm2  = znrm2
    asum  = zasum
    iamax = izamax
    swap  = zswap
    copy  = zcopy
    
    axpy n alpha pX incX pY incY = 
        with alpha $ \pAlpha ->
            zaxpy n pAlpha pX incX pY incY
    
    scal n alpha pX incX =
        with alpha $ \pAlpha ->
            zscal n pAlpha pX incX
            
    rotg  = zrotg

    rot = zdrot

    conj n pX incX =
        let pXI   = (castPtr pX) `advancePtr` 1
            alpha = -1
            incXI = 2 * incX
        in dscal n alpha pXI incXI
    
    acxpy n a pX incX pY incY =
        let pXR   = castPtr pX
            pYR   = castPtr pY
            pXI   = pXR `advancePtr` 1
            pYI   = pYR `advancePtr` 1
            incX' = 2 * incX
            incY' = 2 * incY
        in case a of
            (ra :+  0) -> do
                daxpy n ( ra) pXR incX' pYR incY'
                daxpy n (-ra) pXI incX' pYI incY'
            (0  :+ ia) -> do
                daxpy n ( ia) pXR incX' pYI incY'
                daxpy n ( ia) pXI incX' pYR incY'
            _ -> go n pX pY
        where
            go n' pX' pY'
                | n' `seq` pX' `seq` pY' `seq` False = undefined
                | n' <= 0 =
                    return ()
                | otherwise = do
                    x <- peek pX'
                    y <- peek pY'
                    poke pY' (a * (conjugate x) + y)
                    
                    let n''  = n' - 1
                        pX'' = pX' `advancePtr` incX
                        pY'' = pY' `advancePtr` incY
                        
                    go n'' pX'' pY''
        
    mul n  = ztbmv colMajor upper noTrans   nonUnit n 0
    cmul n = ztbmv colMajor upper conjTrans nonUnit n 0

    div n  = ztbsv colMajor upper noTrans   nonUnit n 0
    cdiv n = ztbsv colMajor upper conjTrans nonUnit n 0
