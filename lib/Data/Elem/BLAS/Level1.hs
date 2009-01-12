{-# LANGUAGE FlexibleInstances #-}
{-# OPTIONS_GHC -fno-excess-precision #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Elem.BLAS.Level1
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Vector operations.
--

module Data.Elem.BLAS.Level1
    where
     
import Control.Monad( liftM )
import Foreign ( Ptr, Storable, advancePtr, castPtr, peek, poke, with )
import Foreign.Storable.Complex ()
import Data.Complex ( Complex(..) )

import Data.Elem.BLAS.Base
import BLAS.Types
import BLAS.CTypes
import Data.Elem.BLAS.Double  
import Data.Elem.BLAS.Zomplex
        
-- | Types with vector-vector operations.
class (Elem a) => BLAS1 a where
    dot   :: ConjEnum -> ConjEnum -> Int -> Ptr a -> Int -> Ptr a -> Int -> IO a
    nrm2  :: Int -> Ptr a -> Int -> IO Double
    asum  :: Int -> Ptr a -> Int -> IO Double
    iamax :: Int -> Ptr a -> Int -> IO Int
    
    scal  :: Int -> a -> Ptr a -> Int -> IO () 

    axpy  :: Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    rotg  :: Ptr a -> Ptr a -> Ptr a -> Ptr a -> IO ()
    rot   :: Int -> Ptr a -> Int -> Ptr a -> Int -> Double -> Double -> IO ()

    -- | Replaces @y@ with @alpha (conj x) + y@
    acxpy :: Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    -- | Replaces @y@ with @x*y@.
    vmul :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    -- | Replaces @y@ with @conj(x)*y@.
    vcmul :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    -- | Replaces @y@ with @y/x@.
    vdiv :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()

    -- | Replaces @y@ with @y/conj(x)@.
    vcdiv :: Int -> Ptr a -> Int -> Ptr a -> Int -> IO ()



instance BLAS1 Double where
    dot _ _ = ddot
    nrm2    = dnrm2
    asum    = dasum
    iamax   = idamax
    axpy    = daxpy
    scal    = dscal
    rotg    = drotg
    rot     = drot
    acxpy   = daxpy
    vmul n  = dtbmv upper noTrans nonUnit n 0
    vcmul   = vmul
    vdiv n  = dtbsv upper noTrans nonUnit n 0
    vcdiv   = vdiv

instance BLAS1 (Complex Double) where
    dot conjX conjY n pX incX pY incY =
        case (conjX, conjY) of
            (NoConj, NoConj) -> dotc
            (Conj  , NoConj) -> dotu
            (Conj  , Conj  ) -> liftM conjugate dotc
            (NoConj, Conj  ) -> liftM conjugate dotu
      where
        dotu = with 0 $ \pDotu -> do
                   zdotu_sub n pX incX pY incY pDotu
                   peek pDotu

        dotc = with 0 $ \pDotc -> do
                   zdotc_sub n pX incX pY incY pDotc
                   peek pDotc
    {-# INLINE dot #-}

    nrm2  = znrm2
    asum  = zasum
    iamax = izamax
    
    axpy n alpha pX incX pY incY = 
        with alpha $ \pAlpha ->
            zaxpy n pAlpha pX incX pY incY
    
    scal n alpha pX incX =
        with alpha $ \pAlpha ->
            zscal n pAlpha pX incX
            
    rotg  = zrotg

    rot = zdrot

    acxpy n a pX incX pY incY =
        let pXR   = castPtr pX
            pYR   = castPtr pY
            pXI   = pXR `advancePtr` 1
            pYI   = pYR `advancePtr` 1
            incX' = 2 * incX
            incY' = 2 * incY
        in case a of
            (ra :+  0) -> do
                io <- daxpy n ( ra) pXR incX' pYR incY'
                io `seq` daxpy n (-ra) pXI incX' pYI incY'
            (0  :+ ia) -> do
                io <- daxpy n ( ia) pXR incX' pYI incY'
                io `seq` daxpy n ( ia) pXI incX' pYR incY'
            _ -> go n pX pY
        where
            go n' pX' pY' | n' <= 0 = return ()
                          | otherwise = do
                                x  <- peek pX'
                                y  <- peek pY'
                                io <- poke pY' (a * (conjugate x) + y)
                                
                                let n''  = n' - 1
                                    pX'' = pX' `advancePtr` incX
                                    pY'' = pY' `advancePtr` incY
                                    
                                io `seq` go n'' pX'' pY''
        
    vmul n  = ztbmv upper noTrans   nonUnit n 0
    vcmul n = ztbmv upper conjTrans nonUnit n 0

    vdiv n  = ztbsv upper noTrans   nonUnit n 0
    vcdiv n = ztbsv upper conjTrans nonUnit n 0
