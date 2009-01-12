{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Elem.BLAS.Base
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Elem.BLAS.Base (
    Elem(..),
    ) where

import Control.Monad
import Data.AEq
import Data.Complex( Complex(..), magnitude )
import qualified Data.Complex as Complex
import Foreign
import Foreign.Storable.Complex ()


import BLAS.Types( ConjEnum(..) )
import Data.Elem.BLAS.Double( dcopy, dswap, dscal )
import Data.Elem.BLAS.Zomplex( zcopy, zswap )


-- | The base class for elements.
class (AEq e, Storable e, Fractional e) => Elem e where
    -- | Get the complex conjugate of a value.
    conjugate :: e -> e
    
    -- | Get the magnitude of a value.
    norm :: e -> Double
    
    -- | Get the l1 norm of a value.
    norm1 :: e -> Double
    
    -- | Convert a double to an element.
    fromReal :: Double -> e

    -- | Try to coerce a value to a double.  This will fail unless the
    -- complex part is zero (according to a comparison by @(~==)@).
    maybeToReal :: e -> Maybe Double

    copy  :: ConjEnum -> ConjEnum -> Int -> Ptr e -> Int -> Ptr e -> Int -> IO ()
    swap  :: ConjEnum -> ConjEnum -> Int -> Ptr e -> Int -> Ptr e -> Int -> IO ()
    vconj :: Int -> Ptr e -> Int -> IO ()


    
instance Elem Double where
    conjugate   = id
    {-# INLINE conjugate #-}
    norm        = abs
    {-# INLINE norm #-}
    norm1       = abs
    {-# INLINE norm1 #-}
    fromReal    = id
    {-# INLINE fromReal #-}
    maybeToReal = Just
    {-# INLINE maybeToReal #-}
    copy _ _ = dcopy
    {-# INLINE copy #-}
    swap _ _ = dswap
    {-# INLINE swap #-}
    vconj _ _ _ = return ()
    {-# INLINE vconj #-}
    
instance Elem (Complex Double) where
    conjugate      = Complex.conjugate
    {-# INLINE conjugate #-}
    norm           = magnitude
    {-# INLINE norm #-}    
    norm1 (x :+ y) = abs x + abs y
    {-# INLINE norm1 #-}    
    fromReal x     = x :+ 0
    {-# INLINE fromReal #-}    
    maybeToReal (x :+ y) | y ~== 0   = Just x
                         | otherwise = Nothing
    {-# INLINE maybeToReal #-}
    
    copy conjX conjY n pX incX pY incY = do
        io <- zcopy n pX incX pY incY
        io `seq` when (conjX /= conjY) $ 
                     vconj n pY incY
    {-# INLINE copy #-}
    
    swap conjX conjY n pX incX pY incY = do
        io <- zswap n pX incX pY incY
        io `seq` when (conjX /= conjY) $ do
                     io1 <- vconj n pX incX
                     io1 `seq` vconj n pY incY
    {-# INLINE swap #-}
    
    vconj n pX incX =
        let pXI   = (castPtr pX) `advancePtr` 1
            alpha = -1
            incXI = 2 * incX
        in dscal n alpha pXI incXI
    {-# INLINE vconj #-}
