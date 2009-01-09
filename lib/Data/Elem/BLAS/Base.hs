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

import Data.AEq
import Data.Complex             ( Complex(..), magnitude )
import qualified Data.Complex as Complex
import Foreign                  ( Storable )
import Foreign.Storable.Complex ()

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
    
    -- | Inicates whether or not the value should be used in tests.  For
    -- 'Double's, @isTestableElem e@ is defined as 
    -- @not (isNaN e || isInfinite e || isDenormalized e)@.
    isTestableElem :: e -> Bool
    
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
    isTestableElem e = not (isNaN e || isInfinite e || isDenormalized e)
    {-# INLINE isTestableElem #-}
    
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
    isTestableElem (x :+ y) = isTestableElem x && isTestableElem y
    {-# INLINE isTestableElem #-}
