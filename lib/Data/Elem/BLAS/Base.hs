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
    module Data.Elem.Conj,
    ) where

import Data.Elem.Conj
import Data.Complex             ( Complex(..), magnitude )
import Foreign                  ( Storable )
import Foreign.Storable.Complex ()

-- | The base class for elements.
class (Storable e, Fractional e, Conj e) => Elem e where
    -- | Get the magnitude of a value.
    norm :: e -> Double
    
    -- | Get the l1 norm of a value.
    norm1 :: e -> Double
    
    -- | Convert a double to an element
    fromReal :: Double -> e

    -- | Coerce an element to a double
    toReal :: e -> Double
    
instance Elem Double where
    norm     = abs
    norm1    = abs
    fromReal = id
    toReal   = id
    
instance Elem (Complex Double) where
    norm             = magnitude
    norm1 (x :+ y)   = abs x + abs y
    fromReal x       = x :+ 0
    toReal  (x :+ _) = x
    