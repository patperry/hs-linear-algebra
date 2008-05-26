{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Elem
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Elem (
    Elem(..)
    ) where

import Data.Complex             ( Complex(..), conjugate, magnitude )
import Foreign                  ( Storable )
import Foreign.Storable.Complex ()

class (Storable e, Fractional e) => Elem e where
    -- | Take the complex conjugate of a value.  For real values
    -- this is equal to @id@.
    conj :: e -> e
    
    -- | Get the magnitude of a value.
    norm :: e -> Double
    
    -- | Convert a double to an element
    fromReal :: Double -> e

    -- | Coerce an element to a double
    toReal :: e -> Double
    

instance Elem Double where
    conj     = id
    norm     = abs
    fromReal = id
    toReal   = id
    
instance Elem (Complex Double) where
    conj             = conjugate
    norm             = magnitude
    fromReal x       = x :+ 0
    toReal  (x :+ _) = x
    