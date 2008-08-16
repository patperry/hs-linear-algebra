{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Conj
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Conj (
    Conj(..)
    ) where

import Data.Complex


class Conj e where
    -- | Take the complex conjugate of a value.  For real values
    -- this is equal to @id@.
    conj :: e -> e

instance Conj Double where
    conj     = id
    
instance Conj (Complex Double) where
    conj             = conjugate
