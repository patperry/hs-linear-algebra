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
import Data.Vector.Dense.Class.Internal.Base( BaseVector, conjVector )

class Conj e where
    -- | Take the complex conjugate of a value.  For real values
    -- this is equal to @id@.
    conj :: e -> e

instance Conj (Complex Double) where
    conj = conjugate
    {-# INLINE conj #-}
    
instance Conj Double where
    conj = id
    {-# INLINE conj #-}

instance (BaseVector x e) => Conj (x n e) where
    conj = conjVector
    {-# INLINE conj #-}
