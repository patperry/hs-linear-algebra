-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Vector
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Vector (
    Vector(..)
    ) where

import BLAS.Elem.Base ( Elem )

-- | A class for vectors.
class Vector x where
    -- | Get the dimension of the vector.
    dim :: x n e -> Int
    
    -- | Get the complex conjugate of a vector.
    conj :: (Elem e) => x n e -> x n e
     