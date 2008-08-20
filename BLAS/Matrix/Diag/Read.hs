-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Diag.Read
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Diag.Read (
    getDiag,
    DiagRead(..),
    ) where

import BLAS.Internal( checkedDiag )
import BLAS.Tensor( shape )

import BLAS.Matrix.Base
import Data.Vector.Dense.Class

class (BaseMatrix a e, WriteVector x e m) => DiagRead a x e m where
    -- | Same as 'getDiag' but index is not range-checked.
    unsafeGetDiag :: a (k,l) e -> Int -> m (x n e)

-- | Get the given diagonal in a matrix.  Negative indices correspond
-- to sub-diagonals.
getDiag :: (DiagRead a x e m) => a (k,l) e -> Int -> m (x n e)
getDiag a = checkedDiag (shape a) (unsafeGetDiag a)
