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
import Data.Vector.Dense.Class( WriteVector )

class (BaseMatrix a e) => DiagRead a e m where
    -- | Same as 'getDiag' but index is not range-checked.
    unsafeGetDiag :: (WriteVector x e m) => a mn e -> Int -> m (x k e)

-- | Get the given diagonal in a matrix.  Negative indices correspond
-- to sub-diagonals.
getDiag :: (DiagRead a e m, WriteVector x e m) => a mn e -> Int -> m (x k e)
getDiag a = checkedDiag (shape a) (unsafeGetDiag a)
