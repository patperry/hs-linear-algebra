-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Diag.Immutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Diag.Immutable (
    diag,
    IDiag(..),
    ) where

import BLAS.Internal( checkedDiag )

import BLAS.Matrix.Base
import Data.Vector.Dense

class (BaseMatrix a e) => IDiag a e where
    -- | Same as 'diag' but index is not range-checked.
    unsafeDiag :: a mn e -> Int -> Vector k e
    
-- | Get the given row in a matrix.
diag :: (IDiag a e) => a mn e -> Int -> Vector k e
diag a = checkedDiag (shape a) (unsafeDiag a)
