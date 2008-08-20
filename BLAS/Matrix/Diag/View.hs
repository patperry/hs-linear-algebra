-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Diag.View
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Diag.View (
    diagView,
    DiagView(..),
    ) where

import BLAS.Internal( checkedDiag )
import BLAS.Tensor( shape )
import BLAS.Matrix.Base
import Data.Vector.Dense.Class

class (BaseMatrix a e, BaseVector x e) => DiagView a x e | a -> x where
    -- | Same as 'diagView' but index is not range-checked.    
    unsafeDiagView :: a (k,l) e -> Int -> x n e

-- | Get a vector view of the given diagonal in a matrix.
diagView :: (DiagView a x e) => a (k,l) e -> Int -> x n e
diagView a = checkedDiag (shape a) (unsafeDiagView a)




