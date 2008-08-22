-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.RowCol.View
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.RowCol.View (
    rowView,
    colView,
    rowViews,
    colViews,
    
    RowColView(..),
    ) where

import BLAS.Internal( checkedRow, checkedCol )
import BLAS.Tensor( shape )
import BLAS.Matrix.Base
import Data.Vector.Dense.Class( BaseVector )

class (BaseMatrix a e, BaseVector x e) => RowColView a x e | a -> x where
    -- | Same as 'rowView' but index is not range-checked.
    unsafeRowView :: a (k,l) e -> Int -> x l e
    
    -- | Same as 'colView' but index is not range-checked.    
    unsafeColView :: a (k,l) e -> Int -> x k e


-- | Get a vector view of the given row in a matrix.
rowView :: (RowColView a x e) => a (m,n) e -> Int -> x n e
rowView a = checkedRow (shape a) (unsafeRowView a)

-- | Get a vector view of the given column in a matrix.
colView :: (RowColView a x e) => a (m,n) e -> Int -> x m e
colView a = checkedCol (shape a) (unsafeColView a)

-- | Get a list of vector views of the rows of the matrix.
rowViews :: (RowColView a x e) => a (m,n) e -> [x n e]
rowViews a = [ unsafeRowView a i | i <- [0..numRows a - 1] ]

-- | Get a list of vector views of the columns of the matrix.
colViews :: (RowColView a x e) => a (m,n) e -> [x m e]
colViews a = [ unsafeColView a j | j <- [0..numCols a - 1] ]
