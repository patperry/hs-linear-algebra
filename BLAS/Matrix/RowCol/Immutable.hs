-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.RowCol.Immutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.RowCol.Immutable (
    row,
    col,
    rows,
    cols,
    
    IRowCol(..),
    ) where

import BLAS.Internal( checkedRow, checkedCol )

import BLAS.Matrix.Base
import Data.Vector.Dense

class (BaseMatrix a e) => IRowCol a e where
    -- | Same as 'row' but index is not range-checked.
    unsafeRow :: a (m,n) e -> Int -> Vector n e
    
    -- | Same as 'col' but index is not range-checked.    
    unsafeCol :: a (m,n) e -> Int -> Vector m e


-- | Get the given row in a matrix.
row :: (IRowCol a e) => a (m,n) e -> Int -> Vector n e
row a = checkedRow (shape a) (unsafeRow a)

-- | Get the given column in a matrix.
col :: (IRowCol a e) => a (m,n) e -> Int -> Vector m e
col a = checkedCol (shape a) (unsafeCol a)

-- | Get a list the row vectors in the matrix.
rows :: (IRowCol a e) => a (m,n) e -> [Vector n e]
rows a = [ unsafeRow a i | i <- [0..numRows a - 1] ]

-- | Get a list the column vectors in the matrix.
cols :: (IRowCol a e) => a (m,n) e -> [Vector m e]
cols a = [ unsafeCol a j | j <- [0..numCols a - 1] ]
