-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Base (
    Matrix(..)
    ) where

import BLAS.Elem.Base ( Elem )

-- | A base class for matrices.
class Matrix a where
    -- | The number of rows in the matrix.
    numRows :: a (m,n) e -> Int
    
    -- | The number of columns in the matrix.
    numCols :: a (m,n) e -> Int
    
    -- | Creates a new matrix view that conjugates and transposes the 
    -- given matrix.
    herm :: Elem e => a (m,n) e -> a (n,m) e
