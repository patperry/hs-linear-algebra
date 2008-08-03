-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Base (
    Matrix(..),
    isSquare,
    isFat,
    isTall,
    ) where

import BLAS.Elem.Base ( Elem )
import BLAS.Tensor.Base

-- | A base class for matrices.
class Matrix a where
    -- | The number of rows in the matrix.
    numRows :: a (m,n) e -> Int
    
    -- | The number of columns in the matrix.
    numCols :: a (m,n) e -> Int
    
    -- | Creates a new matrix view that conjugates and transposes the 
    -- given matrix.
    herm :: Elem e => a (m,n) e -> a (n,m) e


isSquare :: (Matrix a) => a (m,n) e -> Bool
isSquare a = numRows a == numCols a

isFat :: (Matrix a) => a (m,n) e -> Bool
isFat a = numRows a <= numCols a

isTall :: (Matrix a) => a (m,n) e -> Bool
isTall a = numRows a >= numCols a

instance (Matrix a) => Tensor (a (m,n)) (Int,Int) e where
    shape a = (numRows a, numCols a)
    
    bounds a = ((0,0), (m-1,n-1))
      where (m,n) = shape a
      