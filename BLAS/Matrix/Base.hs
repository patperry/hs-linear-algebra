{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Base (
    BaseMatrix(..),
    numRows,
    numCols,
    isSquare,
    isFat,
    isTall,
    ) where

import BLAS.Tensor

-- | A base class for matrices.
class (BaseTensor a (Int,Int) e) => BaseMatrix a e where
    -- | Creates a new matrix view that conjugates and transposes the 
    -- given matrix.
    herm :: a (m,n) e -> a (n,m) e

-- | Get the number of rows in the matrix.
numRows :: (BaseMatrix a e) => a mn e -> Int
numRows = fst . shape
{-# INLINE numRows #-}

-- | Get the number of rows in the matrix.
numCols :: (BaseMatrix a e) => a mn e -> Int
numCols = snd . shape
{-# INLINE numCols #-}

isSquare :: (BaseMatrix a e) => a mn e -> Bool
isSquare a = numRows a == numCols a

isFat :: (BaseMatrix a e) => a mn e -> Bool
isFat a = numRows a <= numCols a

isTall :: (BaseMatrix a e) => a mn e -> Bool
isTall a = numRows a >= numCols a

