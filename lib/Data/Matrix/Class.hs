{-# LANGUAGE TypeFamilies, MultiParamTypeClasses, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Class (
    HasVectorView(..),
    MatrixShaped(..),
    numRows,
    numCols,
    isSquare,
    isFat,
    isTall,
    ) where

import Data.Tensor.Class

class HasVectorView (a :: * -> * -> *) where
    type VectorView a :: * -> * -> *

-- | A base class for matrices.
class (Shaped a (Int,Int) e) => MatrixShaped a e where
    -- | Creates a new matrix view that conjugates and transposes the 
    -- given matrix.
    herm :: a (m,n) e -> a (n,m) e

-- | Get the number of rows in the matrix.
numRows :: (MatrixShaped a e) => a mn e -> Int
numRows = fst . shape
{-# INLINE numRows #-}

-- | Get the number of rows in the matrix.
numCols :: (MatrixShaped a e) => a mn e -> Int
numCols = snd . shape
{-# INLINE numCols #-}

isSquare :: (MatrixShaped a e) => a mn e -> Bool
isSquare a = numRows a == numCols a

isFat :: (MatrixShaped a e) => a mn e -> Bool
isFat a = numRows a <= numCols a

isTall :: (MatrixShaped a e) => a mn e -> Bool
isTall a = numRows a >= numCols a

