{-# LANGUAGE TypeFamilies, MultiParamTypeClasses, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Overloaded interface for matrices.  This module contains the common
-- functionality for the types defined in "Data.Matrix.Dense.Class" and
-- "Data.Matrix.Banded.Class", as well as provides a common base class
-- for the mutable and immutable matrix classes defined in submodules.
--

module Data.Matrix.Class (
    -- * Matrix shape
    MatrixShaped(..),
    numRows,
    numCols,
    isSquare,
    isFat,
    isTall,
    flipShape,

    -- * Associated types for matrices
    HasVectorView(..),
    HasMatrixStorage(..),

    -- * Matrix storage types
    module BLAS.Types,

    ) where

import BLAS.Types
import Data.Tensor.Class

-- | A class for matrices with an associated type for row, column, and
-- diagonal vector views.
class HasVectorView (a :: * -> * -> *) where
    -- | An associated type for a vector view into a matrix type @a@.  
    -- Typically, the view will share the same storage as the matrix,
    -- so that modifying an element in the view will affect the
    -- corresponding element in the matrix, and vice versa.
    type VectorView a :: * -> * -> *

-- | A class for matrix types that use a matrix internally for storage,
-- "Data.Matrix.Banded.Class" for example.
class HasMatrixStorage (a :: * -> * -> *) where
    -- | An associated type for the underlying matrix storage.
    type MatrixStorage a :: * -> * -> *

-- | A base class for objects shaped like matrices.
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

-- | Indicate whether or not a matrix has the same number of rows and columns.
isSquare :: (MatrixShaped a e) => a mn e -> Bool
isSquare a = numRows a == numCols a
{-# INLINE isSquare #-}

-- | Indicate whether or not the number of rows is less than or equal to 
-- the number of columns.
isFat :: (MatrixShaped a e) => a mn e -> Bool
isFat a = numRows a <= numCols a
{-# INLINE isFat #-}

-- | Indicate whether or not the number of rows is greater than or equal to 
-- the number of columns.
isTall :: (MatrixShaped a e) => a mn e -> Bool
isTall a = numRows a >= numCols a
{-# INLINE isTall #-}

-- | Replaces @(m,n)@ with @(n,m)@
flipShape :: (Int,Int) -> (Int,Int)
flipShape (m,n) = (n,m)
{-# INLINE flipShape #-}

