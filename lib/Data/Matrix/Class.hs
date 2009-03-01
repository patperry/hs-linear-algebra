{-# LANGUAGE TypeFamilies, MultiParamTypeClasses, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Common functionality for the types defined in
-- "Data.Matrix.Dense.Class" and "Data.Matrix.Banded.Class", and 
-- a base class for the mutable and immutable matrix
-- classes defined in the submodules of this one.
--

module Data.Matrix.Class (
    -- * Matrix shape
    MatrixShaped,
    numRows,
    numCols,
    isSquare,
    isFat,
    isTall,
    flipShape,
    
    -- * Hermitian transpose
    HasHerm(..),

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
class HasVectorView (a :: * -> *) where
    -- | An associated type for a vector view into a matrix type @a@.  
    -- Typically, the view will share the same storage as the matrix,
    -- so that modifying an element in the view will affect the
    -- corresponding element in the matrix, and vice versa.
    type VectorView a :: * -> *

-- | A class for matrix types that use a matrix internally for storage,
-- "Data.Matrix.Banded.Class" for example.
class HasMatrixStorage (a :: * -> *) where
    -- | An associated type for the underlying matrix storage.
    type MatrixStorage a :: * -> *

-- | A base class for objects shaped like matrices.
class (Shaped a (Int,Int)) => MatrixShaped a
    
-- | A class for objects with a hermitian transpose
class (MatrixShaped a) => HasHerm a where
    -- | Creates a new matrix view that conjugates and transposes the 
    -- given matrix.
    herm :: a e -> a e

-- | Get the number of rows in the matrix.
numRows :: (MatrixShaped a) => a e -> Int
numRows = fst . shape
{-# INLINE numRows #-}

-- | Get the number of rows in the matrix.
numCols :: (MatrixShaped a) => a e -> Int
numCols = snd . shape
{-# INLINE numCols #-}

-- | Indicate whether or not a matrix has the same number of rows and columns.
isSquare :: (MatrixShaped a) => a e -> Bool
isSquare a = numRows a == numCols a
{-# INLINE isSquare #-}

-- | Indicate whether or not the number of rows is less than or equal to 
-- the number of columns.
isFat :: (MatrixShaped a) => a e -> Bool
isFat a = numRows a <= numCols a
{-# INLINE isFat #-}

-- | Indicate whether or not the number of rows is greater than or equal to 
-- the number of columns.
isTall :: (MatrixShaped a) => a e -> Bool
isTall a = numRows a >= numCols a
{-# INLINE isTall #-}

-- | Replaces @(m,n)@ with @(n,m)@
flipShape :: (Int,Int) -> (Int,Int)
flipShape (m,n) = (n,m)
{-# INLINE flipShape #-}

