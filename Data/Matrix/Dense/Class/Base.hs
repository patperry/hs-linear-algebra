-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Base (
    BaseMatrix(..),
    
    submatrix,
    coerceMatrix,
    ) where

import BLAS.Internal( checkedSubmatrix )
import BLAS.Tensor( shape )

import BLAS.Matrix hiding ( BaseMatrix )
import qualified BLAS.Matrix as BLAS

import Foreign ( ForeignPtr, Ptr )
import Unsafe.Coerce

class (BLAS.BaseMatrix a e) => BaseMatrix a e where
    lda :: a mn e -> Int
    isHerm :: a mn e -> Bool
    unsafeSubmatrix :: a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
    withMatrixPtr :: a mn e -> (Ptr e -> IO b) -> IO b

    matrixViewArray :: ForeignPtr e -> Int -> (Int,Int) -> Int -> Bool -> a mn e
    arrayFromMatrix :: a mn e -> (ForeignPtr e, Int, (Int,Int), Int, Bool)
    
-- | @submatrix a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrix :: (BaseMatrix a e) => a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
submatrix a = checkedSubmatrix (shape a) (unsafeSubmatrix a)
{-# INLINE submatrix #-}

-- | Cast the shape type of the matrix.
coerceMatrix :: (BaseMatrix a e) => a mn e -> a mn' e
coerceMatrix = unsafeCoerce
{-# INLINE coerceMatrix #-}
