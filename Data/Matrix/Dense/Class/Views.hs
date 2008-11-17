-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Views
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Views (
    -- * Matrix views
    submatrixView,
    splitRowsAt,
    splitColsAt,
    unsafeSubmatrixView,
    unsafeSplitRowsAt,
    unsafeSplitColsAt,

    -- * Row and Column views
    rowViews,
    colViews,
    rowView,
    colView,
    diagView,
    unsafeRowView,
    unsafeColView,
    unsafeDiagView,
    
    -- * Getting rows and columns
    getDiag,
    unsafeGetDiag,
    
    ) where

import BLAS.Internal( checkedSubmatrix, checkedRow, checkedCol, checkedDiag )
import BLAS.Tensor( shape )
import BLAS.Matrix.Shaped( herm )

import Data.Matrix.Dense.Class.Internal
import Data.Vector.Dense.Class.Internal( WriteVector, newCopyVector )

import Foreign


-- | @submatrixView a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrixView :: (BaseMatrix a e) => a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
submatrixView a = checkedSubmatrix (shape a) (unsafeSubmatrixView a)
{-# INLINE submatrixView #-}

-- | Same as 'submatrixView' but indices are not range-checked.
unsafeSubmatrixView :: (BaseMatrix a e) => 
    a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
unsafeSubmatrixView a (i,j) (m,n)
    | isHermMatrix a  = 
        coerceMatrix $ herm $ 
            unsafeSubmatrixView (herm $ coerceMatrix a) (j,i) (n,m)
    | otherwise =
        let (fp,p,_,_,ld,_) = arrayFromMatrix a
            o  = indexOfMatrix a (i,j)
            p' = p `advancePtr` o
        in matrixViewArray fp p' m n ld False

splitRowsAt :: (BaseMatrix a e) =>
    Int -> a (m,n) e -> (a (m1,n) e, a (m2,n) e)
splitRowsAt m1 a = ( submatrixView a (0,0)  (m1,n)
                   , submatrixView a (m1,0) (m2,n)
                   )
  where 
    (m,n) = shape a
    m2    = m - m1

unsafeSplitRowsAt :: (BaseMatrix a e) =>
    Int -> a (m,n) e -> (a (m1,n) e, a (m2,n) e)
unsafeSplitRowsAt m1 a = ( unsafeSubmatrixView a (0,0)  (m1,n)
                         , unsafeSubmatrixView a (m1,0) (m2,n)
                         )
  where 
    (m,n) = shape a
    m2    = m - m1

splitColsAt :: (BaseMatrix a e) =>
    Int -> a (m,n) e -> (a (m,n1) e, a (m,n2) e)
splitColsAt n1 a = ( submatrixView a (0,0)  (m,n1)
                   , submatrixView a (0,n1) (m,n2)
                   )
  where
    (m,n) = shape a
    n2    = n - n1

unsafeSplitColsAt :: (BaseMatrix a e) =>
    Int -> a (m,n) e -> (a (m,n1) e, a (m,n2) e)
unsafeSplitColsAt n1 a = ( unsafeSubmatrixView a (0,0)  (m,n1)
                         , unsafeSubmatrixView a (0,n1) (m,n2)
                         )
  where
    (m,n) = shape a
    n2    = n - n1


-- | Get a vector view of the given diagonal in a matrix.
diagView :: (BaseMatrix a e) => a mn e -> Int -> VectorView a k e
diagView a = checkedDiag (shape a) (unsafeDiagView a)

-- | Get a vector view of the given row in a matrix.
rowView :: (BaseMatrix a e) => a (m,n) e -> Int -> VectorView a n e
rowView a = checkedRow (shape a) (unsafeRowView a)

-- | Get a vector view of the given column in a matrix.
colView :: (BaseMatrix a e) => a (m,n) e -> Int -> VectorView a m e
colView a = checkedCol (shape a) (unsafeColView a)

-- | Get the given diagonal in a matrix.  Negative indices correspond
-- to sub-diagonals.
getDiag :: (ReadMatrix a e m, WriteVector y e m) => 
    a mn e -> Int -> m (y k e)
getDiag a = checkedDiag (shape a) (unsafeGetDiag a)

-- | Same as 'getDiag' but not range-checked.
unsafeGetDiag :: (ReadMatrix a e m, WriteVector y e m) => 
    a mn e -> Int -> m (y k e)
unsafeGetDiag a i = newCopyVector (unsafeDiagView a i)
