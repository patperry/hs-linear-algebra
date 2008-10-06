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
    unsafeSubmatrixView,

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

import BLAS.Elem( BLAS1 )
import BLAS.Internal( checkedSubmatrix, checkedRow, checkedCol, checkedDiag )
import BLAS.Tensor( shape )
import BLAS.Matrix.Base( herm )

import Data.Matrix.Dense.Class.Internal
import Data.Vector.Dense.Class.Internal( WriteVector, newCopyVector )


-- | @submatrixView a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrixView :: (BaseMatrix a x e) => a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
submatrixView a = checkedSubmatrix (shape a) (unsafeSubmatrixView a)
{-# INLINE submatrixView #-}

-- | Same as 'submatrixView' but indices are not range-checked.
unsafeSubmatrixView :: (BaseMatrix a x e) => 
    a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
unsafeSubmatrixView a (i,j) (m,n)
    | isHermMatrix a  = 
        coerceMatrix $ herm $ 
            unsafeSubmatrixView (herm $ coerceMatrix a) (j,i) (n,m)
    | otherwise =
        let f = fptrOfMatrix a
            o = indexOfMatrix a (i,j)
            l = ldaOfMatrix a
        in matrixViewArray f o (m,n) l False


-- | Get a vector view of the given diagonal in a matrix.
diagView :: (BaseMatrix a x e) => a mn e -> Int -> x k e
diagView a = checkedDiag (shape a) (unsafeDiagView a)

-- | Get a vector view of the given row in a matrix.
rowView :: (BaseMatrix a x e) => a (m,n) e -> Int -> x n e
rowView a = checkedRow (shape a) (unsafeRowView a)

-- | Get a vector view of the given column in a matrix.
colView :: (BaseMatrix a x e) => a (m,n) e -> Int -> x m e
colView a = checkedCol (shape a) (unsafeColView a)

-- | Get the given diagonal in a matrix.  Negative indices correspond
-- to sub-diagonals.
getDiag :: (ReadMatrix a x e m, WriteVector y e m, BLAS1 e) => 
    a mn e -> Int -> m (y k e)
getDiag a = checkedDiag (shape a) (unsafeGetDiag a)

-- | Same as 'getDiag' but not range-checked.
unsafeGetDiag :: (ReadMatrix a x e m, WriteVector y e m, BLAS1 e) => 
    a mn e -> Int -> m (y k e)
unsafeGetDiag a i = newCopyVector (unsafeDiagView a i)
