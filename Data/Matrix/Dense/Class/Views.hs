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
    submatrix,
    unsafeSubmatrix,

    -- * Row and Column views
    module BLAS.Matrix.RowCol.View,
    module BLAS.Matrix.RowCol.Read,
    module BLAS.Matrix.Diag.View,
    module BLAS.Matrix.Diag.Read,
    
    ) where

import BLAS.Internal( checkedSubmatrix )
import BLAS.Tensor( shape )
import BLAS.Matrix( herm )

import BLAS.Matrix.RowCol.View
import BLAS.Matrix.RowCol.Read
import BLAS.Matrix.Diag.View
import BLAS.Matrix.Diag.Read

import Data.Matrix.Dense.Class.Internal


-- | @submatrix a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrix :: (BaseMatrix a x e) => a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
submatrix a = checkedSubmatrix (shape a) (unsafeSubmatrix a)
{-# INLINE submatrix #-}

-- | Same as 'submatrix' but indices are not range-checked.
unsafeSubmatrix :: (BaseMatrix a x e) => 
    a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
unsafeSubmatrix a (i,j) (m,n)
    | isHerm a  = 
        coerceMatrix $ herm $ 
            unsafeSubmatrix (herm $ coerceMatrix a) (j,i) (n,m)
    | otherwise =
        let f = fptrOfMatrix a
            o = indexOfMatrix a (i,j)
            l = lda a
        in matrixViewArray f o (m,n) l False
