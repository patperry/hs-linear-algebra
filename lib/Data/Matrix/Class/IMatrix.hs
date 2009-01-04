{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class.IMatrix
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Class.IMatrix (
    -- * Rows and columns
    row,
    col,
    rows,
    cols,

    -- * Multiplication
    (<*>),
    (<**>),
    sapply,
    sapplyMat,
    
    IMatrix(..),
    unsafeApply,
    unsafeApplyMat,
    ) where

import Data.Elem.BLAS( BLAS1, BLAS3 )
import BLAS.Internal ( checkedRow, checkedCol, checkMatVecMult, 
    checkMatMatMult )

import Data.Matrix.Class.MMatrix( unsafeGetSApply, unsafeGetSApplyMat )

import Data.Vector.Dense
import Data.Vector.Dense.ST( runSTVector )
import Data.Matrix.Herm
import Data.Matrix.Tri.Internal
import Data.Matrix.Dense.Internal
import Data.Matrix.Dense.Class( unsafeRowView, unsafeColView )
import Data.Matrix.Dense.ST( runSTMatrix )

infixr 7 <*>, <**>

class (MatrixShaped a e, BLAS1 e) => IMatrix a e where
    unsafeSApply :: e -> a (m,n) e -> Vector n e -> Vector m e
    unsafeSApplyMat :: e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e

    -- | Same as 'row' but index is not range-checked.
    unsafeRow :: a (m,n) e -> Int -> Vector n e
    unsafeRow a i = let
        e = basisVector (numRows a) i
        in conj $ unsafeApply (herm a) e
    
    -- | Same as 'col' but index is not range-checked.    
    unsafeCol :: a (m,n) e -> Int -> Vector m e
    unsafeCol a j = let
        e = basisVector (numCols a) j
        in unsafeApply a e


-- | Get the given row in a matrix.
row :: (IMatrix a e) => a (m,n) e -> Int -> Vector n e
row a = checkedRow (shape a) (unsafeRow a)

-- | Get the given column in a matrix.
col :: (IMatrix a e) => a (m,n) e -> Int -> Vector m e
col a = checkedCol (shape a) (unsafeCol a)

-- | Get a list the row vectors in the matrix.
rows :: (IMatrix a e) => a (m,n) e -> [Vector n e]
rows a = [ unsafeRow a i | i <- [0..numRows a - 1] ]

-- | Get a list the column vectors in the matrix.
cols :: (IMatrix a e) => a (m,n) e -> [Vector m e]
cols a = [ unsafeCol a j | j <- [0..numCols a - 1] ]


-- | Apply to a vector
(<*>) :: (IMatrix a e) => a (m,n) e -> Vector n e -> Vector m e
(<*>) a x = checkMatVecMult (shape a) (dim x) $ unsafeApply a x
    
-- | Apply to a matrix
(<**>) :: (IMatrix a e) => a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
(<**>) a b = checkMatMatMult (shape a) (shape b) $ unsafeApplyMat a b

sapply :: (IMatrix a e) => e -> a (m,n) e -> Vector n e -> Vector m e
sapply k a x = checkMatVecMult (shape a) (dim x) $ unsafeSApply k a x
    
sapplyMat :: (IMatrix a e) => e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e    
sapplyMat k a b = checkMatMatMult (shape a) (shape b) $ unsafeSApplyMat k a b

unsafeApply :: (IMatrix a e) => a (m,n) e -> Vector n e -> Vector m e
unsafeApply = unsafeSApply 1

unsafeApplyMat :: (IMatrix a e) => a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
unsafeApplyMat = unsafeSApplyMat 1

instance (BLAS3 e) => IMatrix Matrix e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b
    unsafeRow                 = unsafeRowView
    unsafeCol                 = unsafeColView

instance (BLAS3 e) => IMatrix (Herm Matrix) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    

instance (BLAS3 e) => IMatrix (Tri Matrix) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    

{-# RULES
"scale.apply/sapply"       forall k a x. (<*>) (k *> a) x = sapply k a x
"scale.applyMat/sapplyMat" forall k a b. (<**>) (k *> a) b = sapplyMat k a b
  #-}
