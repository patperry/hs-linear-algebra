{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class.IMatrixBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface for immutable matrices.  The matrices provide
-- access to rows and columns, and can operate via multiplication on 
-- immutable dense vectors and matrices.
--

module Data.Matrix.Class.IMatrixBase (
    -- * The IMatrix type class
    IMatrix(..),

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
    
    -- * Unsafe operations
    unsafeApply,
    unsafeApplyMat,
    ) where

import Data.Elem.BLAS
import BLAS.Internal ( checkedRow, checkedCol, checkMatVecMult, 
    checkMatMatMult )

import Data.Tensor.Class

import Data.Vector.Dense
import Data.Matrix.Dense.Base

infixr 7 <*>, <**>

-- | Get the given row in a matrix.
row :: (IMatrix a, Elem e) => a (m,n) e -> Int -> Vector n e
row a = checkedRow (shape a) (unsafeRow a)
{-# INLINE row #-}

-- | Get the given column in a matrix.
col :: (IMatrix a, Elem e) => a (m,n) e -> Int -> Vector m e
col a = checkedCol (shape a) (unsafeCol a)
{-# INLINE col #-}

-- | Matrix multiplication by a vector.
(<*>) :: (IMatrix a, BLAS3 e) => a (m,n) e -> Vector n e -> Vector m e
(<*>) a x = checkMatVecMult (shape a) (dim x) $ unsafeApply a x
{-# INLINE (<*>) #-}

-- | Matrix multiplication by a matrix.
(<**>) :: (IMatrix a, BLAS3 e) => a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
(<**>) a b = checkMatMatMult (shape a) (shape b) $ unsafeApplyMat a b
{-# INLINE (<**>) #-}

-- | Scale and multiply by a vector.  
-- @sapply k a x@ is equal to @a \<*> (k *> x)@, and often it is faster.
sapply :: (IMatrix a, BLAS3 e) => e -> a (m,n) e -> Vector n e -> Vector m e
sapply k a x = checkMatVecMult (shape a) (dim x) $ unsafeSApply k a x
{-# INLINE sapply #-}
    
-- | Scale and multiply by a matrix.
-- @sapplyMat k a b@ is equal to @a \<**> (k *> b)@, and often it is faster.
sapplyMat :: (IMatrix a, BLAS3 e) => e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e    
sapplyMat k a b = checkMatMatMult (shape a) (shape b) $ unsafeSApplyMat k a b
{-# INLINE sapplyMat #-}

unsafeApply :: (IMatrix a, BLAS3 e) => a (m,n) e -> Vector n e -> Vector m e
unsafeApply = unsafeSApply 1
{-# INLINE unsafeApply #-}

unsafeApplyMat :: (IMatrix a, BLAS3 e) => a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
unsafeApplyMat = unsafeSApplyMat 1
{-# INLINE unsafeApplyMat #-}

{-# RULES
"scale.apply/sapply"       forall k a x. a <*>  (k *> x) = sapply k a x
"scale.applyMat/sapplyMat" forall k a b. a <**> (k *> b) = sapplyMat k a b
  #-}
