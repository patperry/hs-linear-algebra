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

    -- * Operators
    (<*>),
    (<**>),

    -- * Rows and columns
    row,
    col,
    rows,
    cols,

    -- * Multiplication
    applyVector,
    applyMatrix,
    sapplyVector,
    sapplyMatrix,
    
    -- * Unsafe operations
    unsafeApplyVector,
    unsafeApplyMatrix,
    ) where

import Data.Elem.BLAS
import BLAS.Internal ( checkedRow, checkedCol, checkMatVecMult, 
    checkMatMatMult )

import Data.Tensor.Class

import Data.Vector.Dense
import Data.Matrix.Dense.Base

infixr 7 <*>, <**>

-- | Get the given row in a matrix.
row :: (IMatrix a, Elem e) => a e -> Int -> Vector e
row a = checkedRow (shape a) (unsafeRow a)
{-# INLINE row #-}

-- | Get the given column in a matrix.
col :: (IMatrix a, Elem e) => a e -> Int -> Vector e
col a = checkedCol (shape a) (unsafeCol a)
{-# INLINE col #-}

-- | Matrix multiplication by a vector.
applyVector :: (IMatrix a, BLAS2 e) => a e -> Vector e -> Vector e
applyVector a x = checkMatVecMult (shape a) (dim x) $ unsafeApplyVector a x
{-# INLINE applyVector #-}

-- | Matrix multiplication by a matrix.
applyMatrix :: (IMatrix a, BLAS3 e) => a e -> Matrix e -> Matrix e
applyMatrix a b = checkMatMatMult (shape a) (shape b) $ unsafeApplyMatrix a b
{-# INLINE applyMatrix #-}

-- | Scale and multiply by a vector.  
-- @sapplyVector k a x@ is equal to @a \<*> (k *> x)@, and often it is faster.
sapplyVector :: (IMatrix a, BLAS2 e) => e -> a e -> Vector e -> Vector e
sapplyVector k a x = checkMatVecMult (shape a) (dim x) $ unsafeSApplyVector k a x
{-# INLINE sapplyVector #-}
    
-- | Scale and multiply by a matrix.
-- @sapplyMatrix k a b@ is equal to @a \<**> (k *> b)@, and often it is faster.
sapplyMatrix :: (IMatrix a, BLAS3 e) => e -> a e -> Matrix e -> Matrix e    
sapplyMatrix k a b = checkMatMatMult (shape a) (shape b) $ unsafeSApplyMatrix k a b
{-# INLINE sapplyMatrix #-}

unsafeApplyVector :: (IMatrix a, BLAS2 e) => a e -> Vector e -> Vector e
unsafeApplyVector = unsafeSApplyVector 1
{-# INLINE unsafeApplyVector #-}

unsafeApplyMatrix :: (IMatrix a, BLAS3 e) => a e -> Matrix e -> Matrix e
unsafeApplyMatrix = unsafeSApplyMatrix 1
{-# INLINE unsafeApplyMatrix #-}

-- | Operator form of matrix multiplication by a vector.
(<*>) :: (IMatrix a, BLAS2 e) => a e -> Vector e -> Vector e
(<*>) = applyVector
{-# INLINE (<*>) #-}

-- | Operator form of matrix multiplication by a matrix.
(<**>) :: (IMatrix a, BLAS3 e) => a e -> Matrix e -> Matrix e
(<**>) = applyMatrix
{-# INLINE (<**>) #-}

{-# RULES
"scale.apply/sapply" forall k a x. a <*> (k *> x) = sapplyVector k a x
"scale.applyMat/sapplyMat" forall k a b. a <**> (k *> b) = sapplyMatrix k a b
 #-}
