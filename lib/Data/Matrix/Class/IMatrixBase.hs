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

import Data.Elem.BLAS( BLAS3 )
import BLAS.Internal ( checkedRow, checkedCol, checkMatVecMult, 
    checkMatMatMult )

import Data.Matrix.Class
import Data.Tensor.Class

import Data.Vector.Dense
import Data.Vector.Dense.ST( runSTVector )
import Data.Matrix.Herm
import Data.Matrix.TriBase
import Data.Matrix.Dense.Base
import Data.Matrix.Dense.ST( runSTMatrix )

infixr 7 <*>, <**>

-- | A type class for immutable matrices.  The member functions of the
-- type class do not perform any checks on the validity of shapes or
-- indices, so in general their safe counterparts should be preferred.
class (MatrixShaped a, BLAS3 e) => IMatrix a e where
    unsafeSApply :: e -> a (m,n) e -> Vector n e -> Vector m e
    unsafeSApplyMat :: e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e

    unsafeRow :: a (m,n) e -> Int -> Vector n e
    unsafeRow a i = let
        e = basisVector (numRows a) i
        in conj $ unsafeApply (herm a) e
    {-# INLINE unsafeRow #-}
    
    unsafeCol :: a (m,n) e -> Int -> Vector m e
    unsafeCol a j = let
        e = basisVector (numCols a) j
        in unsafeApply a e
    {-# INLINE unsafeCol #-}

-- | Get the given row in a matrix.
row :: (IMatrix a e) => a (m,n) e -> Int -> Vector n e
row a = checkedRow (shape a) (unsafeRow a)
{-# INLINE row #-}

-- | Get the given column in a matrix.
col :: (IMatrix a e) => a (m,n) e -> Int -> Vector m e
col a = checkedCol (shape a) (unsafeCol a)
{-# INLINE col #-}

-- | Get a list the row vectors in the matrix.
rows :: (IMatrix a e) => a (m,n) e -> [Vector n e]
rows a = [ unsafeRow a i | i <- [0..numRows a - 1] ]
{-# INLINE rows #-}

-- | Get a list the column vectors in the matrix.
cols :: (IMatrix a e) => a (m,n) e -> [Vector m e]
cols a = [ unsafeCol a j | j <- [0..numCols a - 1] ]
{-# INLINE cols #-}

-- | Matrix multiplication by a vector.
(<*>) :: (IMatrix a e) => a (m,n) e -> Vector n e -> Vector m e
(<*>) a x = checkMatVecMult (shape a) (dim x) $ unsafeApply a x
{-# INLINE (<*>) #-}

-- | Matrix multiplication by a matrix.
(<**>) :: (IMatrix a e) => a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
(<**>) a b = checkMatMatMult (shape a) (shape b) $ unsafeApplyMat a b
{-# INLINE (<**>) #-}

-- | Scale and multiply by a vector.  
-- @sapply k a x@ is equal to @a \<*> (k *> x)@, and often it is faster.
sapply :: (IMatrix a e) => e -> a (m,n) e -> Vector n e -> Vector m e
sapply k a x = checkMatVecMult (shape a) (dim x) $ unsafeSApply k a x
{-# INLINE sapply #-}
    
-- | Scale and multiply by a matrix.
-- @sapplyMat k a b@ is equal to @a \<**> (k *> b)@, and often it is faster.
sapplyMat :: (IMatrix a e) => e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e    
sapplyMat k a b = checkMatMatMult (shape a) (shape b) $ unsafeSApplyMat k a b
{-# INLINE sapplyMat #-}

unsafeApply :: (IMatrix a e) => a (m,n) e -> Vector n e -> Vector m e
unsafeApply = unsafeSApply 1
{-# INLINE unsafeApply #-}

unsafeApplyMat :: (IMatrix a e) => a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
unsafeApplyMat = unsafeSApplyMat 1
{-# INLINE unsafeApplyMat #-}

instance (BLAS3 e) => IMatrix Matrix e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    {-# INLINE unsafeSApply #-}    
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b
    {-# INLINE unsafeSApplyMat #-}   
    unsafeRow                 = unsafeRowView
    {-# INLINE unsafeRow #-}   
    unsafeCol                 = unsafeColView
    {-# INLINE unsafeCol #-}

instance (BLAS3 e) => IMatrix (Herm Matrix) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    {-# INLINE unsafeSApply #-}    
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    
    {-# INLINE unsafeSApplyMat #-}   
    unsafeRow a i = runSTVector $ unsafeGetRow a i
    {-# INLINE unsafeRow #-}
    unsafeCol a j = runSTVector $ unsafeGetCol a j
    {-# INLINE unsafeCol #-}

instance (BLAS3 e) => IMatrix (Tri Matrix) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    {-# INLINE unsafeSApply #-}
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b
    {-# INLINE unsafeSApplyMat #-}   
    unsafeRow a i = runSTVector $ unsafeGetRow a i
    {-# INLINE unsafeRow #-}
    unsafeCol a j = runSTVector $ unsafeGetCol a j
    {-# INLINE unsafeCol #-}

{-# RULES
"scale.apply/sapply"       forall k a x. a <*>  (k *> x) = sapply k a x
"scale.applyMat/sapplyMat" forall k a b. a <**> (k *> b) = sapplyMat k a b
  #-}
