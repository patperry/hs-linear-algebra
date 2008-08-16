{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Base (
    BaseVector(..),
    
    module BLAS.Conj,
    
    dim,
    subvector,
    subvectorWithStride,
    coerceVector,
    ) where

import BLAS.Internal ( checkedSubvector, checkedSubvectorWithStride )
import BLAS.Conj
import BLAS.Tensor
import Foreign ( Ptr )
import Unsafe.Coerce

class (BaseTensor x Int e) => BaseVector x e where
    stride :: x n e -> Int
    isConj :: x n e -> Bool
    
    conjVector :: x n e -> x n e
    
    unsafeSubvector :: x n e -> Int -> Int -> x n' e
    unsafeSubvector = unsafeSubvectorWithStride 1
    
    unsafeSubvectorWithStride :: Int -> x n e -> Int -> Int -> x n' e
    
    withVectorPtr :: x n e -> (Ptr e -> IO a) -> IO a

instance (BaseVector x e) => Conj (x n e) where
    conj = conjVector
    
-- | Get the dimension of the vector.
dim :: (BaseVector x e) => x n e -> Int
dim = shape
{-# INLINE dim #-}
    
-- | @subvector x o n@ creates a subvector view of @x@ starting at index @o@ 
-- and having length @n@.
subvector :: (BaseVector x e) => x n e -> Int -> Int -> x n' e
subvector x = checkedSubvector (dim x) (unsafeSubvector x)

-- | @subvectorWithStride s x o n@ creates a subvector view of @x@ starting 
-- at index @o@, having length @n@ and stride @s@.
subvectorWithStride :: (BaseVector x e) => Int -> x n e -> Int -> Int -> x n' e
subvectorWithStride s x = 
    checkedSubvectorWithStride s (dim x) (unsafeSubvectorWithStride s x)

-- | Cast the shape type of the vector.
coerceVector :: (BaseVector x e) => x n e -> x n' e
coerceVector = unsafeCoerce
{-# INLINE coerceVector #-}
