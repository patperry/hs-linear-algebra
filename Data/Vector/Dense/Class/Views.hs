-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Views
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Views (
    -- * Vector views
    subvector,
    subvectorWithStride,
    unsafeSubvector,
    unsafeSubvectorWithStride,

    ) where

import BLAS.Internal ( checkedSubvector, checkedSubvectorWithStride )
import Data.Vector.Dense.Class.Internal


-- | @subvector x o n@ creates a subvector view of @x@ starting at index @o@ 
-- and having length @n@.
subvector :: (BaseVector x e) => x n e -> Int -> Int -> x n' e
subvector x = checkedSubvector (dim x) (unsafeSubvector x)
{-# INLINE subvector #-}

unsafeSubvector :: (BaseVector x e) => x n e -> Int -> Int -> x n' e
unsafeSubvector = unsafeSubvectorWithStride 1
{-# INLINE unsafeSubvector #-}

-- | @subvectorWithStride s x o n@ creates a subvector view of @x@ starting 
-- at index @o@, having length @n@ and stride @s@.
subvectorWithStride :: (BaseVector x e) => Int -> x n e -> Int -> Int -> x n' e
subvectorWithStride s x = 
    checkedSubvectorWithStride s (dim x) (unsafeSubvectorWithStride s x)
{-# INLINE subvectorWithStride #-}

unsafeSubvectorWithStride :: (BaseVector x e) => 
    Int -> x n e -> Int -> Int -> x n' e
unsafeSubvectorWithStride s' x o' n' =
    let (f,o,_,s,c) = arrayFromVector x
    in vectorViewArray f (o + s*o') n' (s*s') c
