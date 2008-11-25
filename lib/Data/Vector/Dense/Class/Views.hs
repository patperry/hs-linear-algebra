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
    subvectorView,
    subvectorViewWithStride,
    unsafeSubvectorView,
    unsafeSubvectorViewWithStride,

    ) where

import BLAS.Internal ( checkedSubvector, checkedSubvectorWithStride )
import Data.Vector.Dense.Class.Internal
import Foreign


-- | @subvectorView x o n@ creates a subvector view of @x@ starting at index @o@ 
-- and having length @n@.
subvectorView :: (BaseVector x e) => 
    x n e -> Int -> Int -> x n' e
subvectorView x = checkedSubvector (dim x) (unsafeSubvectorView x)
{-# INLINE subvectorView #-}

unsafeSubvectorView :: (BaseVector x e) => 
    x n e -> Int -> Int -> x n' e
unsafeSubvectorView = unsafeSubvectorViewWithStride 1
{-# INLINE unsafeSubvectorView #-}

-- | @subvectorViewWithStride s x o n@ creates a subvector view of @x@ starting 
-- at index @o@, having length @n@ and stride @s@.
subvectorViewWithStride :: (BaseVector x e) => 
    Int -> x n e -> Int -> Int -> x n' e
subvectorViewWithStride s x = 
    checkedSubvectorWithStride s (dim x) (unsafeSubvectorViewWithStride s x)
{-# INLINE subvectorViewWithStride #-}

unsafeSubvectorViewWithStride :: (BaseVector x e) => 
    Int -> x n e -> Int -> Int -> x n' e
unsafeSubvectorViewWithStride s' x o' n' =
    let (f,p,_,s,c) = arrayFromVector x
    in vectorViewArray f (p `advancePtr` (s*o')) n' (s*s') c
