-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Copying
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Copying (
    -- * Copying vectors
    newCopyVector,
    copyVector,
    swapVector,
    unsafeCopyVector,
    unsafeSwapVector,
    

    ) where

import BLAS.Internal( checkBinaryOp )

import BLAS.Tensor

import Data.Vector.Dense.Class.Internal


-- | @copyVector dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyVector :: (WriteVector y e m, ReadVector x e m) =>
    y n e -> x n e -> m ()
copyVector y x = checkBinaryOp (shape x) (shape y) $ unsafeCopyVector y x
{-# INLINE copyVector #-}


-- | Swap the values stored in two vectors.
swapVector :: (WriteVector y e m) => 
    y n e -> y n e -> m ()
swapVector x y = checkBinaryOp (shape x) (shape y) $ unsafeSwapVector x y
{-# INLINE swapVector #-}

