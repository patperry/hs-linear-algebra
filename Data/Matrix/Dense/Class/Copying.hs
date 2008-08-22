-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Copying
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Copying (
    -- * Copying matrices
    newCopyMatrix,
    copyMatrix,
    swapMatrix,
    unsafeCopyMatrix,
    unsafeSwapMatrix,
    
    ) where

import BLAS.Elem
import BLAS.Internal( checkBinaryOp )
import BLAS.Tensor( shape )

import Data.Matrix.Dense.Class.Internal


-- | @copyMatrix dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyMatrix :: (BLAS1 e, WriteMatrix b y e m,  ReadMatrix a x e m) => 
    b mn e -> a mn e -> m ()
copyMatrix b a = checkBinaryOp (shape b) (shape a) $ unsafeCopyMatrix b a
{-# INLINE copyMatrix #-}

-- | @swapMatrix x y@ swaps the values stored in two matrices.
swapMatrix :: (BLAS1 e, WriteMatrix a x e m) => 
    a mn e -> a mn e -> m ()
swapMatrix a b = checkBinaryOp (shape b) (shape a) $ unsafeSwapMatrix a b
{-# INLINE swapMatrix #-}
