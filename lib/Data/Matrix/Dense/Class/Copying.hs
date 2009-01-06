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
    
    -- * Swapping rows and columns
    swapRows,
    swapCols,
    unsafeSwapRows,
    unsafeSwapCols,
    
    ) where

import BLAS.Internal( checkBinaryOp )

import Control.Monad( when )

import Data.Tensor.Class( shape )
import Data.Matrix.Dense.Class.Internal
import Data.Matrix.Dense.Class.Views
import Data.Vector.Dense.Base


-- | @copyMatrix dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyMatrix :: (WriteMatrix b e m, ReadMatrix a e m) => 
    b mn e -> a mn e -> m ()
copyMatrix b a = checkBinaryOp (shape b) (shape a) $ unsafeCopyMatrix b a
{-# INLINE copyMatrix #-}

-- | @swapMatrix x y@ swaps the values stored in two matrices.
swapMatrix :: (WriteMatrix a e m) => 
    a mn e -> a mn e -> m ()
swapMatrix a b = checkBinaryOp (shape b) (shape a) $ unsafeSwapMatrix a b
{-# INLINE swapMatrix #-}

swapRows :: (WriteMatrix a e m) => a (r,s) e -> Int -> Int -> m ()
swapRows a i j = 
    when (i /= j) $ unsafeSwapVector (rowView a i) (rowView a j)

swapCols :: (WriteMatrix a e m) => a (r,s) e -> Int -> Int -> m ()
swapCols a i j = 
    when (i /= j) $ unsafeSwapVector (colView a i) (colView a j)

unsafeSwapRows :: (WriteMatrix a e m) => a (r,s) e -> Int -> Int -> m ()
unsafeSwapRows a i j = 
    when (i /= j) $ unsafeSwapVector (unsafeRowView a i) (unsafeRowView a j)

unsafeSwapCols :: (WriteMatrix a e m) => a (r,s) e -> Int -> Int -> m ()
unsafeSwapCols a i j = 
    when (i /= j) $ unsafeSwapVector (unsafeColView a i) (unsafeColView a j)
