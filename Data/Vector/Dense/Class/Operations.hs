-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Operations
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Operations (
    -- * Vector operations
    getConjVector,
    getScaledVector,
    getShiftedVector,
    
    axpyVector,
    mulVector,
    divVector,
    
    unsafeAxpyVector,
    unsafeMulVector,
    unsafeDivVector,
    
    module BLAS.Numeric,
    ) where

import BLAS.Elem( BLAS1 )
import BLAS.Internal( checkBinaryOp )
import BLAS.Numeric
import BLAS.Tensor( shape )
import Data.Vector.Dense.Class.Internal

-- | Get a new vector with elements with the conjugates of the elements
-- of the given vector
getConjVector :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    x n e -> m (y n e)
getConjVector = getUnaryOp doConjVector
{-# INLINE getConjVector #-}

-- | Get a new vector by scaling the elements of another vector
-- by a given value.
getScaledVector :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    e -> x n e -> m (y n e)
getScaledVector e = getUnaryOp (scaleByVector e)
{-# INLINE getScaledVector #-}

-- | Get a new vector by shifting the elements of another vector
-- by a given value.
getShiftedVector :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    e -> x n e -> m (y n e)
getShiftedVector e = getUnaryOp (shiftByVector e)
{-# INLINE getShiftedVector #-}



axpyVector :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    e -> x n e -> y n e -> m ()
axpyVector alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyVector alpha x y
{-# INLINE axpyVector #-}

mulVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
mulVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeMulVector y x
{-# INLINE mulVector #-}

divVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
divVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeDivVector y x
{-# INLINE divVector #-}


getUnaryOp :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    (y n e -> m ()) -> x n e -> m (y n e)
getUnaryOp f x = do
    y <- newCopyVector x
    f y
    return y
{-# INLINE getUnaryOp #-}

