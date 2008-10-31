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
    -- ** Unary
    getConjVector,
    getScaledVector,
    getShiftedVector,
    
    -- ** Binary
    getAddVector,
    getSubVector,
    getMulVector,
    getDivVector,
    addVector,
    subVector,
    axpyVector,
    mulVector,
    divVector,

    -- ** Unsafe
    unsafeGetAddVector,
    unsafeGetSubVector,
    unsafeGetMulVector,
    unsafeGetDivVector,
    unsafeAddVector,
    unsafeSubVector,
    unsafeAxpyVector,
    unsafeMulVector,
    unsafeDivVector,
    ) where

import BLAS.Elem( BLAS1 )
import BLAS.Internal( checkBinaryOp )
import BLAS.Tensor( BaseTensor(..) )
import Data.Vector.Dense.Class.Internal

---------------------------- Unary Operations -----------------------------

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


---------------------------- Binary Operations -----------------------------

-- | @getAddVector x y@ creates a new vector equal to the sum @x+y@.  The 
-- operands must have the same dimension.
getAddVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m, BLAS1 e) => 
    x n e -> y n e -> m (z n e)
getAddVector = checkTensorOp2 unsafeGetAddVector
{-# INLINE getAddVector #-}

unsafeGetAddVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m, BLAS1 e) => 
    x n e -> y n e -> m (z n e)
unsafeGetAddVector = unsafeGetBinaryOp unsafeAddVector
{-# INLINE unsafeGetAddVector #-}

-- | @getSubVector x y@ creates a new tensor equal to the difference @x-y@.  
-- The operands must have the same dimension.
getSubVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m, BLAS1 e) => 
    x n e -> y n e -> m (z n e)
getSubVector = checkTensorOp2 unsafeGetSubVector
{-# INLINE getSubVector #-}

unsafeGetSubVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m, BLAS1 e) => 
    x n e -> y n e -> m (z n e)
unsafeGetSubVector = unsafeGetBinaryOp unsafeSubVector
{-# INLINE unsafeGetSubVector #-}

-- | @getMulVector x y@ creates a new vector equal to the elementwise product 
-- @x*y@.  The operands must have the same dimensino
getMulVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m, BLAS1 e) => 
    x n e -> y n e -> m (z n e)
getMulVector = checkTensorOp2 unsafeGetMulVector
{-# INLINE getMulVector #-}

unsafeGetMulVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m, BLAS1 e) => 
    x n e -> y n e -> m (z n e)
unsafeGetMulVector = unsafeGetBinaryOp unsafeMulVector
{-# INLINE unsafeGetMulVector #-}

-- | @getDivVector x y@ creates a new vector equal to the elementwise 
-- ratio @x/y@.  The operands must have the same shape.
getDivVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m, BLAS1 e) => 
    x n e -> y n e -> m (z n e)
getDivVector = checkTensorOp2 unsafeGetDivVector
{-# INLINE getDivVector #-}

unsafeGetDivVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m, BLAS1 e) => 
    x n e -> y n e -> m (z n e)
unsafeGetDivVector = unsafeGetBinaryOp unsafeDivVector
{-# INLINE unsafeGetDivVector #-}

-- | @axpyVector alpha x y@ replaces @y@ with @alpha * x + y@.
axpyVector :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    e -> x n e -> y n e -> m ()
axpyVector alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyVector alpha x y
{-# INLINE axpyVector #-}

-- | @addVector y x@ replaces @y@ with @y+x@.
addVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
addVector y x = checkBinaryOp (dim y) (dim x) $ unsafeAddVector y x
{-# INLINE addVector #-}

unsafeAddVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
unsafeAddVector y x = unsafeAxpyVector 1 x y

-- | @subVector y x@ replaces @y@ with @y-x@.
subVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
subVector y x = checkBinaryOp (dim y) (dim x) $ unsafeSubVector y x
{-# INLINE subVector #-}

unsafeSubVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
unsafeSubVector y x = unsafeAxpyVector (-1) x y

-- | @mulVector y x@ replaces @y@ with @y*x@.
mulVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
mulVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeMulVector y x
{-# INLINE mulVector #-}

-- | @divVector y x@ replaces @y@ with @y/x@.
divVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
divVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeDivVector y x
{-# INLINE divVector #-}


checkTensorOp2 :: (BaseTensor x i, BaseTensor y i) => 
    (x n e -> y n e -> a) ->
        x n e -> y n e -> a
checkTensorOp2 f x y = 
    checkBinaryOp (shape x) (shape y) $ f x y
{-# INLINE checkTensorOp2 #-}

getUnaryOp :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    (y n e -> m ()) -> x n e -> m (y n e)
getUnaryOp f x = do
    y <- newCopyVector x
    f y
    return y
{-# INLINE getUnaryOp #-}

unsafeGetBinaryOp :: 
    (WriteVector z e m, ReadVector x e m, ReadVector y e m, BLAS1 e) => 
    (z n e -> y n e -> m ()) ->
        x n e -> y n e -> m (z n e)
unsafeGetBinaryOp f x y = do
    z <- newCopyVector x
    f z y
    return z
