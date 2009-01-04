-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Operations
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Operations (
    -- * Matrix operations
    -- ** Multiplication
    module Data.Matrix.Class.MMatrix,
    
    -- ** Unary
    getConjMatrix,
    getScaledMatrix,
    getShiftedMatrix,
    
    -- ** Binary
    getAddMatrix,
    getSubMatrix,
    getMulMatrix,
    getDivMatrix,
    addMatrix,
    subMatrix,
    axpyMatrix,
    mulMatrix,
    divMatrix,

    -- ** Unsafe
    unsafeGetAddMatrix,
    unsafeGetSubMatrix,
    unsafeGetMulMatrix,
    unsafeGetDivMatrix,
    unsafeAddMatrix,
    unsafeSubMatrix,
    unsafeAxpyMatrix,
    unsafeMulMatrix,
    unsafeDivMatrix,
    
    ) where

import BLAS.Internal( checkBinaryOp )
import BLAS.Tensor( BaseTensor(..) )

import Data.Matrix.Dense.Class.Internal

import Data.Matrix.Class.MMatrix

---------------------------- Unary Operations -----------------------------

-- | Get a new matrix with elements with the conjugates of the elements
-- of the given matrix.
getConjMatrix :: (ReadMatrix a e m, WriteMatrix b e m) =>
    a mn e -> m (b mn e)
getConjMatrix = getUnaryOp doConjMatrix
{-# INLINE getConjMatrix #-}

-- | Get a new matrix by scaling the elements of another matrix
-- by a given value.
getScaledMatrix :: (ReadMatrix a e m, WriteMatrix b e m) =>
    e -> a mn e -> m (b mn e)
getScaledMatrix e = getUnaryOp (scaleByMatrix e)
{-# INLINE getScaledMatrix #-}

-- | Get a new matrix by shifting the elements of another matrix
-- by a given value.
getShiftedMatrix :: (ReadMatrix a e m, WriteMatrix b e m) =>
    e -> a mn e -> m (b mn e)
getShiftedMatrix e = getUnaryOp (shiftByMatrix e)
{-# INLINE getShiftedMatrix #-}


---------------------------- Binary Operations -----------------------------


-- | @getAddMatrix a b@ creates a new matrix equal to the sum @a+b@.  The 
-- operands must have the same shape.
getAddMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    a mn e -> b mn e -> m (c mn e)
getAddMatrix = checkTensorOp2 unsafeGetAddMatrix
{-# INLINE getAddMatrix #-}

unsafeGetAddMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    a mn e -> b mn e -> m (c mn e)
unsafeGetAddMatrix = unsafeGetBinaryOp unsafeAddMatrix
{-# INLINE unsafeGetAddMatrix #-}

-- | @getSubMatrix a b@ creates a new matrix equal to the difference @a-b@.  The 
-- operands must have the same shape.
getSubMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    a mn e -> b mn e -> m (c mn e)
getSubMatrix = checkTensorOp2 unsafeGetSubMatrix
{-# INLINE getSubMatrix #-}

unsafeGetSubMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    a mn e -> b mn e -> m (c mn e)
unsafeGetSubMatrix = unsafeGetBinaryOp unsafeSubMatrix
{-# INLINE unsafeGetSubMatrix #-}

-- | @getMulMatrix a b@ creates a new matrix equal to the elementwise product 
-- @a*b@.  The operands must have the same shape.
getMulMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    a mn e -> b mn e -> m (c mn e)
getMulMatrix = checkTensorOp2 unsafeGetMulMatrix
{-# INLINE getMulMatrix #-}

unsafeGetMulMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    a mn e -> b mn e -> m (c mn e)
unsafeGetMulMatrix = unsafeGetBinaryOp unsafeMulMatrix
{-# INLINE unsafeGetMulMatrix #-}

-- | @getDivMatrix a b@ creates a new matrix equal to the elementwise ratio
-- @a/b@.  The operands must have the same shape.
getDivMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    a mn e -> b mn e -> m (c mn e)
getDivMatrix = checkTensorOp2 unsafeGetDivMatrix
{-# INLINE getDivMatrix #-}

unsafeGetDivMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    a mn e -> b mn e -> m (c mn e)
unsafeGetDivMatrix = unsafeGetBinaryOp unsafeDivMatrix
{-# INLINE unsafeGetDivMatrix #-}


axpyMatrix :: (ReadMatrix a e m, WriteMatrix b e m) =>
    e -> a n e -> b n e -> m ()
axpyMatrix alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyMatrix alpha x y
{-# INLINE axpyMatrix #-}

addMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b n e -> a n e -> m ()
addMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeAddMatrix b a
{-# INLINE addMatrix #-}

unsafeAddMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b n e -> a n e -> m ()
unsafeAddMatrix b a = unsafeAxpyMatrix 1 a b

subMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b n e -> a n e -> m ()
subMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeSubMatrix b a
{-# INLINE subMatrix #-}

unsafeSubMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b n e -> a n e -> m ()
unsafeSubMatrix b a = unsafeAxpyMatrix (-1) a b

mulMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b n e -> a n e -> m ()
mulMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeMulMatrix b a
{-# INLINE mulMatrix #-}

divMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b n e -> a n e -> m ()
divMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeDivMatrix b a
{-# INLINE divMatrix #-}


checkTensorOp2 :: (BaseTensor x i e, BaseTensor y i e) => 
    (x n e -> y n e -> a) ->
        x n e -> y n e -> a
checkTensorOp2 f x y = 
    checkBinaryOp (shape x) (shape y) $ f x y
{-# INLINE checkTensorOp2 #-}

getUnaryOp :: (ReadMatrix a e m, WriteMatrix b e m) =>
    (b mn e -> m ()) -> a mn e -> m (b mn e)
getUnaryOp f a = do
    b <- newCopyMatrix a
    f b
    return b
{-# INLINE getUnaryOp #-}

unsafeGetBinaryOp :: 
    (WriteMatrix c e m, ReadMatrix a e m, ReadMatrix b f m) => 
    (c n e -> b n f -> m ()) ->
        a n e -> b n f -> m (c n e)
unsafeGetBinaryOp f a b = do
    c <- newCopyMatrix a
    f c b
    return c
