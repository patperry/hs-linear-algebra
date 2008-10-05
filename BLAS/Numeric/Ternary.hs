{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Numeric.Ternary
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Numeric.Ternary (
    -- * Elementwise operations with specified destination.    
    doAdd,
    doSub,
    doMul,
    doDiv,
    
    -- * Elementwise operations that create new objects to store the result.    
    getAdd,
    getSub,
    getMul,
    getDiv,

    -- * The Numeric3 type class
    Numeric3(..),
    ) where

import BLAS.Internal ( checkBinaryOp, checkTernaryOp )
import BLAS.Tensor
import BLAS.Numeric.Read
import BLAS.Numeric.Write

-- | Class for mutable dense read-only tensors.
class (ReadNumeric x i e m, ReadNumeric y i e m, WriteNumeric z i e m) => Numeric3 x y z i e m where
    -- | Same as 'doAdd' but does not check that the arguments are the same
    -- shape.
    unsafeDoAdd :: x n e -> y n e -> z n e -> m ()

    -- | Same as 'doSub' but does not check that the arguments are the same
    -- shape.
    unsafeDoSub :: x n e -> y n e -> z n e -> m ()

    -- | Same as 'doMul' but does not check that the arguments are the same
    -- shape.
    unsafeDoMul :: x n e -> y n e -> z n e -> m ()

    -- | Same as 'doDiv' but does not check that the arguments are the same
    -- shape.
    unsafeDoDiv :: x n e -> y n e -> z n e -> m ()

    -- | Same as 'getAdd' but does not check that the arguments are the same 
    -- shape.
    unsafeGetAdd :: x n e -> y n e -> m (z n e)
    unsafeGetAdd = unsafeGetBinaryOp unsafeDoAdd
    {-# INLINE unsafeGetAdd #-}

    -- | Same as 'getSub' but does not check that the arguments are the same
    -- shape.
    unsafeGetSub :: x n e -> y n e -> m (z n e)
    unsafeGetSub = unsafeGetBinaryOp unsafeDoSub
    {-# INLINE unsafeGetSub #-}

    -- | Same as 'getMul' but does not check that the arguments are the same
    -- shape.
    unsafeGetMul :: x n e -> y n e -> m (z n e)
    unsafeGetMul = unsafeGetBinaryOp unsafeDoMul
    {-# INLINE unsafeGetMul #-}

    -- | Same as 'getDiv' but does not check that the arguments are the same
    -- shape.
    unsafeGetDiv :: x n e -> y n e -> m (z n e)
    unsafeGetDiv = unsafeGetBinaryOp unsafeDoDiv
    {-# INLINE unsafeGetDiv #-}


unsafeGetBinaryOp :: (Numeric3 x y z i e m) => 
    (x n e -> y n e -> z n e -> m ()) ->
        x n e -> y n e -> m (z n e)
unsafeGetBinaryOp f x y = do
    z <- newZero (shape x)
    f x y z
    return z
    

-- | @doAdd x y z@ replaces @z@ with @x+y@.  The operands must be the
-- same shape.
doAdd :: (Numeric3 x y z i e m) => x n e -> y n e -> z n e -> m ()
doAdd = checkTensorOp3 unsafeDoAdd
{-# INLINE doAdd #-}

-- | @doSub x y z@ replaces @z@ with @x-y@.  The operands must be the
-- same shape.
doSub :: (Numeric3 x y z i e m) => x n e -> y n e -> z n e -> m ()
doSub = checkTensorOp3 unsafeDoSub
{-# INLINE doSub #-}

-- | @doMul x y z@ replaces @z@ with @x*y@.  The operands must be the
-- same shape.
doMul :: (Numeric3 x y z i e m) => x n e -> y n e -> z n e -> m ()
doMul = checkTensorOp3 unsafeDoMul
{-# INLINE doMul #-}

-- | @doDiv x y z@ replaces @z@ with @x/y@.  The operands must be the
-- same shape.
doDiv :: (Numeric3 x y z i e m) => x n e -> y n e -> z n e -> m ()
doDiv = checkTensorOp3 unsafeDoDiv
{-# INLINE doDiv #-}

-- | @getAdd x y@ creates a new tensor equal to the sum @x+y@.  The operands
-- must have the same shape.
getAdd :: (Numeric3 x y z i e m) => x n e -> y n e -> m (z n e)
getAdd = checkTensorOp2 unsafeGetAdd
{-# INLINE getAdd #-}

-- | @getSub x y@ creates a new tensor equal to the difference @x-y@.  The 
-- operands must have the same shape.
getSub :: (Numeric3 x y z i e m) => x n e -> y n e -> m (z n e)
getSub = checkTensorOp2 unsafeGetSub
{-# INLINE getSub #-}

-- | @getMul x y@ creates a new tensor equal to the product @x*y@.  The 
-- operands must have the same shape.
getMul :: (Numeric3 x y z i e m) => x n e -> y n e -> m (z n e)
getMul = checkTensorOp2 unsafeGetMul
{-# INLINE getMul #-}

-- | @getDiv x y@ creates a new tensor equal to the ratio @x/y@.  The 
-- operands must have the same shape.
getDiv :: (Numeric3 x y z i e m) => x n e -> y n e -> m (z n e)
getDiv = checkTensorOp2 unsafeGetDiv
{-# INLINE getDiv #-}

checkTensorOp2 :: (BaseTensor x i e, BaseTensor y i e) => 
    (x n e -> y n e -> a) ->
        x n e -> y n e -> a
checkTensorOp2 f x y = 
    checkBinaryOp (shape x) (shape y) $ f x y
{-# INLINE checkTensorOp2 #-}

checkTensorOp3 :: (BaseTensor x i e, BaseTensor y i e, BaseTensor z i e) => 
    (x n e -> y n e -> z n e -> a) ->
    x n e -> y n e -> z n e -> a
checkTensorOp3 f x y z = 
    checkTernaryOp (shape x) (shape y) (shape z) $ f x y z
{-# INLINE checkTensorOp3 #-}
