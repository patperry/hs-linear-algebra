-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Numeric.Binary
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Numeric.Binary (
    -- * In-place binary operations
    axpy,
    (+=),
    (-=),
    (*=),
    (//=),
    
    -- * The Numeric2 type class
    Numeric2(..),
    ) where

import BLAS.Internal ( checkBinaryOp )
import BLAS.Tensor

infixl 1 +=, -=, *=, //=

-- | Class for mutable dense read-only tensors.
class (ReadTensor x i e m, WriteTensor y i e m) => Numeric2 x y i e m where
    -- | Same as 'axpy' but does not check that the argument shapes match.    
    unsafeAxpy :: e -> x n e -> y n e -> m ()

    -- | Same as '(+=)' but does not check that the argument shapes match.    
    unsafeAdd :: y n e -> x n e -> m ()
    unsafeAdd y x = unsafeAxpy 1 x y

    -- | Same as '(-=)' but does not check that the argument shapes match.    
    unsafeSub :: y n e -> x n e -> m ()
    unsafeSub y x = unsafeAxpy (-1) x y

    -- | Same as '(*=)' but does not check that the argument shapes match.
    unsafeMul :: y n e -> x n e -> m ()

    -- | Same as '(//=)' but does not check that the argument shapes match.            
    unsafeDiv :: y n e -> x n e -> m ()

-- | @axpy alpha x y@ replaces @y@ with @alpha x + y@.  The operands
-- must be the same shape.
axpy :: (Numeric2 x y i e m) => e -> x n e -> y n e -> m ()
axpy k = checkTensorOp2 $ unsafeAxpy k
{-# INLINE axpy #-}

-- | @y += x@ replaces @y@ with @y + x@.  The operands
-- must be the same shape.
(+=) :: (Numeric2 x y i e m) => y n e -> x n e -> m ()
(+=) = checkTensorOp2 unsafeAdd
{-# INLINE (+=) #-}

-- | @y -= x@ replaces @y@ with @y - x@.  The operands
-- must be the same shape.
(-=) :: (Numeric2 x y i e m) => y n e -> x n e -> m ()
(-=) = checkTensorOp2 unsafeSub
{-# INLINE (-=) #-}

-- | @y *= x@ replaces @y@ with @y * x@.  The operands
-- must be the same shape.
(*=) :: (Numeric2 x y i e m) => y n e -> x n e -> m ()
(*=) = checkTensorOp2 unsafeMul
{-# INLINE (*=) #-}

-- | @y //= x@ replaces @y@ with @y / x@.  The operands
-- must be the same shape.
(//=) :: (Numeric2 x y i e m) => y n e -> x n e -> m ()
(//=) = checkTensorOp2 unsafeDiv
{-# INLINE (//=) #-}


checkTensorOp2 :: (BaseTensor x i e, BaseTensor y i e) => 
    (x n e -> y n e -> a) ->
        x n e -> y n e -> a
checkTensorOp2 f x y = 
    checkBinaryOp (shape x) (shape y) $ f x y
{-# INLINE checkTensorOp2 #-}
