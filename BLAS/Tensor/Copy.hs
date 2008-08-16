{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Copy
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Copy (
    CopyTensor(..),
    copy,
    ) where

import BLAS.Internal( checkBinaryOp )
import BLAS.Tensor.Write

-- | Class for mutable dense read-only tensors.
class (ReadTensor x i e m, WriteTensor y i e m) => CopyTensor x y i e m | m -> y where
    -- | Create a copy of a tensor.
    newCopy :: x n e -> m (y n e)
    
    -- | Same as 'copy' but does not check that the argument shapes match.
    unsafeCopy :: y n e -> x n e -> m ()

-- | @copy dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copy :: (CopyTensor x y i e m) => y n e -> x n e -> m ()
copy y x = checkBinaryOp (shape x) (shape y) $ unsafeCopy y x
{-# INLINE copy #-}
