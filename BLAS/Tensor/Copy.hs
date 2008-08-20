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
    copyTensor,
    ) where

import BLAS.Internal( checkBinaryOp )
import BLAS.Tensor.Write

-- | Class for mutable dense read-only tensors.
class (ReadTensor x i e m, WriteTensor y i e m) => CopyTensor x y i e m | y -> m where
    -- | Create a copy of a tensor.
    newCopyTensor :: x n e -> m (y n e)
    
    -- | Same as 'copy' but does not check that the argument shapes match.
    unsafeCopyTensor :: y n e -> x n e -> m ()

-- | @copyTensor dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyTensor :: (CopyTensor x y i e m) => y n e -> x n e -> m ()
copyTensor y x = checkBinaryOp (shape x) (shape y) $ unsafeCopyTensor y x
{-# INLINE copyTensor #-}
