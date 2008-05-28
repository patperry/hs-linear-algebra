{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Dense.ReadOnly
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Dense.ReadOnly (
    RDTensor(..),
    ) where

import BLAS.Tensor.ReadOnly

-- | Class for mutable dense read-only tensors.
class (RTensor x i e m) => RDTensor x i e m where
    -- | Creates a new tensor with elements all initialized to zero.
    newZero :: i -> m (x e)
    
    -- | Creates a new tensor with elements all initialized to the 
    -- given value.
    newConstant :: i -> e -> m (x e)
    