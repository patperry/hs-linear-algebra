{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Dense.Immutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Dense.Immutable (
    IDTensor(..),
    ) where

import BLAS.Tensor.Immutable

-- | Class for immutable dense tensors.
class (ITensor x i e) => (IDTensor x i e) where
    -- | Get a zero tensor of the given shape.
    zero :: i -> x e
    
    -- | Get a new constant tensor of the given shape.
    constant :: i -> e -> x e
