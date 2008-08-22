{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Swap
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Swap (
    swapTensor,
    SwapTensor(..),
    ) where

import BLAS.Tensor.Base
import BLAS.Tensor.Write

class (WriteTensor x i e m) => SwapTensor x i e m where
    -- | Same as 'swapTensor' but arguments are not range-checed.
    unsafeSwapTensor :: x n e -> x n e -> m ()


-- | Swap the values stored in two tensors of the same shape.
swapTensor :: (SwapTensor x i e m) => x n e -> x n e -> m ()
swapTensor x y
    | n1 /= n2 =
        fail $ "tried to swap tensors of shapes `" ++ show n1 ++ "'"
               ++ " and `" ++ show n2 ++ "'"
    | otherwise =
        unsafeSwapTensor x y
  where
    n1 = shape x
    n2 = shape y
