{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Scalable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Scalable (
    Scalable(..)
    ) where

infixl 7 *>

-- | A class for scalable tensors.
class (Num e) => Scalable x e where
    -- | Scale a tensor by the given value.
    (*>) :: e -> x e -> x e
