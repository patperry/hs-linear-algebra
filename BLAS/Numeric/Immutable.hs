{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Numeric.Immutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Numeric.Immutable (
    INumeric(..)
    ) where

import BLAS.Tensor.Immutable

infixl 7 *>
infixl 5 `shift`

class (ITensor x i e) => INumeric x i e where
    -- | Scale every element by the given value.
    (*>) :: e -> x n e -> x n e
    
    -- | Add a constant to every element.
    shift :: e -> x n e -> x n e
