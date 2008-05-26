{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor (
    Tensor(..),
    ) where

import Data.Ix

infixl *>

-- | A class of mutable or immutable tensor (i.e. Vector, Matrix, etc.)
-- having a shape and index bounds.
class (Ix i) => Tensor x i e | x -> i where
    -- | Get the shape of the tensor.  For vectors this is the dimension.
    -- For matrices, this will be a pair @(m,n)@ of the number of rows
    -- and columns.
    shape :: x e -> i
    
    -- | Get the range of valid indices in the tensor.
    bounds :: x e -> (i,i)
