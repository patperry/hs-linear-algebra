{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Tensor.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Tensor.Class (
    BaseTensor(..),
    ) where

import Data.Ix

-- | The base class for tensors (i.e. Vector, Matrix, etc.).
class (Ix i, Eq i, Show i) => BaseTensor x i e | x -> i where
    -- | Get the shape of the tensor.  For vectors this is the dimension.
    -- For matrices, this will be a pair @(m,n)@ of the number of rows
    -- and columns.
    shape :: x n e -> i
    
    -- | Get the range of valid indices in the tensor.
    bounds :: x n e -> (i,i)