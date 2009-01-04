{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Tensor.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Overloaded interface for mutable and immutable tensors.  This module
-- contains the common functionality for the classes in 
-- "Data.Tensor.Class.ITensor" and "Data.Tensor.Class.MTensor".
--

module Data.Tensor.Class (
    Shaped(..),
    ) where

import Data.Ix

-- | The base class for objects with shapes and indices 
-- (i.e. Vector, Matrix, etc.).
class (Ix i, Show i) => Shaped x i e | x -> i where
    -- | Get the shape of the tensor.  For vectors this is the dimension.
    -- For matrices, this will be a pair @(m,n)@ of the number of rows
    -- and columns.
    shape :: x n e -> i
    
    -- | Get the range of valid indices in the tensor.
    bounds :: x n e -> (i,i)
