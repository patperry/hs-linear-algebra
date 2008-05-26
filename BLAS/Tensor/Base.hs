{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Base (
    Tensor(..),
    ) where

import Data.Ix

infixl *>

-- | The base class for tensors (i.e. Vector, Matrix, etc.).
class (Ix i) => Tensor x i e | x -> i where
    -- | Get the shape of the tensor.  For vectors this is the dimension.
    -- For matrices, this will be a pair @(m,n)@ of the number of rows
    -- and columns.
    shape :: x e -> i
    
    -- | Get the range of valid indices in the tensor.
    bounds :: x e -> (i,i)
