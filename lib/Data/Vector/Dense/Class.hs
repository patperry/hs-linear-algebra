-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class (
    -- * The dense vector type classes
    BaseVector(..),
    ReadVector,
    WriteVector,
    
    -- * Vector shape
    dim,
    coerceVector,
    module Data.Tensor.Class,

    -- * Conjugating vectors
    module BLAS.Conj,
    
    module Data.Vector.Dense.Class.Creating,
    module Data.Vector.Dense.Class.Elements,
    module Data.Vector.Dense.Class.Special,
    module Data.Vector.Dense.Class.Views,
    module Data.Vector.Dense.Class.Copying,
    module Data.Vector.Dense.Class.Properties,
    module Data.Vector.Dense.Class.Operations,
    
    -- * Low-level functions
    stride,
    isConj,
    withVectorPtr,
    
    ) where

import Data.Vector.Dense.Class.Internal( BaseVector(..), dim, stride, isConj,
    ReadVector, WriteVector, coerceVector, withVectorPtr )
import Data.Tensor.Class
import BLAS.Conj
import Data.Vector.Dense.Class.Creating
import Data.Vector.Dense.Class.Elements
import Data.Vector.Dense.Class.Special
import Data.Vector.Dense.Class.Views
import Data.Vector.Dense.Class.Copying
import Data.Vector.Dense.Class.Properties
import Data.Vector.Dense.Class.Operations
