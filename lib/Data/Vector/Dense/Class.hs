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
    module Data.Elem.Conj,
    
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

    -- * Conversions between mutable and immutable vectors
    UnsafeFreezeVector(..),
    UnsafeThawVector(..),
    freezeVector,
    thawVector,
    
    ) where

import Data.Vector.Dense.Class.Internal( IOVector, BaseVector(..), ReadVector, 
    WriteVector, dim, stride, isConj, coerceVector, withVectorPtr )
import Data.Tensor.Class
import Data.Elem.Conj
import Data.Vector.Dense.Class.Creating
import Data.Vector.Dense.Class.Elements
import Data.Vector.Dense.Class.Special
import Data.Vector.Dense.Class.Views
import Data.Vector.Dense.Class.Copying
import Data.Vector.Dense.Class.Properties
import Data.Vector.Dense.Class.Operations

import Data.Vector.Dense.Internal hiding ( V )
import qualified Data.Vector.Dense.Internal as I
import Data.Vector.Dense.Internal (  )
import Data.Vector.Dense.STBase

class UnsafeFreezeVector x where
    unsafeFreezeVector :: x n e -> Vector n e
instance UnsafeFreezeVector IOVector where
    unsafeFreezeVector = I.V
instance UnsafeFreezeVector (STVector s) where
    unsafeFreezeVector = unsafeFreezeVector . unsafeSTVectorToIOVector

class UnsafeThawVector x where
    unsafeThawVector :: Vector n e -> x n e
instance UnsafeThawVector IOVector where
    unsafeThawVector (I.V x) = x
instance UnsafeThawVector (STVector s) where
    unsafeThawVector = unsafeIOVectorToSTVector . unsafeThawVector
    
freezeVector :: (ReadVector x e m, WriteVector y e m, UnsafeFreezeVector y) =>
    x n e -> m (Vector n e)
freezeVector x = do
    x' <- newCopyVector x
    return (unsafeFreezeVector x')

thawVector :: (WriteVector y e m) =>
    Vector n e -> m (y n e)
thawVector = newCopyVector
