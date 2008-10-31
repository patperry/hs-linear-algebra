{-# OPTIONS_HADDOCK prune #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense (
    module Data.Vector.Dense.Internal,

    -- * Converting between mutable and immutable vectors
    UnsafeFreezeVector(..),
    UnsafeThawVector(..),
    freezeVector,
    thawVector,
    
    ) where

import BLAS.Elem
import Data.Vector.Dense.Internal hiding ( V )
import qualified Data.Vector.Dense.Internal as I
import Data.Vector.Dense.IO
import Data.Vector.Dense.ST

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
    
freezeVector :: (ReadVector x m, WriteVector y m, UnsafeFreezeVector y, BLAS1 e) =>
    x n e -> m (Vector n e)
freezeVector x = do
    x' <- newCopyVector x
    return (unsafeFreezeVector x')

thawVector :: (WriteVector y m, BLAS1 e) =>
    Vector n e -> m (y n e)
thawVector = newCopyVector
