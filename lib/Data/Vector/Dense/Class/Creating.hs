{-# OPTIONS_HADDOCK hide, prune #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Creating
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Creating (
    -- * Creating new vectors
    newVector_,
    newVector,
    newListVector,
    unsafeNewVector,

    ) where

import Foreign

import BLAS.UnsafeIOToM

import Data.Tensor.Class.Write
import Data.Vector.Dense.Class.Internal


-- | Creates a new vector with the given association list.  Unspecified
-- indices will get initialized to zero.
newVector :: (WriteVector x e m) => Int -> [(Int,e)] -> m (x n e)
newVector = newVectorHelp writeElem

-- | Same as 'newVector' but indices are not range-checked.
unsafeNewVector :: (WriteVector x e m) => Int -> [(Int,e)] -> m (x n e)
unsafeNewVector = newVectorHelp unsafeWriteElem

newVectorHelp :: (WriteVector x e m) => 
    (x n e -> Int -> e -> m ()) -> Int -> [(Int,e)] -> m (x n e)
newVectorHelp set n ies = do
    x <- newZeroVector n
    mapM_ (uncurry $ set x) ies
    return x

-- | Creates a new vector of the given dimension with the given elements.
-- If the list has length less than the passed-in dimenson, the tail of
-- the vector will be uninitialized.
newListVector :: (WriteVector x e m) => Int -> [e] -> m (x n e)
newListVector n es = do
    x <- newVector_ n
    unsafeIOToM $ withVectorPtr x $ flip pokeArray $ take n $ es ++ (repeat 0)
    return x
