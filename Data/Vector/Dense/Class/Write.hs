{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies,
        FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Write
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Write (
    WriteVector(..),
    
    -- * Creating new vectors
    newVector,
    newListVector,
    unsafeNewVector,
    
    -- * Special vectors
    newZeroVector,
    newConstantVector,
    newBasisVector,
    setBasis,
        
    module Data.Vector.Dense.Class.Read
    ) where

import Foreign

import BLAS.Tensor
import Data.Vector.Dense.Class.Read


class (WriteTensor x Int e m, ReadVector x e m) => WriteVector x e m | x -> m, m -> x where
    -- | Creates a new vector of the given length.  The elements will be 
    -- uninitialized.
    newVector_ :: Int -> m (x n e)

---------------------------  Creating Vectors --------------------------------

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
    x <- newZero n
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

---------------------------  Special Vectors --------------------------------

-- | @newBasisVector n i@ creates a vector of length @n@ that is all zero 
-- except for at position @i@, where it equal to one.
newBasisVector :: (WriteVector x e m) => Int -> Int -> m (x n e)
newBasisVector n i = do
    x <- newVector_ n
    setBasis i x
    return x


-- | @setBasis x i@ sets the @i@th coordinate of @x@ to @1@, and all other
-- coordinates to @0@.  If the vector has been scaled, it is possible that
-- @readVector x i@ will not return exactly @1@.  See 'setElem'.
setBasis :: (WriteVector x e m) => Int -> x n e -> m ()
setBasis i x
    | i < 0 || i >= dim x =
        fail $ "tried to set a vector of dimension `" ++ show (dim x) ++ "'"
               ++ " to basis vector `" ++ show i ++ "'"
    | otherwise = do
        setZero x
        unsafeWriteElem x i 1 


-- | Create a zero vector of the specified length.
newZeroVector :: (WriteVector x e m) => Int -> m (x n e)
newZeroVector = newZero
    
-- | Create a constant vector of the specified length.
newConstantVector :: (WriteVector x e m) => Int -> e -> m (x n e)
newConstantVector = newConstant
