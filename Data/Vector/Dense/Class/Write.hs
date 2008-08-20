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
    
    -- * Copying vectors
    newCopyVector,
    copyVector,
    unsafeCopyVector,
    
    -- * Special vectors
    newZeroVector,
    newConstantVector,
    newBasisVector,
    setBasis,
        
    module Data.Vector.Dense.Class.Read
    ) where

import Control.Monad( forM_ )
import Foreign

import BLAS.C( BLAS1 )
import BLAS.Internal( checkBinaryOp )
import qualified BLAS.C as BLAS

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


---------------------------  Copying Vectors --------------------------------

-- | Creats a new vector by copying another one.
newCopyVector :: (BLAS1 e, ReadVector x e m, WriteVector y e m) => 
    x n e -> m (y n e)    
newCopyVector x
    | isConj x = 
        newCopyVector (conj x) >>= return . conj
    | otherwise = do
        y <- newVector_ (dim x)
        unsafeCopyVector y x
        return y

-- | @copyVector dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyVector :: (BLAS1 e, WriteVector y e m, ReadVector x e m) =>
    y n e -> x n e -> m ()
copyVector y x = checkBinaryOp (shape x) (shape y) $ unsafeCopyVector y x
{-# INLINE copyVector #-}

-- | Same as 'copyVector' but the sizes of the arguments are not checked.
unsafeCopyVector :: (BLAS1 e, WriteVector y e m, ReadVector x e m) =>
    y n e -> x n e -> m ()
unsafeCopyVector y x
    | isConj x && isConj y =
        unsafeCopyVector (conj y) (conj x)
    | isConj x || isConj y =
        forM_ [0..(dim x - 1)] $ \i -> do
            unsafeReadElem x i >>= unsafeWriteElem y i
    | otherwise =
        call2 BLAS.copy
  where
    call2 f =
        let n    = dim x
            incX = stride x
            incY = stride y
        in unsafeIOToM $
               withVectorPtr x $ \pX ->
                   withVectorPtr y $ \pY ->
                       f n pX incX pY incY


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
