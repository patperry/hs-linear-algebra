-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Mutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Mutable (
    MTensor(..),
    writeElem,
    modifyElem,
    ) where

import BLAS.Tensor.Base
import BLAS.Tensor.ReadOnly

-- | Class for modifiable mutable tensors.
class (RTensor x i e m) => (MTensor x i e m) where
    -- | Get the maximum number of elements that can be stored in the tensor.
    getMaxSize :: x e -> m Int
    getMaxSize = getSize
    
    -- | Sets all stored elements to zero.
    setZero :: x e -> m ()
    
    -- | Sets all stored elements to the given value.
    setConstant :: e -> x e -> m ()

    -- | True if the value at a given index can be changed
    canModifyElem :: x e -> i -> m Bool
    
    -- | Set the value of the element at the given index, without doing any
    -- range checking.
    unsafeWriteElem :: x e -> i -> e -> m ()
    
    -- | Modify the value of the element at the given index, without doing
    -- any range checking.
    unsafeModifyElem :: x e -> i -> (e -> e) -> m ()
    unsafeModifyElem x i f = do
        e <- unsafeReadElem x i
        unsafeWriteElem x i (f e)
    
    -- | Replace each element by a function applied to it
    modifyWith :: (e -> e) -> x e -> m ()
    
-- | Set the value of the element at the given index.
writeElem :: (MTensor x i e m, Show i) => x e -> i -> e -> m ()
writeElem x i e = do
    ok <- canModifyElem x i
    case ok of
        False -> 
            fail $ "tried to set element at index `" ++ show i ++ "'"
                   ++ " in an object with shape `" ++ show s ++ "'"
                   ++ " but that element cannot be modified"
        True ->
            unsafeWriteElem x i e
  where
    s = shape x

-- | Update the value of the element at the given index.
modifyElem :: (MTensor x i e m, Show i) => x e -> i -> (e -> e) -> m ()
modifyElem x i f = do
    ok <- canModifyElem x i
    case ok of
        False -> 
            fail $ "tried to modify element at index `" ++ show i ++ "'"
                   ++ " in an object with shape `" ++ show s ++ "'"
                   ++ " but that element cannot be modified"
        True ->
            unsafeModifyElem x i f
  where
    s = shape x
    