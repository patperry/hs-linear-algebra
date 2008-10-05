{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Write
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Write (
    WriteTensor(..),
    writeElem,
    modifyElem,
    swapElem,
    ) where

import Data.Ix( inRange )
import BLAS.Tensor.Base
import BLAS.Tensor.Read

-- | Class for modifiable mutable tensors.
class (Num e, ReadTensor x i e m) => WriteTensor x i e m | x -> m where
    -- | Get the maximum number of elements that can be stored in the tensor.
    getMaxSize :: x n e -> m Int
    getMaxSize = getSize
    
    -- | Sets all stored elements to zero.
    setZero :: x n e -> m ()
    setZero = setConstant 0
    
    -- | Sets all stored elements to the given value.
    setConstant :: e -> x n e -> m ()

    -- | True if the value at a given index can be changed
    canModifyElem :: x n e -> i -> m Bool
    
    -- | Set the value of the element at the given index, without doing any
    -- range checking.
    unsafeWriteElem :: x n e -> i -> e -> m ()
    
    -- | Modify the value of the element at the given index, without doing
    -- any range checking.
    unsafeModifyElem :: x n e -> i -> (e -> e) -> m ()
    unsafeModifyElem x i f = do
        e <- unsafeReadElem x i
        unsafeWriteElem x i (f e)
    
    -- | Replace each element by a function applied to it
    modifyWith :: (e -> e) -> x n e -> m ()

    -- | Same as 'swapElem' but arguments are not range-checked.
    unsafeSwapElem :: x n e -> i -> i -> m ()
    unsafeSwapElem x i j = do
        e <- unsafeReadElem x i
        f <- unsafeReadElem x j
        unsafeWriteElem x j e
        unsafeWriteElem x i f
    

-- | Set the value of the element at the given index.
writeElem :: (WriteTensor x i e m) => x n e -> i -> e -> m ()
writeElem x i e = do
    ok <- canModifyElem x i
    case ok && inRange (bounds x) i of
        False -> 
            fail $ "tried to set element at index `" ++ show i ++ "'"
                   ++ " in an object with shape `" ++ show s ++ "'"
                   ++ " but that element cannot be modified"
        True ->
            unsafeWriteElem x i e
  where
    s = shape x

-- | Update the value of the element at the given index.
modifyElem :: (WriteTensor x i e m) => x n e -> i -> (e -> e) -> m ()
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

-- | Swap the values stored at two positions in the tensor.
swapElem :: (WriteTensor x i e m) => x n e -> i -> i -> m ()
swapElem x i j
    | not ((inRange (bounds x) i) && (inRange (bounds x) j)) = 
        fail $ "Tried to swap elements `" ++ show i ++ "' and `"
               ++ show j ++ "' in a tensor of shape `" ++ show (shape x) 
               ++ "'."
    | otherwise =
        unsafeSwapElem x i j
