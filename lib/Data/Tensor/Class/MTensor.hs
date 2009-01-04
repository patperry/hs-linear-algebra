{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Tensor.Class.MTensor
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Tensor.Class.MTensor (
    -- * Read-only tensors
    ReadTensor(..),
    readElem,
    
    -- * Modifiable tensors
    WriteTensor(..),
    writeElem,
    modifyElem,
    swapElems,
    ) where

import Data.Elem.BLAS( Elem, BLAS1 )
import Data.Elem.Conj
import Data.Ix
import Data.Tensor.Class

-- | Class for mutable read-only tensors.
class (BaseTensor x i e, Monad m) => ReadTensor x i e m | x -> i where
    -- | Get the number of elements stored in the tensor.
    getSize :: x n e -> m Int
    
    -- | Get the value at the specified index, without doing any 
    -- range-checking.
    unsafeReadElem :: (Elem e) => x n e -> i -> m e

    -- | Returns a lazy list of the indices in the tensor.  
    -- Because of the laziness, this function should be used with care.
    -- See also "getIndices'".
    getIndices :: x n e -> m [i]

    -- | Returns a list of the indices in the tensor.  See also
    -- 'getIndices'.
    getIndices' :: x n e -> m [i]

    -- | Returns a lazy list of the elements in the tensor.  
    -- Because of the laziness, this function should be used with care.
    -- See also "getElems'".    
    getElems :: (Elem e) => x n e -> m [e]
    getElems x = getAssocs x >>= return . snd . unzip

    -- | Returns a list of the elements in the tensor.  See also
    -- 'getElems'.
    getElems' :: (Elem e) => x n e -> m [e]
    getElems' x = getAssocs' x >>= return . snd . unzip
    
    -- | Returns a lazy list of the elements-index pairs in the tensor.  
    -- Because of the laziness, this function should be used with care.
    -- See also "getAssocs'".        
    getAssocs :: (Elem e) => x n e -> m [(i,e)]
    getAssocs x = do
        is <- getIndices x
        es <- getElems x
        return $ zip is es

    -- | Returns a list of the index-elements pairs in the tensor.  See also
    -- 'getAssocs'.
    getAssocs' :: (Elem e) => x n e -> m [(i,e)]
    getAssocs' x = do
        is <- getIndices' x
        es <- getElems' x
        return $ zip is es


-- | Gets the value at the specified index after checking that the argument
-- is in bounds.
readElem :: (ReadTensor x i e m, Elem e) => x n e -> i -> m e
readElem x i =
    case (inRange b i) of
        False -> 
            fail $ "tried to get element at a index `" ++ show i ++ "'"
                   ++ " in an object with shape `" ++ show s ++ "'"
        True -> 
            unsafeReadElem x i
  where
      b = bounds x
      s = shape x


-- | Class for modifiable mutable tensors.
class (ReadTensor x i e m) => WriteTensor x i e m | x -> m where
    -- | Get the maximum number of elements that can be stored in the tensor.
    getMaxSize :: x n e -> m Int
    getMaxSize = getSize
    
    -- | Sets all stored elements to zero.
    setZero :: (Elem e) => x n e -> m ()
    setZero = setConstant 0
    
    -- | Sets all stored elements to the given value.
    setConstant :: (Elem e) => e -> x n e -> m ()

    -- | True if the value at a given index can be changed
    canModifyElem :: x n e -> i -> m Bool
    
    -- | Set the value of the element at the given index, without doing any
    -- range checking.
    unsafeWriteElem :: (Elem e) => x n e -> i -> e -> m ()
    
    -- | Modify the value of the element at the given index, without doing
    -- any range checking.
    unsafeModifyElem :: (Elem e) => x n e -> i -> (e -> e) -> m ()
    unsafeModifyElem x i f = do
        e <- unsafeReadElem x i
        unsafeWriteElem x i (f e)
    
    -- | Replace each element by a function applied to it
    modifyWith :: (Elem e) => (e -> e) -> x n e -> m ()

    -- | Same as 'swapElem' but arguments are not range-checked.
    unsafeSwapElems :: (Elem e) => x n e -> i -> i -> m ()
    unsafeSwapElems x i j = do
        e <- unsafeReadElem x i
        f <- unsafeReadElem x j
        unsafeWriteElem x j e
        unsafeWriteElem x i f
    
    -- | Replace every element with its complex conjugate.
    doConj :: (BLAS1 e) => x n e -> m ()
    doConj = modifyWith conj

    -- | Scale every element in the vector by the given value.
    scaleBy :: (BLAS1 e) => e -> x n e -> m ()
    scaleBy 1 = const $ return ()
    scaleBy k = modifyWith (k*)

    -- | Add a value to every element in a vector.
    shiftBy :: (BLAS1 e) => e -> x n e -> m ()
    shiftBy 0 = const $ return ()
    shiftBy k = modifyWith (k+)


-- | Set the value of the element at the given index.
writeElem :: (WriteTensor x i e m, Elem e) => x n e -> i -> e -> m ()
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
modifyElem :: (WriteTensor x i e m, Elem e) => x n e -> i -> (e -> e) -> m ()
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
swapElems :: (WriteTensor x i e m, Elem e) => x n e -> i -> i -> m ()
swapElems x i j
    | not ((inRange (bounds x) i) && (inRange (bounds x) j)) = 
        fail $ "Tried to swap elements `" ++ show i ++ "' and `"
               ++ show j ++ "' in a tensor of shape `" ++ show (shape x) 
               ++ "'."
    | otherwise =
        unsafeSwapElems x i j
