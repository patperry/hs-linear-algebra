{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Tensor.Class.MTensor
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Overloaded interface for mutable tensors.  This modules includes tensors
-- which can be /read/ in a monad, 'ReadTensor', as well as tensors which
-- can be /modified/ in a monad, 'WriteTensor'.
--

module Data.Tensor.Class.MTensor (
    -- * Read-only tensor type class
    ReadTensor(..),
    readElem,
    
    -- * Modifiable tensor type class
    WriteTensor(..),
    writeElem,
    modifyElem,
    swapElems,
    ) where

import Data.Elem.BLAS
import Data.Ix
import Data.Tensor.Class

-- | Class for mutable read-only tensors.
class (Shaped x i, Monad m) => ReadTensor x i m where
    -- | Get the number of elements stored in the tensor.
    getSize :: x e -> m Int
    
    -- | Get the value at the specified index, without doing any 
    -- range-checking.
    unsafeReadElem :: x e -> i -> m e

    -- | Returns a lazy list of the indices in the tensor.  
    -- Because of the laziness, this function should be used with care.
    -- See also "getIndices'".
    getIndices :: x e -> m [i]

    -- | Returns a list of the indices in the tensor.  See also
    -- 'getIndices'.
    getIndices' :: x e -> m [i]

    -- | Returns a lazy list of the elements in the tensor.  
    -- Because of the laziness, this function should be used with care.
    -- See also "getElems'".    
    getElems :: x e -> m [e]
    getElems x = getAssocs x >>= return . snd . unzip
    {-# INLINE getElems #-}

    -- | Returns a list of the elements in the tensor.  See also
    -- 'getElems'.
    getElems' :: x e -> m [e]
    getElems' x = getAssocs' x >>= return . snd . unzip
    {-# INLINE getElems' #-}
    
    -- | Returns a lazy list of the elements-index pairs in the tensor.  
    -- Because of the laziness, this function should be used with care.
    -- See also "getAssocs'".        
    getAssocs :: x e -> m [(i,e)]
    getAssocs x = do
        is <- getIndices x
        es <- getElems x
        return $ zip is es
    {-# INLINE getAssocs #-}

    -- | Returns a list of the index-elements pairs in the tensor.  See also
    -- 'getAssocs'.
    getAssocs' :: x e -> m [(i,e)]
    getAssocs' x = do
        is <- getIndices' x
        es <- getElems' x
        return $ zip is es
    {-# INLINE getAssocs' #-}

-- | Gets the value at the specified index after checking that the argument
-- is in bounds.
readElem :: (ReadTensor x i m) => x e -> i -> m e
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
{-# INLINE readElem #-}

-- | Class for modifiable mutable tensors.
class (ReadTensor x i m) => WriteTensor x i m where
    -- | Get the maximum number of elements that can be stored in the tensor.
    getMaxSize :: x e -> m Int
    getMaxSize = getSize
    {-# INLINE getMaxSize #-}
    
    -- | Sets all stored elements to zero.
    setZero :: (Num e) => x e -> m ()
    setZero = setConstant 0
    {-# INLINE setZero #-}
    
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
    {-# INLINE unsafeModifyElem #-}
    
    -- | Replace each element by a function applied to it
    modifyWith :: (e -> e) -> x e -> m ()

    -- | Same as 'swapElem' but arguments are not range-checked.
    unsafeSwapElems :: x e -> i -> i -> m ()
    unsafeSwapElems x i j = do
        e <- unsafeReadElem x i
        f <- unsafeReadElem x j
        unsafeWriteElem x j e
        unsafeWriteElem x i f
    {-# INLINE unsafeSwapElems #-}
    
    -- | Set all of the elements of the tensor.  It almost equivalent to:
    --
    -- @
    --     setElems x es = do
    --          is <- getIndices x
    --          zipWithM_ (unsafeWriteElem x) is es
    -- @
    --
    -- execpt that it will call 'fail' if the input list does not have size
    -- equal to the size of the tensor.
    setElems :: x e -> [e] -> m ()
    setElems x = 
        let go _ []     []     = return ()
            go n (i:is) (e:es) = unsafeWriteElem x i e >> go (n+1) is es
            go n []     _      = 
                   fail $ "setElems <tensor of size " ++ show n ++ ">"
                        ++ " <list of length at least " ++ show (n+1) ++ ">:"
                        ++ " too many elements given."
            go n _      []     = do
                   s <- getSize x
                   fail $ "setElems <tensor of size " ++ show s ++ ">"
                        ++ " <list of length " ++ show n ++ ">:"
                        ++ " not enough elements given."
        in \es -> do
            is <- getIndices x
            go (0 :: Int) is es
    {-# INLINE setElems #-}
    
    -- | Replace every element with its complex conjugate.
    doConj :: (BLAS1 e) => x e -> m ()
    doConj = modifyWith conjugate
    {-# INLINE doConj #-}

    -- | Scale every element in the vector by the given value.
    scaleBy :: (BLAS1 e) => e -> x e -> m ()
    scaleBy 1 = const $ return ()
    scaleBy k = modifyWith (k*)
    {-# INLINE scaleBy #-}

    -- | Add a value to every element in a vector.
    shiftBy :: (BLAS1 e) => e -> x e -> m ()
    shiftBy 0 = const $ return ()
    shiftBy k = modifyWith (k+)
    {-# INLINE shiftBy #-}

-- | Set the value of the element at the given index.
writeElem :: (WriteTensor x i m) => x e -> i -> e -> m ()
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
{-# INLINE writeElem #-}

-- | Update the value of the element at the given index.
modifyElem :: (WriteTensor x i m) => x e -> i -> (e -> e) -> m ()
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
{-# INLINE modifyElem #-}

-- | Swap the values stored at two positions in the tensor.
swapElems :: (WriteTensor x i m) => x e -> i -> i -> m ()
swapElems x i j
    | not ((inRange (bounds x) i) && (inRange (bounds x) j)) = 
        fail $ "Tried to swap elements `" ++ show i ++ "' and `"
               ++ show j ++ "' in a tensor of shape `" ++ show (shape x) 
               ++ "'."
    | otherwise =
        unsafeSwapElems x i j
{-# INLINE swapElems #-}
