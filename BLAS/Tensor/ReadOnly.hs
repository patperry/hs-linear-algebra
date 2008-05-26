{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.ReadOnly
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.ReadOnly (
    RTensor(..),
    readElem,
    ) where

import Data.Ix
import BLAS.Tensor.Base

-- | Class for mutable read-only tensors.
class (Tensor x i e, Monad m) => RTensor x i e m where
    -- | Get the number of elements stored in the tensor.
    getSize :: x e -> m Int
    
    -- | Get a copy of the tensor.
    newCopy :: x e -> m (x e)
    
    -- | Creates a new tensor with elements all initialized to zero.
    newZero :: i -> m (x e)
    
    -- | Creates a new tensor with elements all initialized to the 
    -- given value.
    newConstant :: i -> e -> m (x e)
    
    -- | Get the value at the specified index, without doing any 
    -- range-checking.
    unsafeReadElem :: x e -> i -> m e

    -- | Returns a lazy list of the indices in the tensor.  
    -- Because of the laziness, this function should be used with care.
    getIndices :: x e -> m [i]
    getIndices x = getAssocs x >>= return . fst . unzip

    -- | Returns a lazy list of the elements in the tensor.  
    -- Because of the laziness, this function should be used with care.
    getElems :: x e -> m [e]
    getElems x = getAssocs x >>= return . snd . unzip
    
    -- | Returns a lazy list of the elements-index pairs in the tensor.  
    -- Because of the laziness, this function should be used with care.
    getAssocs :: x e -> m [(i,e)]
    getAssocs x = do
        is <- getIndices x
        es <- getElems x
        return $ zip is es


-- | Gets the value at the specified index after checking that the argument
-- is in bounds.
readElem :: (RTensor x i e m, Show i) => x e -> i -> m e
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
