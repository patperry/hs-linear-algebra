{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Read
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Read (
    ReadTensor(..),
    readElem,
    ) where

import Data.Ix
import BLAS.Elem
import BLAS.Tensor.Base

-- | Class for mutable read-only tensors.
class (BaseTensor x i, Monad m) => ReadTensor x i m | x -> i where
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
readElem :: (ReadTensor x i m, Elem e) => x n e -> i -> m e
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
