{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Immutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor.Immutable (
    ITensor(..),
    (!),
    ) where

import BLAS.Tensor.Base
import BLAS.Elem
import Data.Ix

infixl 9 !

-- | A class for immutable tensors.
class (BLAS1 e, Tensor x i e) => ITensor x i e where
    -- | Get the numer of elements stored in the tensor.
    size :: x e -> Int
    
    -- | Get a new tensor by replacing the elements at the given indices.
    (//) :: x e -> [(i,e)] -> x e

    -- | Get the value at the given index, without doing any bounds-checking.
    unsafeAt :: x e -> i -> e
    
    -- | Same as '(//)' but doesn't do any bounds-checking.
    unsafeReplace :: x e -> [(i,e)] -> x e
    
    -- | Get the indices of the elements stored in the tensor.
    indices :: x e -> [i]
    indices = fst . unzip . assocs
    
    -- | Get the elements stored in the tensor.
    elems :: x e -> [e]
    elems = snd . unzip . assocs

    -- | Get the list of @(@index@,@ element@)@ pairs stored in the tensor.
    assocs :: x e -> [(i,e)]

    -- accum :: (e -> e' -> e) -> x e -> [(i,e')] -> x e
    
    -- | Apply a function elementwise to a tensor.
    amap :: (ITensor x i e') => (e -> e') -> x e -> x e'
    
    -- | Apply a function to pairs of elements of tensors that are the 
    -- same shape.
    azipWith :: (ITensor x i f, ITensor x i g) => (e -> f -> g) -> x e -> x f -> x g
    
    -- ixmap :: i -> (i -> i) -> x e -> x e
    -- unsafeIxMap

-- | Get the value at the given index.  Range-checks the argument.
(!) :: (ITensor x i e, Show i) => x e -> i -> e
(!) x i =
    case (inRange b i) of
        False -> 
            error $ "tried to get element at a index `" ++ show i ++ "'"
                    ++ " in an object with shape `" ++ show s ++ "'"
        True -> 
            unsafeAt x i
  where
    b = bounds x
    s = shape x
