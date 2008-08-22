-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Internal.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Internal.Base (
    BaseVector(..),
    fptrOfVector,
    offsetOfVector,
    dim,
    stride,
    isConj,
    conjVector
    ) where

import Foreign
import BLAS.Tensor.Base

class (BaseTensor x Int e) => BaseVector x e where
    
    -- | Given a pointer, an offset, a length, a stride, and a conjugacy
    -- flag, create a vector view of the underlying memory.
    vectorViewArray :: ForeignPtr e -> Int -> Int -> Int -> Bool -> x n e

    arrayFromVector :: x n e -> (ForeignPtr e, Int, Int, Int, Bool)

conjVector :: (BaseVector x e) => x n e -> x n e
conjVector x = let (f,o,n,s,c) = arrayFromVector x
               in vectorViewArray f o n s (not c)
{-# INLINE conjVector #-} 

fptrOfVector :: (BaseVector x e) => x n e -> ForeignPtr e
fptrOfVector x = let (f,_,_,_,_) = arrayFromVector x in f
{-# INLINE fptrOfVector #-}

offsetOfVector :: (BaseVector x e) => x n e -> Int
offsetOfVector x = let (_,o,_,_,_) = arrayFromVector x in o
{-# INLINE offsetOfVector #-}

dim :: (BaseVector x e) => x n e -> Int
dim x = let (_,_,n,_,_) = arrayFromVector x in n
{-# INLINE dim #-}

stride :: (BaseVector x e) => x n e -> Int
stride x = let (_,_,_,s,_) = arrayFromVector x in s
{-# INLINE stride #-}

isConj :: (BaseVector x e) => x n e -> Bool
isConj x = let (_,_,_,_,c) = arrayFromVector x in c
{-# INLINE isConj #-}
