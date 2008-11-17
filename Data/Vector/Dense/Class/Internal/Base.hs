{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts #-}
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
    withVectorPtr,
    dim,
    stride,
    isConj,
    conjVector
    ) where

import Foreign
import BLAS.Tensor.Base

class (BaseTensor x Int e, Storable e) => BaseVector x e where
    
    -- | Give a storage region, a base pointer, a length, a stride, and a conjugacy
    -- flag, create a vector view of the underlying memory.
    vectorViewArray :: ForeignPtr e -> Ptr e -> Int -> Int -> Bool -> x n e

    arrayFromVector :: x n e -> (ForeignPtr e, Ptr e, Int, Int, Bool)

conjVector :: (BaseVector x e) => x n e -> x n e
conjVector x = let (f,o,n,s,c) = arrayFromVector x
               in vectorViewArray f o n s (not c)
{-# INLINE conjVector #-} 

withVectorPtr :: (BaseVector x e) => 
    x n e -> (Ptr e -> IO a) -> IO a
withVectorPtr x f = 
    let (fp,p,_,_,_) = arrayFromVector x
    in do
        a <- f p
        touchForeignPtr fp
        return a
{-# INLINE withVectorPtr #-}

dim :: (BaseVector x e) => x n e -> Int
dim x = let (_,_,n,_,_) = arrayFromVector x in n
{-# INLINE dim #-}

stride :: (BaseVector x e) => x n e -> Int
stride x = let (_,_,_,s,_) = arrayFromVector x in s
{-# INLINE stride #-}

isConj :: (BaseVector x e) => x n e -> Bool
isConj x = let (_,_,_,_,c) = arrayFromVector x in c
{-# INLINE isConj #-}
