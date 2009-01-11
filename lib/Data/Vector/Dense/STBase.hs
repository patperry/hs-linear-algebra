{-# LANGUAGE Rank2Types, FlexibleInstances, MultiParamTypeClasses #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Vector.Dense.STBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.STBase
    where

import Control.Monad
import Control.Monad.ST

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Vector.Dense.Base
import Data.Vector.Dense.IOBase

-- | Dense vectors in the 'ST' monad.  The type arguments are as follows:
--
--     * @s@: the state variable argument for the 'ST' type
--
--     * @n@: a phantom type for the dimension of the vector
--
--     * @e@: the element type of the vector.  Only certain element types
--       are supported.
--
newtype STVector s n e = STVector (IOVector n e)

-- | A safe way to create and work with a mutable vector before returning 
-- an immutable vector for later perusal. This function avoids copying
-- the vector before returning it - it uses unsafeFreezeVector internally,
-- but this wrapper is a safe interface to that function. 
runSTVector :: (forall s . ST s (STVector s n e)) -> Vector n e
runSTVector mx = 
    runST $ mx >>= \(STVector x) -> return (Vector x)

instance Shaped (STVector s) Int where
    shape (STVector x) = shapeIOVector x
    {-# INLINE shape #-}
    bounds (STVector x) = boundsIOVector x
    {-# INLINE bounds #-}

instance ReadTensor (STVector s) Int (ST s) where
    getSize (STVector x) = unsafeIOToST $ getSizeIOVector x
    {-# INLINE getSize #-}
    unsafeReadElem (STVector x) i = unsafeIOToST $ unsafeReadElemIOVector x i
    {-# INLINE unsafeReadElem #-}
    getIndices (STVector x) = unsafeIOToST $ getIndicesIOVector x
    {-# INLINE getIndices #-}
    getIndices' (STVector x) = unsafeIOToST $ getIndicesIOVector' x
    {-# INLINE getIndices' #-}
    getElems (STVector x) = unsafeIOToST $ getElemsIOVector x
    {-# INLINE getElems #-}
    getElems' (STVector x) = unsafeIOToST $ getElemsIOVector' x
    {-# INLINE getElems' #-}
    getAssocs (STVector x) = unsafeIOToST $ getAssocsIOVector x
    {-# INLINE getAssocs #-}
    getAssocs' (STVector x) = unsafeIOToST $ getAssocsIOVector' x
    {-# INLINE getAssocs' #-}

instance WriteTensor (STVector s) Int (ST s) where
    getMaxSize (STVector x) = unsafeIOToST $ getMaxSizeIOVector x
    {-# INLINE getMaxSize #-}
    setZero (STVector x) = unsafeIOToST $ setZeroIOVector x
    {-# INLINE setZero #-}
    setConstant e (STVector x) = unsafeIOToST $ setConstantIOVector e x
    {-# INLINE setConstant #-}
    canModifyElem (STVector x) i = unsafeIOToST $ canModifyElemIOVector x i
    {-# INLINE canModifyElem #-}
    unsafeWriteElem (STVector x) i e= unsafeIOToST $ unsafeWriteElemIOVector x i e
    {-# INLINE unsafeWriteElem #-}
    unsafeModifyElem (STVector x) i f = unsafeIOToST $ unsafeModifyElemIOVector x i f
    {-# INLINE unsafeModifyElem #-}
    modifyWith f (STVector x) = unsafeIOToST $ modifyWithIOVector f x
    {-# INLINE modifyWith #-}
    doConj (STVector x) = unsafeIOToST $ doConjIOVector x
    {-# INLINE doConj #-}
    scaleBy k (STVector x) = unsafeIOToST $ scaleByIOVector k x
    {-# INLINE scaleBy #-}
    shiftBy k (STVector x) = unsafeIOToST $ shiftByIOVector k x
    {-# INLINE shiftBy #-}

instance BaseVector (STVector s) where
    dim (STVector x) = dimIOVector x
    {-# INLINE dim #-}
    stride (STVector x) = strideIOVector x
    {-# INLINE stride #-}
    conjEnum (STVector x) = conjEnumIOVector x
    {-# INLINE conjEnum #-}
    conj (STVector x) = STVector (conjIOVector x)
    {-# INLINE conj #-}
    unsafeSubvectorViewWithStride s (STVector x) o n = 
        STVector (unsafeSubvectorViewWithStrideIOVector s x o n)
    {-# INLINE unsafeSubvectorViewWithStride #-}    
    unsafeVectorToIOVector (STVector x) = x
    {-# INLINE unsafeVectorToIOVector #-}
    unsafeIOVectorToVector = STVector
    {-# INLINE unsafeIOVectorToVector #-}

instance ReadVector (STVector s) (ST s) where
    unsafePerformIOWithVector (STVector x) f = unsafeIOToST $ f x
    {-# INLINE unsafePerformIOWithVector #-}
    freezeVector (STVector x) = unsafeIOToST $ freezeIOVector x
    {-# INLINE freezeVector #-}
    unsafeFreezeVector (STVector x) = unsafeIOToST $ unsafeFreezeIOVector x
    {-# INLINE unsafeFreezeVector #-}

instance WriteVector (STVector s) (ST s) where
    newVector_ = liftM STVector . unsafeIOToST . newIOVector_
    {-# INLINE newVector_ #-}
    unsafeConvertIOVector = unsafeIOToST . liftM STVector
    {-# NOINLINE unsafeConvertIOVector #-}
    thawVector = liftM STVector . unsafeIOToST . thawIOVector
    {-# INLINE thawVector #-}
    unsafeThawVector = liftM STVector . unsafeIOToST . unsafeThawIOVector
    {-# INLINE unsafeThawVector #-}
