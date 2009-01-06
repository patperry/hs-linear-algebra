{-# LANGUAGE Rank2Types, FlexibleInstances, MultiParamTypeClasses #-}
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

import Data.Elem.BLAS( Elem, BLAS1 )

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Vector.Dense.Base
import Data.Vector.Dense.IOBase


newtype STVector s n e = STVector (IOVector n e)

runSTVector :: (forall s . ST s (STVector s n e)) -> Vector n e
runSTVector mx = 
    runST $ mx >>= \(STVector x) -> return (Vector x)

instance Shaped (STVector s) Int e where
    shape (STVector x) = shapeIOVector x
    {-# INLINE shape #-}
    bounds (STVector x) = boundsIOVector x
    {-# INLINE bounds #-}

instance (Elem e) => ReadTensor (STVector s) Int e (ST s) where
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

instance (BLAS1 e) => WriteTensor (STVector s) Int e (ST s) where
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

instance (Elem e) => BaseVector (STVector s) e where
    dim (STVector x) = dimIOVector x
    {-# INLINE dim #-}
    stride (STVector x) = strideIOVector x
    {-# INLINE stride #-}
    isConj (STVector x) = isConjIOVector x
    {-# INLINE isConj #-}
    conj (STVector x) = STVector (conjIOVector x)
    {-# INLINE conj #-}
    unsafeSubvectorViewWithStride s (STVector x) o n = 
        STVector (unsafeSubvectorViewWithStrideIOVector s x o n)
    {-# INLINE unsafeSubvectorViewWithStride #-}    
    unsafeIOVectorToVector = STVector
    {-# INLINE unsafeIOVectorToVector #-}
    unsafeVectorToIOVector (STVector x) = x
    {-# INLINE unsafeVectorToIOVector #-}
    withVectorPtrIO (STVector x) = withIOVectorPtr x
    {-# INLINE withVectorPtrIO #-}

instance (BLAS1 e) => ReadVector (STVector s) e (ST s) where
    newVector_ = liftM STVector . unsafeIOToST . newIOVector_
    {-# INLINE newVector_ #-}
    withVectorPtr (STVector x) f = unsafeIOToST $ withIOVectorPtr x f
    {-# INLINE withVectorPtr #-}
    freezeVector (STVector x) = unsafeIOToST $ freezeIOVector x
    {-# INLINE freezeVector #-}
    unsafeFreezeVector (STVector x) = unsafeIOToST $ unsafeFreezeIOVector x
    {-# INLINE unsafeFreezeVector #-}
    thawVector = liftM STVector . unsafeIOToST . thawIOVector
    {-# INLINE thawVector #-}
    unsafeThawVector = liftM STVector . unsafeIOToST . unsafeThawIOVector
    {-# INLINE unsafeThawVector #-}

instance (BLAS1 e) => WriteVector (STVector s) e (ST s) where
