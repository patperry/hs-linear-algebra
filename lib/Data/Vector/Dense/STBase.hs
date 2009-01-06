{-# LANGUAGE Rank2Types, FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Vector.Dense.STBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.STBase (
    -- * The @STVector@ data type
    STVector,
    runSTVector,

    unsafeIOVectorToSTVector,
    unsafeSTVectorToIOVector,
    ) where

import Control.Monad.ST
import Data.Vector.Dense.Class.Internal( STVector, unsafeIOVectorToSTVector,
    unsafeSTVectorToIOVector )
import Data.Vector.Dense.Internal ( Vector(..) )

newtype STVector s n e = ST (IOVector n e)

unsafeIOVectorToSTVector :: IOVector n e -> STVector s n e
unsafeIOVectorToSTVector = ST

unsafeSTVectorToIOVector :: STVector s n e -> IOVector n e
unsafeSTVectorToIOVector (ST x) = x


runSTVector :: (forall s . ST s (STVector s n e)) -> Vector n e
runSTVector x = 
    runST $ x >>= return . V . unsafeSTVectorToIOVector

instance (Storable e) => BaseVector (STVector s) e where
    vectorViewArray f p n s c = ST $ DV f p n s c
    {-# INLINE vectorViewArray #-}    
    
    arrayFromVector (ST x)    = arrayFromVector x
    {-# INLINE arrayFromVector #-}

instance (Storable e) => Shaped (STVector s) Int e where
    bounds = boundsVector
    shape  = shapeVector
        
instance (BLAS1 e) => ReadTensor (STVector s) Int e (ST s) where
    getSize        = getSizeVector
    getAssocs      = getAssocsVector
    getIndices     = getIndicesVector
    getElems       = getElemsVector
    getAssocs'     = getAssocsVector'
    getIndices'    = getIndicesVector'
    getElems'      = getElemsVector'
    unsafeReadElem = unsafeReadElemVector

instance (BLAS1 e) => ReadVector (STVector s) e (ST s) where    

instance (BLAS1 e) => WriteTensor (STVector s) Int e (ST s) where
    setConstant     = setConstantVector
    setZero         = setZeroVector
    canModifyElem   = canModifyElemVector
    unsafeWriteElem = unsafeWriteElemVector
    modifyWith      = modifyWithVector
    doConj          = doConjVector
    scaleBy         = scaleByVector
    shiftBy         = shiftByVector

instance (BLAS1 e) => WriteVector (STVector s) e (ST s) where
