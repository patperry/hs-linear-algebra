{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances, TypeFamilies,
        Rank2Types #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.STBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.STBase
    where

import Control.Monad
import Control.Monad.ST

import Data.Elem.BLAS( Elem, BLAS1, BLAS3 )

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Matrix.Class
import Data.Matrix.Dense.Base
import Data.Matrix.Dense.IOBase

import Data.Matrix.Herm
import Data.Matrix.TriBase

import Data.Vector.Dense.STBase

-- | Dense matrix in the 'ST' monad.  The type arguments are as follows:
--
--     * @s@: the state variable argument for the 'ST' type
--
--     * @np@: a phantom type for the shape of the matrix.  Most functions
--       will demand that this be specified as a pair.  When writing a function
--       signature, you should always prefer @STMatrix s (n,p) e@ to
--       @STMatrix s np e@.
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
newtype STMatrix s np e = STMatrix (IOMatrix np e)

-- | A safe way to create and work with a mutable matrix before returning 
-- an immutable matrix for later perusal. This function avoids copying
-- the matrix before returning it - it uses unsafeFreezeMatrix internally,
-- but this wrapper is a safe interface to that function. 
runSTMatrix :: (forall s . ST s (STMatrix s n e)) -> Matrix n e
runSTMatrix mx = 
    runST $ mx >>= \(STMatrix x) -> return (Matrix x)

instance HasVectorView (STMatrix s) where
    type VectorView (STMatrix s) = STVector s

instance Shaped (STMatrix s) (Int,Int) where
    shape (STMatrix a) = shapeIOMatrix a
    {-# INLINE shape #-}
    bounds (STMatrix a) = boundsIOMatrix a
    {-# INLINE bounds #-}

instance (Elem e) => ReadTensor (STMatrix s) (Int,Int) e (ST s) where
    getSize (STMatrix a) = unsafeIOToST $ getSizeIOMatrix a
    {-# INLINE getSize #-}
    unsafeReadElem (STMatrix a) i = unsafeIOToST $ unsafeReadElemIOMatrix a i
    {-# INLINE unsafeReadElem #-}
    getIndices (STMatrix a) = unsafeIOToST $ getIndicesIOMatrix a
    {-# INLINE getIndices #-}
    getIndices' (STMatrix a) = unsafeIOToST $ getIndicesIOMatrix' a
    {-# INLINE getIndices' #-}
    getElems (STMatrix a) = unsafeIOToST $ getElemsIOMatrix a
    {-# INLINE getElems #-}
    getElems' (STMatrix a) = unsafeIOToST $ getElemsIOMatrix' a
    {-# INLINE getElems' #-}
    getAssocs (STMatrix a) = unsafeIOToST $ getAssocsIOMatrix a
    {-# INLINE getAssocs #-}
    getAssocs' (STMatrix a) = unsafeIOToST $ getAssocsIOMatrix' a
    {-# INLINE getAssocs' #-}

instance (BLAS1 e) => WriteTensor (STMatrix s) (Int,Int) e (ST s) where
    getMaxSize (STMatrix a) =  unsafeIOToST $ getMaxSizeIOMatrix a
    {-# INLINE getMaxSize #-}
    setZero (STMatrix a) = unsafeIOToST $ setZeroIOMatrix a
    {-# INLINE setZero #-}
    setConstant e (STMatrix a) = unsafeIOToST $ setConstantIOMatrix e a
    {-# INLINE setConstant #-}
    canModifyElem (STMatrix a) i = unsafeIOToST $ canModifyElemIOMatrix a i
    {-# INLINE canModifyElem #-}
    unsafeWriteElem (STMatrix a) i e = unsafeIOToST $ unsafeWriteElemIOMatrix a i e
    {-# INLINE unsafeWriteElem #-}
    unsafeModifyElem (STMatrix a) i f = unsafeIOToST $ unsafeModifyElemIOMatrix a i f
    {-# INLINE unsafeModifyElem #-}
    modifyWith f (STMatrix a) = unsafeIOToST $ modifyWithIOMatrix f a
    {-# INLINE modifyWith #-}
    doConj (STMatrix a) = unsafeIOToST $ doConjIOMatrix a
    {-# INLINE doConj #-}
    scaleBy k (STMatrix a) = unsafeIOToST $ scaleByIOMatrix k a
    {-# INLINE scaleBy #-}
    shiftBy k (STMatrix a) = unsafeIOToST $ shiftByIOMatrix k a
    {-# INLINE shiftBy #-}

instance MatrixShaped (STMatrix s) where
    herm (STMatrix a) = STMatrix (herm a)
    {-# INLINE herm #-}
    
instance (BLAS3 e) => MMatrix (STMatrix s) e (ST s) where
    unsafeDoSApplyAdd = gemv
    {-# INLINE unsafeDoSApplyAdd #-}
    unsafeDoSApplyAddMat = gemm
    {-# INLINE unsafeDoSApplyAddMat #-}
    unsafeGetRow = unsafeGetRowMatrix
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColMatrix
    {-# INLINE unsafeGetCol #-}
    getRows = getRowsST
    {-# INLINE getRows #-}
    getCols = getColsST
    {-# INLINE getCols #-}

instance (BLAS3 e) => MMatrix (Herm (STMatrix s)) e (ST s) where
    unsafeDoSApplyAdd = hemv'
    {-# INLINE unsafeDoSApplyAdd #-}
    unsafeDoSApplyAddMat = hemm'
    {-# INLINE unsafeDoSApplyAddMat #-}    
    getRows = getRowsST
    {-# INLINE getRows #-}
    getCols = getColsST
    {-# INLINE getCols #-}

instance (BLAS3 e) => MMatrix (Tri (STMatrix s)) e (ST s) where
    unsafeDoSApplyAdd = unsafeDoSApplyAddTriMatrix
    {-# INLINE unsafeDoSApplyAdd #-}
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatTriMatrix
    {-# INLINE unsafeDoSApplyAddMat #-}
    unsafeDoSApply_ = trmv
    {-# INLINE unsafeDoSApply_ #-}
    unsafeDoSApplyMat_ = trmm
    {-# INLINE unsafeDoSApplyMat_ #-}
    getRows = getRowsST
    {-# INLINE getRows #-}
    getCols = getColsST
    {-# INLINE getCols #-}

instance (BLAS3 e) => MSolve  (Tri (STMatrix s)) e (ST s) where
    unsafeDoSSolve = unsafeDoSSolveTriMatrix
    {-# INLINE unsafeDoSSolve #-}
    unsafeDoSSolveMat = unsafeDoSSolveMatTriMatrix
    {-# INLINE unsafeDoSSolveMat #-}    
    unsafeDoSSolve_ = trsv
    {-# INLINE unsafeDoSSolve_ #-}
    unsafeDoSSolveMat_ = trsm
    {-# INLINE unsafeDoSSolveMat_ #-}

instance BaseMatrix (STMatrix s) where
    ldaMatrix (STMatrix a) = ldaIOMatrix a
    {-# INLINE ldaMatrix #-}
    transEnumMatrix (STMatrix a) = transEnumIOMatrix a
    {-# INLINE transEnumMatrix #-}
    unsafeSubmatrixView (STMatrix a) ij mn =
        STMatrix (unsafeSubmatrixViewIOMatrix a ij mn)
    {-# INLINE unsafeSubmatrixView #-}
    unsafeDiagView (STMatrix a) i = STVector (unsafeDiagViewIOMatrix a i)
    {-# INLINE unsafeDiagView #-}
    unsafeRowView (STMatrix a) i = STVector (unsafeRowViewIOMatrix a i)
    {-# INLINE unsafeRowView #-}
    unsafeColView (STMatrix a) i = STVector (unsafeColViewIOMatrix a i)
    {-# INLINE unsafeColView #-}
    maybeViewMatrixAsVector (STMatrix a) = liftM STVector (maybeViewMatrixAsVector a)
    {-# INLINE maybeViewMatrixAsVector #-}
    maybeViewVectorAsMatrix mn (STVector x) = 
        liftM STMatrix $ maybeViewVectorAsIOMatrix mn x
    {-# INLINE maybeViewVectorAsMatrix #-}    
    maybeViewVectorAsRow (STVector x) = liftM STMatrix (maybeViewVectorAsRow x)
    {-# INLINE maybeViewVectorAsRow #-}
    maybeViewVectorAsCol (STVector x) = liftM STMatrix (maybeViewVectorAsCol x)
    {-# INLINE maybeViewVectorAsCol #-}
    unsafeIOMatrixToMatrix = STMatrix
    {-# INLINE unsafeIOMatrixToMatrix #-}
    unsafeMatrixToIOMatrix (STMatrix a) = a
    {-# INLINE unsafeMatrixToIOMatrix #-}

instance (BLAS3 e) => ReadMatrix (STMatrix s) e (ST s) where
    unsafePerformIOWithMatrix (STMatrix a) f = unsafeIOToST $ f a
    {-# INLINE unsafePerformIOWithMatrix #-}
    freezeMatrix (STMatrix a) = unsafeIOToST $ freezeIOMatrix a
    {-# INLINE freezeMatrix #-}
    unsafeFreezeMatrix (STMatrix a) = unsafeIOToST $ unsafeFreezeIOMatrix a
    {-# INLINE unsafeFreezeMatrix #-}

instance (BLAS3 e) => WriteMatrix (STMatrix s) e (ST s) where
    newMatrix_ = unsafeIOToST . liftM STMatrix . newIOMatrix_
    {-# INLINE newMatrix_ #-}
    unsafeConvertIOMatrix = unsafeIOToST . liftM STMatrix
    {-# INLINE unsafeConvertIOMatrix #-}
    thawMatrix = unsafeIOToST . liftM STMatrix . thawIOMatrix
    {-# INLINE thawMatrix #-}
    unsafeThawMatrix = unsafeIOToST . liftM STMatrix . unsafeThawIOMatrix
    {-# INLINE unsafeThawMatrix #-}
    