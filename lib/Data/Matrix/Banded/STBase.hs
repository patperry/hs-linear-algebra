{-# LANGUAGE Rank2Types #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.STBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.STBase
    where

import Control.Monad
import Control.Monad.ST

import Data.Matrix.Class
import Data.Matrix.Class.MMatrixBase
import Data.Matrix.Class.MSolveBase

import Data.Matrix.Herm
import Data.Matrix.Tri

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Matrix.Dense.STBase( STMatrix(..) )
import Data.Vector.Dense.STBase( STVector(..) )

import Data.Matrix.Banded.Base
import Data.Matrix.Banded.IOBase( IOBanded )
import qualified Data.Matrix.Banded.IOBase as IO

-- | Banded matrix in the 'ST' monad.  The type arguments are as follows:
--
--     * @s@: the state variable argument for the 'ST' type
--
--     * @np@: a phantom type for the shape of the matrix.  Most functions
--       will demand that this be specified as a pair.  When writing a function
--       signature, you should always prefer @STBanded s (n,p) e@ to
--       @STBanded s np e@.
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
newtype STBanded s np e = STBanded (IOBanded np e)

-- | A safe way to create and work with a mutable banded matrix before returning 
-- an immutable one for later perusal. This function avoids copying
-- the matrix before returning it - it uses unsafeFreezeBanded internally,
-- but this wrapper is a safe interface to that function. 
runSTBanded :: (forall s . ST s (STBanded s (n,p) e)) -> Banded (n,p) e
runSTBanded mx = 
    runST $ mx >>= \(STBanded x) -> return (Banded x)

instance HasVectorView (STBanded s) where
    type VectorView (STBanded s) = STVector s

instance HasMatrixStorage (STBanded s) where
    type MatrixStorage (STBanded s) = (STMatrix s)

instance Shaped (STBanded s) (Int,Int) where
    shape (STBanded a) = IO.shapeIOBanded a
    {-# INLINE shape #-}
    bounds (STBanded a) = IO.boundsIOBanded a
    {-# INLINE bounds #-}

instance MatrixShaped (STBanded s) where
    herm (STBanded a) = STBanded $ IO.hermIOBanded a
    {-# INLINE herm #-}

instance ReadTensor (STBanded s) (Int,Int) (ST s) where
    getSize (STBanded a) = unsafeIOToST $ IO.getSizeIOBanded a
    {-# INLINE getSize #-}
    getAssocs (STBanded a) = unsafeIOToST $ IO.getAssocsIOBanded a
    {-# INLINE getAssocs #-}
    getIndices (STBanded a) = unsafeIOToST $ IO.getIndicesIOBanded a
    {-# INLINE getIndices #-}
    getElems (STBanded a) = unsafeIOToST $ IO.getElemsIOBanded a
    {-# INLINE getElems #-}
    getAssocs' (STBanded a) = unsafeIOToST $ IO.getAssocsIOBanded' a
    {-# INLINE getAssocs' #-}
    getIndices' (STBanded a) = unsafeIOToST $ IO.getIndicesIOBanded' a
    {-# INLINE getIndices' #-}
    getElems' (STBanded a) = unsafeIOToST $ IO.getElemsIOBanded' a
    {-# INLINE getElems' #-}
    unsafeReadElem (STBanded a) i = unsafeIOToST $ IO.unsafeReadElemIOBanded a i
    {-# INLINE unsafeReadElem #-}

instance WriteTensor (STBanded s) (Int,Int) (ST s) where
    setConstant k (STBanded a) = unsafeIOToST $ IO.setConstantIOBanded k a
    {-# INLINE setConstant #-}
    setZero (STBanded a) = unsafeIOToST $ IO.setZeroIOBanded a
    {-# INLINE setZero #-}
    modifyWith f (STBanded a) = unsafeIOToST $ IO.modifyWithIOBanded f a
    {-# INLINE modifyWith #-}
    unsafeWriteElem (STBanded a) i e = unsafeIOToST $ IO.unsafeWriteElemIOBanded a i e
    {-# INLINE unsafeWriteElem #-}
    canModifyElem (STBanded a) i = unsafeIOToST $ IO.canModifyElemIOBanded a i
    {-# INLINE canModifyElem #-}

instance BaseBanded (STBanded s) where
    numLower (STBanded a) = IO.numLowerIOBanded a
    {-# INLINE numLower #-}
    numUpper (STBanded a) = IO.numUpperIOBanded a
    {-# INLINE numUpper #-}
    bandwidths (STBanded a) = IO.bandwidthsIOBanded a
    {-# INLINE bandwidths #-}
    ldaBanded (STBanded a) = IO.ldaIOBanded a
    {-# INLINE ldaBanded #-}
    transEnumBanded (STBanded a) = IO.transEnumIOBanded a
    {-# INLINE transEnumBanded #-}
    maybeMatrixStorageFromBanded (STBanded a) = 
        liftM STMatrix $ IO.maybeMatrixStorageFromIOBanded a
    {-# INLINE maybeMatrixStorageFromBanded #-}
    maybeBandedFromMatrixStorage mn kl (STMatrix a) = 
        liftM STBanded $ IO.maybeIOBandedFromMatrixStorage mn kl a
    {-# INLINE maybeBandedFromMatrixStorage #-}
    viewVectorAsBanded mn (STVector x) = STBanded $ IO.viewVectorAsIOBanded mn x
    {-# INLINE viewVectorAsBanded #-}
    maybeViewBandedAsVector (STBanded a) = 
        liftM STVector $ IO.maybeViewIOBandedAsVector a
    {-# INLINE maybeViewBandedAsVector #-}        
    unsafeDiagViewBanded (STBanded a) i = 
        STVector $ IO.unsafeDiagViewIOBanded a i
    {-# INLINE unsafeDiagViewBanded #-}
    unsafeRowViewBanded (STBanded a) i = 
        case IO.unsafeRowViewIOBanded a i of (nb,x,na) -> (nb, STVector x, na)
    {-# INLINE unsafeRowViewBanded #-}
    unsafeColViewBanded (STBanded a) j = 
        case IO.unsafeColViewIOBanded a j of (nb,x,na) -> (nb, STVector x, na)
    {-# INLINE unsafeColViewBanded #-}
    unsafeIOBandedToBanded = STBanded
    {-# INLINE unsafeIOBandedToBanded #-}
    unsafeBandedToIOBanded (STBanded a) = a
    {-# INLINE unsafeBandedToIOBanded #-}
    
instance ReadBanded (STBanded s) (ST s) where
    unsafePerformIOWithBanded (STBanded a) f = unsafeIOToST $ f a
    {-# INLINE unsafePerformIOWithBanded #-}
    freezeBanded (STBanded a) = unsafeIOToST $ freezeIOBanded a
    {-# INLINE freezeBanded #-}
    unsafeFreezeBanded (STBanded a) = unsafeIOToST $ unsafeFreezeIOBanded a
    {-# INLINE unsafeFreezeBanded #-}

instance WriteBanded (STBanded s) (ST s) where
    newBanded_ mn bw = unsafeIOToST $ liftM STBanded $ IO.newIOBanded_ mn bw
    {-# INLINE newBanded_ #-}
    thawBanded = unsafeIOToST . liftM STBanded . thawIOBanded
    {-# INLINE thawBanded #-}
    unsafeThawBanded = unsafeIOToST . liftM STBanded . unsafeThawIOBanded
    {-# INLINE unsafeThawBanded #-}
    unsafeConvertIOBanded ioa = unsafeIOToST $ liftM STBanded ioa
    {-# INLINE unsafeConvertIOBanded #-}
        

instance MMatrix (STBanded s) (ST s) where
    unsafeDoSApplyAdd    = gbmv
    {-# INLINE unsafeDoSApplyAdd #-}
    unsafeDoSApplyAddMat = gbmm
    {-# INLINE unsafeDoSApplyAddMat #-}
    unsafeGetRow         = unsafeGetRowBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol         = unsafeGetColBanded
    {-# INLINE unsafeGetCol #-}
    getRows = getRowsST
    {-# INLINE getRows #-}
    getCols = getColsST
    {-# INLINE getCols #-}

instance MMatrix (Herm (STBanded s)) (ST s) where
    unsafeDoSApplyAdd    = hbmv
    {-# INLINE unsafeDoSApplyAdd #-}    
    unsafeDoSApplyAddMat = hbmm
    {-# INLINE unsafeDoSApplyAddMat #-}    
    getRows = getRowsST
    {-# INLINE getRows #-}
    getCols = getColsST
    {-# INLINE getCols #-}
    unsafeGetRow = unsafeGetRowHermBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColHermBanded
    {-# INLINE unsafeGetCol #-}

instance MMatrix (Tri (STBanded s)) (ST s) where
    unsafeDoSApply_      = tbmv
    {-# INLINE unsafeDoSApply_ #-}        
    unsafeDoSApplyMat_   = tbmm
    {-# INLINE unsafeDoSApplyMat_ #-}    
    unsafeDoSApplyAdd    = tbmv'
    {-# INLINE unsafeDoSApplyAdd #-}    
    unsafeDoSApplyAddMat = tbmm'
    {-# INLINE unsafeDoSApplyAddMat #-}    
    getRows = getRowsST
    {-# INLINE getRows #-}
    getCols = getColsST
    {-# INLINE getCols #-}
    unsafeGetRow = unsafeGetRowTriBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColTriBanded
    {-# INLINE unsafeGetCol #-}

instance MSolve (Tri (STBanded s)) (ST s) where
    unsafeDoSSolve_    = tbsv
    {-# INLINE unsafeDoSSolve_ #-}    
    unsafeDoSSolveMat_ = tbsm
    {-# INLINE unsafeDoSSolveMat_ #-}    
    unsafeDoSSolve     = tbsv'
    {-# INLINE unsafeDoSSolve #-}
    unsafeDoSSolveMat  = tbsm'
    {-# INLINE unsafeDoSSolveMat #-}

