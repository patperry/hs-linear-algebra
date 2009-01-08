{-# LANGUAGE FlexibleInstances, FlexibleContexts, MultiParamTypeClasses
    , FunctionalDependencies, TypeFamilies, ScopedTypeVariables #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Class.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Class.Internal (
    -- * Banded matrix types
    IOBanded,
    STBanded,
    unsafeIOBandedToSTBanded,
    unsafeSTBandedToIOBanded,

    -- * Banded type classes
    BaseBanded_(..),
    BaseBanded,
    ReadBanded,
    WriteBanded,

    -- * Low-level Banded properties
    bandedViewMatrix,
    matrixFromBanded,
    ldaOfBanded,
    isHermBanded,
    hermBanded,

    -- * Bandwidth properties
    bandwidth,
    numLower,
    numUpper,

    -- * Coercing the Banded shape
    coerceBanded,

    -- * WriteTensor functions
    newBanded_,
    newZeroBanded,
    setZeroBanded,
    newConstantBanded,
    setConstantBanded,
    modifyWithBanded,
    canModifyElemBanded,
    unsafeWriteElemBanded,

    -- * Vector views
    unsafeRowViewBanded,
    unsafeColViewBanded,
    unsafeGetRowBanded,
    unsafeGetColBanded,
    
    -- * Utility functions
    shapeBanded,
    boundsBanded,
    withBandedPtr,
    withBandedElemPtr,
    indexOfBanded,
    indicesBanded,
    gbmv,
    gbmm,
    hbmv',
    hbmm',
    tbmv,
    tbmm,
    tbmv',
    tbmm',
    unsafeDoSSolveTriBanded,
    unsafeDoSSolveMatTriBanded,
    tbsv,
    tbsm,
    
    ) where

import Control.Monad
import Control.Monad.ST
import Data.Ix
import Data.List( foldl' )
import Foreign
import Unsafe.Coerce

import Data.Elem.BLAS( Elem, BLAS2, BLAS3, conjugate )
import qualified Data.Elem.BLAS as BLAS
import BLAS.Internal( diagLen )
import BLAS.Types
import BLAS.UnsafeIOToM

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Matrix.Class

import Data.Vector.Dense.Base
import Data.Vector.Dense.IOBase
import Data.Vector.Dense.STBase

import Data.Matrix.Herm
import Data.Matrix.Tri.Internal

import Data.Matrix.Dense.IOBase
import Data.Matrix.Dense.Base





------------------------- Basic Banded Properties ---------------------------



-------------------------- ReadTensor functions -----------------------------


------------------------- WriteTensor functions -----------------------------

------------------------------ Vector views ---------------------------------




--------------------------- Utility functions -------------------------------

------------------------------------ Instances ------------------------------

newtype STBanded s mn e = ST (IOBanded mn e)

instance HasVectorView (STBanded s) where
    type VectorView (STBanded s) = (STVector s)

instance (Elem e) => BaseBanded_ IOBanded e where
    bandedViewArray f p m n kl ku ld h      = BM f p m n kl ku ld h
    arrayFromBanded (BM f p m n kl ku ld h) = (f,p,m,n,kl,ku,ld,h)

instance (Elem e) => BaseBanded_ (STBanded s) e where
    bandedViewArray f p m n kl ku ld h           = ST (BM f p m n kl ku ld h)
    arrayFromBanded (ST (BM f p m n kl ku ld h)) = (f,p,m,n,kl,ku,ld,h)

instance (Elem e) => BaseBanded IOBanded e
instance (Elem e) => BaseBanded (STBanded s) e

instance (Elem e) => Shaped (STBanded s) (Int,Int) e where
    shape  = shapeBanded
    bounds = boundsBanded

    
instance (Elem e) => MatrixShaped (STBanded s) e where
    herm = hermBanded

instance (BLAS3 e) => ReadBanded IOBanded     e IO
instance (BLAS3 e) => ReadBanded (STBanded s) e (ST s)

instance (BLAS3 e) => ReadTensor (STBanded s) (Int,Int) e (ST s) where
    getSize        = getSizeBanded
    getAssocs      = getAssocsBanded
    getIndices     = getIndicesBanded
    getElems       = getElemsBanded
    getAssocs'     = getAssocsBanded'
    getIndices'    = getIndicesBanded'
    getElems'      = getElemsBanded'
    unsafeReadElem = unsafeReadElemBanded

instance (BLAS3 e) => WriteBanded IOBanded     e IO where
instance (BLAS3 e) => WriteBanded (STBanded s) e (ST s) where

instance (BLAS3 e) => WriteTensor (STBanded s) (Int,Int) e (ST s) where
    setConstant     = setConstantBanded
    setZero         = setZeroBanded
    modifyWith      = modifyWithBanded
    unsafeWriteElem = unsafeWriteElemBanded
    canModifyElem   = canModifyElemBanded

instance (BLAS3 e) => MMatrix (STBanded s) e (ST s) where
    unsafeDoSApplyAdd    = gbmv
    unsafeDoSApplyAddMat = gbmm
    unsafeGetRow         = unsafeGetRowBanded
    unsafeGetCol         = unsafeGetColBanded

instance (BLAS3 e) => MMatrix (Herm (STBanded s)) e (ST s) where
    unsafeDoSApplyAdd    = hbmv'
    unsafeDoSApplyAddMat = hbmm'

instance (BLAS3 e) => MMatrix (Tri (STBanded s)) e (ST s) where
    unsafeDoSApply_      = tbmv
    unsafeDoSApplyMat_   = tbmm
    unsafeDoSApplyAdd    = tbmv'
    unsafeDoSApplyAddMat = tbmm'

instance (BLAS3 e) => MSolve (Tri (STBanded s)) e (ST s) where
    unsafeDoSSolve     = unsafeDoSSolveTriBanded
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriBanded
    unsafeDoSSolve_    = tbsv
    unsafeDoSSolveMat_ = tbsm
