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


class (MatrixShaped a e, HasVectorView a, Elem e) => BaseBanded_ a e where
    bandedViewArray :: ForeignPtr e -> Ptr e -> Int -> Int -> Int -> Int -> Int -> Bool -> a mn e
    arrayFromBanded :: a mn e -> (ForeignPtr e, Ptr e, Int, Int, Int, Int, Int, Bool)

class (BaseBanded_ a e, BaseVector (VectorView a) e) => BaseBanded a e

class ( BaseBanded a e, BLAS2 e, UnsafeIOToM m, ReadTensor a (Int,Int) e m
      , MMatrix a e m, MMatrix (Herm a) e m, MMatrix (Tri a) e m
      , MSolve (Tri a) e m
      , ReadVector (VectorView a) e m) => 
    ReadBanded a e m where

class (ReadBanded a e m, WriteTensor a (Int,Int) e m,
           WriteVector (VectorView a) e m) => 
    WriteBanded a e m | m -> a where


------------------------- Basic Banded Properties ---------------------------



-------------------------- ReadTensor functions -----------------------------


------------------------- WriteTensor functions -----------------------------

------------------------------ Vector views ---------------------------------


hbmv :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) => 
    e -> Herm a (k,k) e -> x k e -> e -> y k e -> m ()
hbmv alpha h (x :: x k e) beta (y :: y k e)
    | numRows h == 0 =
        return ()
    | isConj y = do
        doConjVector y
        hbmv alpha h x beta (conj y)
        doConjVector y
    | isConj x = do
        (x' :: y k e) <- newCopyVector' x
        hbmv alpha h x' beta y
    | otherwise =
        let (u,a) = hermToBase h
            n     = numCols a
            k     = case u of 
                        Upper -> numUpper a
                        Lower -> numLower a      
            u'    = case (isHermBanded a) of
                        True  -> flipUpLo u
                        False -> u
            uploA = u'
            ldA   = ldaOfBanded a
            incX  = stride x
            incY  = stride y
            withPtrA 
                  = case u' of Upper -> withBandedPtr a
                               Lower -> withBandedElemPtr a (0,0)
            x'    = unsafeVectorToIOVector x
            y'    = unsafeVectorToIOVector y
        in unsafeIOToM $
               withPtrA $ \pA ->
               withIOVector x' $ \pX ->
               withIOVector y' $ \pY -> do
                   BLAS.hbmv uploA n k alpha pA ldA pX incX beta pY incY

hbmm :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> Herm a (k,k) e -> b (k,l) e -> e -> c (k,l) e -> m ()
hbmm alpha h b beta c =
    zipWithM_ (\x y -> hbmv alpha h x beta y) (colViews b) (colViews c)

hbmv' :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) => 
    e -> Herm a (r,s) e -> x s e -> e -> y r e -> m ()
hbmv' alpha a x beta y = 
    hbmv alpha (coerceHerm a) x beta (coerceVector y)

hbmm' :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> Herm a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
hbmm' alpha a b beta c = 
    hbmm alpha (coerceHerm a) b beta (coerceMatrix c)

tbmv :: (ReadBanded a e m, WriteVector y e m) => 
    e -> Tri a (k,k) e -> y n e -> m ()
tbmv alpha t x | isConj x = do
    doConjVector x
    tbmv alpha t (conj x)
    doConjVector x

tbmv alpha t x =
    let (u,d,a) = triToBase t
        (transA,u') 
                  = if isHermBanded a then (ConjTrans, flipUpLo u) 
                                      else (NoTrans  , u)
        uploA     = u'
        diagA     = d
        n         = numCols a
        k         = case u of Upper -> numUpper a 
                              Lower -> numLower a
        ldA       = ldaOfBanded a
        incX      = stride x
        withPtrA  = case u' of 
                        Upper -> withBandedPtr a
                        Lower -> withBandedElemPtr a (0,0)
    in do
        scaleByVector alpha x
        unsafeIOToM $
            withPtrA $ \pA ->
            withVectorPtrIO x $ \pX -> do
                BLAS.tbmv uploA transA diagA n k pA ldA pX incX
  where withVectorPtrIO = withIOVector . unsafeVectorToIOVector

tbmm :: (ReadBanded a e m, WriteMatrix b e m) =>
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
tbmm 1     t b = mapM_ (\x -> tbmv 1 t x) (colViews b)
tbmm alpha t b = scaleByMatrix alpha b >> tbmm 1 t b

tbmv' :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) => 
    e -> Tri a (r,s) e -> x s e -> e -> y r e -> m ()
tbmv' alpha a (x :: x s e) beta (y  :: y r e)
    | beta /= 0 = do
        (x' :: y s e) <- newCopyVector x
        tbmv alpha (coerceTri a) x'
        scaleByVector beta y
        unsafeAxpyVector 1 x' (coerceVector y)
    | otherwise = do
        unsafeCopyVector (coerceVector y) x
        tbmv alpha (coerceTri a) (coerceVector y)

tbmm' :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> Tri a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
tbmm' alpha a (b :: b (s,t) e) beta (c :: c (r,t) e)
    | beta /= 0 = do
        (b' :: c (s,t) e) <- newCopyMatrix b
        tbmm alpha (coerceTri a) b'
        scaleByMatrix beta c
        unsafeAxpyMatrix 1 b' (coerceMatrix c)
    | otherwise = do
        unsafeCopyMatrix (coerceMatrix c) b
        tbmm alpha (coerceTri a) (coerceMatrix c)

tbsv :: (ReadBanded a e m, WriteVector y e m) => 
    e -> Tri a (k,k) e -> y n e -> m ()
tbsv alpha t x | isConj x = do
    doConjVector x
    tbsv alpha t (conj x)
    doConjVector x
tbsv alpha t x = 
    let (u,d,a) = triToBase t
        (transA,u') = if isHermBanded a then (ConjTrans, flipUpLo u)
                                        else (NoTrans  , u)
        uploA     = u'
        diagA     = d
        n         = numCols a
        k         = case u of Upper -> numUpper a 
                              Lower -> numLower a        
        ldA       = ldaOfBanded a
        incX      = stride x
        withPtrA  = case u' of 
                        Upper -> withBandedPtr a
                        Lower -> withBandedElemPtr a (0,0)
    in do
        scaleByVector alpha x
        unsafeIOToM $
            withPtrA $ \pA ->
            withVectorPtrIO x $ \pX -> do
                BLAS.tbsv uploA transA diagA n k pA ldA pX incX
  where withVectorPtrIO = withIOVector . unsafeVectorToIOVector

tbsm :: (ReadBanded a e m, WriteMatrix b e m) => 
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
tbsm 1     t b = mapM_ (\x -> tbsv 1 t x) (colViews b)
tbsm alpha t b = scaleByMatrix alpha b >> tbsm 1 t b

unsafeDoSSolveTriBanded :: (ReadBanded a e m,
    ReadVector y e m, WriteVector x e m) =>
        e -> Tri a (k,l) e -> y k e -> x l e -> m ()
unsafeDoSSolveTriBanded alpha a y x = do
    unsafeCopyVector (coerceVector x) y
    tbsv alpha (coerceTri a) (coerceVector x)

unsafeDoSSolveMatTriBanded :: (ReadBanded a e m,
    ReadMatrix c e m, WriteMatrix b e m) =>
        e -> Tri a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
unsafeDoSSolveMatTriBanded alpha a c b = do
    unsafeCopyMatrix (coerceMatrix b) c
    tbsm alpha (coerceTri a) b


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

instance (BLAS3 e) => MMatrix (Herm IOBanded) e IO where
    unsafeDoSApplyAdd    = hbmv'
    unsafeDoSApplyAddMat = hbmm'

instance (BLAS3 e) => MMatrix (Tri (STBanded s)) e (ST s) where
    unsafeDoSApply_      = tbmv
    unsafeDoSApplyMat_   = tbmm
    unsafeDoSApplyAdd    = tbmv'
    unsafeDoSApplyAddMat = tbmm'

instance (BLAS3 e) => MMatrix (Tri IOBanded) e IO where
    unsafeDoSApply_      = tbmv
    unsafeDoSApplyMat_   = tbmm
    unsafeDoSApplyAdd    = tbmv'
    unsafeDoSApplyAddMat = tbmm'

instance (BLAS3 e) => MSolve (Tri IOBanded) e IO where
    unsafeDoSSolve     = unsafeDoSSolveTriBanded
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriBanded
    unsafeDoSSolve_    = tbsv
    unsafeDoSSolveMat_ = tbsm

instance (BLAS3 e) => MSolve (Tri (STBanded s)) e (ST s) where
    unsafeDoSSolve     = unsafeDoSSolveTriBanded
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriBanded
    unsafeDoSSolve_    = tbsv
    unsafeDoSSolveMat_ = tbsm
