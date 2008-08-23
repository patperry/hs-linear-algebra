{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Apply.Read
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Apply.Read (
    -- * Matrix and vector multiplication
    getApply,
    getSApply,
    
    getApplyMat,
    getSApplyMat,

    -- * In-place multiplication
    doApply,
    doSApplyAdd,
    doApply_,
    doSApply_,
    
    doApplyMat,
    doSApplyAddMat,
    doApplyMat_,
    doSApplyMat_,

    -- * The ReadApply type class
    ReadApply(..),

    -- * Unsafe operations
    unsafeGetApply,
    unsafeDoApply,
    unsafeDoApply_,

    unsafeGetApplyMat,
    unsafeDoApplyMat,
    unsafeDoApplyMat_,

    ) where

import Control.Monad.ST( ST )

import BLAS.Elem
import BLAS.Internal( checkSquare, checkMatVecMult, checkMatVecMultAdd,
    checkMatMatMult, checkMatMatMultAdd )
import BLAS.UnsafeIOToM
import BLAS.UnsafeInterleaveM

import BLAS.Matrix.Base

import Data.Vector.Dense.Class

import Data.Matrix.Dense.Internal( Matrix )
import Data.Matrix.Dense.Class.Internal hiding ( BaseMatrix )

-- | Minimal complete definition: (unsafeDoSApplyAdd, unsafeDoSApplyAddMat)
class (BaseMatrix a e, BLAS1 e, Monad m) => ReadApply a e m where
    unsafeGetSApply :: (ReadVector x e m, WriteVector y e m) =>
        e -> a (k,l) e -> x l e -> m (y k e)
    unsafeGetSApply alpha a x = do
        y <- newVector_ (numRows a)
        unsafeDoSApplyAdd alpha a x 0 y
        return y

    unsafeGetSApplyMat :: (ReadMatrix b x e m, WriteMatrix c y e m) =>
        e -> a (r,s) e -> b (s,t) e -> m (c (r,t) e)
    unsafeGetSApplyMat alpha a b = do
        c <- newMatrix_ (numRows a, numCols b)
        unsafeDoSApplyAddMat alpha a b 0 c
        return c

    unsafeDoSApplyAdd :: (ReadVector x e m, WriteVector y e m) =>
        e -> a (k,l) e -> x l e -> e -> y k e -> m ()
    unsafeDoSApplyAdd alpha a x beta y = do
        y' <- unsafeGetSApply alpha a x
        scaleBy beta y
        unsafeAxpyVector 1 y' y

    unsafeDoSApplyAddMat :: (ReadMatrix b x e m, WriteMatrix c y e m) =>
        e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
    unsafeDoSApplyAddMat alpha a b beta c = do
        c' <- unsafeGetSApplyMat alpha a b
        scaleBy beta c
        unsafeAxpyMatrix 1 c' c

    unsafeDoSApply_ :: (WriteVector y e m) =>
        e -> a (n,n) e -> y n e -> m ()
    unsafeDoSApply_ alpha a x = do
        y <- newVector_ (dim x)
        unsafeDoSApplyAdd alpha a x 0 y
        unsafeCopyVector x y

    unsafeDoSApplyMat_ :: (WriteMatrix b y e m) =>
        e -> a (k,k) e -> b (k,l) e -> m ()
    unsafeDoSApplyMat_ alpha a b = do
        c <- newMatrix_ (shape b)
        unsafeDoSApplyAddMat alpha a b 0 c
        unsafeCopyMatrix b c


-- | Scale and apply to a vector
getSApply :: (ReadApply a e m, ReadVector x e m, WriteVector y e m) =>
    e -> a (k,l) e -> x l e -> m (y k e)
getSApply k a x =
    checkMatVecMult (shape a) (dim x) $ 
        unsafeGetSApply k a x

-- | Scale and apply to a matrix
getSApplyMat :: (ReadApply a e m, ReadMatrix b x e m, WriteMatrix c y e m) =>
    e -> a (r,s) e -> b (s,t) e -> m (c (r,t) e)
getSApplyMat k a b =
    checkMatMatMult (shape a) (shape b) $
        unsafeGetSApplyMat k a b
    
-- | @y := alpha a x + beta y@    
doSApplyAdd :: (ReadApply a e m, ReadVector x e m, WriteVector y e m) =>
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
doSApplyAdd alpha a x beta y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $
        unsafeDoSApplyAdd alpha a x beta y

-- | @c := alpha a b + beta c@
doSApplyAddMat :: (ReadApply a e m, ReadMatrix b x e m, WriteMatrix c y e m) =>
    e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
doSApplyAddMat alpha a b beta c =
    checkMatMatMultAdd (shape a) (shape b) (shape c)
        unsafeDoSApplyAddMat alpha a b beta c

-- | Apply to a vector
getApply :: (ReadApply a e m, ReadVector x e m, WriteVector y e m) =>
    a (k,l) e -> x l e -> m (y k e)
getApply a x =
    checkMatVecMult (shape a) (dim x) $ do
        unsafeGetApply a x

-- | Apply to a matrix
getApplyMat :: (ReadApply a e m, ReadMatrix b x e m, WriteMatrix c y e m) =>
    a (r,s) e -> b (s,t) e -> m (c (r,t) e)
getApplyMat a b =
    checkMatMatMult (shape a) (shape b) $
        unsafeGetApplyMat a b

-- | @ x := alpha a x@        
doSApply_ :: (ReadApply a e m, WriteVector y e m) =>
    e -> a (n,n) e -> y n e -> m ()
doSApply_ alpha a x =
    checkSquare (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoSApply_ alpha a x

-- | @ b := alpha a b@
doSApplyMat_ :: (ReadApply a e m, WriteMatrix b y e m) =>
    e -> a (s,s) e -> b (s,t) e -> m ()
doSApplyMat_ alpha a b =
    checkSquare (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoSApplyMat_ alpha a b

unsafeGetApply :: (ReadApply a e m, ReadVector x e m, WriteVector y e m) =>
    a (k,l) e -> x l e -> m (y k e)
unsafeGetApply = unsafeGetSApply 1

unsafeGetApplyMat :: (ReadApply a e m, ReadMatrix b x e m, WriteMatrix c y e m) =>
    a (r,s) e -> b (s,t) e -> m (c (r,t) e)
unsafeGetApplyMat = unsafeGetSApplyMat 1

-- | Apply to a vector and store the result in another vector
doApply :: (ReadApply a e m, ReadVector x e m, WriteVector y e m) =>
    a (k,l) e -> x l e -> y k e -> m ()
doApply a x y =
    checkMatVecMultAdd (numRows a, numCols a) (dim x) (dim y) $
        unsafeDoApply a x y

-- | Apply to a matrix and store the result in another matrix
doApplyMat :: (ReadApply a e m, ReadMatrix b x e m, WriteMatrix c y e m) =>
    a (r,s) e -> b (s,t) e -> c (r,t) e -> m ()
doApplyMat a b c =
    checkMatMatMultAdd (shape a) (shape b) (shape c) $
        unsafeDoApplyMat a b c
        
unsafeDoApply :: (ReadApply a e m, ReadVector x e m, WriteVector y e m) =>
    a (k,l) e -> x l e -> y k e -> m ()
unsafeDoApply a x y = unsafeDoSApplyAdd 1 a x 0 y

unsafeDoApplyMat :: (ReadApply a e m, ReadMatrix b x e m, WriteMatrix c y e m) =>
    a (r,s) e -> b (s,t) e -> c (r,t) e -> m ()
unsafeDoApplyMat a b c = unsafeDoSApplyAddMat 1 a b 0 c

-- | @x := a x@    
doApply_ :: (ReadApply a e m, WriteVector y e m) =>
    a (n,n) e -> y n e -> m ()
doApply_ a x =
    checkSquare (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoApply_ a x

-- | @ b := a b@
doApplyMat_ :: (ReadApply a e m, WriteMatrix b y e m) =>
    a (s,s) e -> b (s,t) e -> m ()
doApplyMat_ a b =
    checkSquare (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoApplyMat_ a b
  
unsafeDoApply_ :: (ReadApply a e m, WriteVector y e m) => 
    a (n,n) e -> y n e -> m ()
unsafeDoApply_ a x =
    unsafeDoSApply_ 1 a x

unsafeDoApplyMat_ :: (ReadApply a e m, WriteMatrix b y e m) =>
    a (s,s) e -> b (s,t) e -> m ()
unsafeDoApplyMat_ a b = 
    unsafeDoSApplyMat_ 1 a b

instance (BLAS3 e) => ReadApply IOMatrix e IO where
    unsafeDoSApplyAdd    = gemv
    unsafeDoSApplyAddMat = gemm

instance (BLAS3 e) => ReadApply (STMatrix s) e (ST s) where
    unsafeDoSApplyAdd    = gemv
    unsafeDoSApplyAddMat = gemm

instance (BLAS3 e, UnsafeIOToM m, UnsafeInterleaveM m) => ReadApply Matrix e m where
    unsafeDoSApplyAdd    = gemv
    unsafeDoSApplyAddMat = gemm
