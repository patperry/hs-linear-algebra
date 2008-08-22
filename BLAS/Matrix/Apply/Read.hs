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
    ReadApply(..),

{-    
    -- * Matrix and vector multiplication
    getApply,
    getApplyMat,
    getSApply,
    getSApplyMat,
    
    -- * In-place multiplication
    doApply,
    doApplyAdd,
    doApplyMat,
    doApplyAddMat,
    doSApply,
    doSApplyAdd,
    doSApplyMat,
    doSApplyAddMat,
    doApply_,
    doApplyMat_,
    doSApply_,
    doSApplyMat_,
    
    -- * Unsafe operations
    unsafeGetApply,
    unsafeGetApplyMat,
    unsafeGetSApply,
    unsafeGetSApplyMat,
-}    
    ) where

import BLAS.Elem ( BLAS1 )
--import BLAS.Internal ( checkMatVecMult, checkMatMatMult, checkMatVecMultAdd,
--    checkMatMatMultAdd, checkSquare )
import BLAS.Matrix.Base

import Data.Vector.Dense.Class
-- import Data.Vector.Dense.Operations

import Data.Matrix.Dense.Class hiding ( BaseMatrix )

-- | Minimal complete definition: (unsafeDoSApplyAdd, unsafeDoSApplyAddMat)
-- or (unsafeDoSApply, unsafeDoSApplyMat)
class (BaseMatrix a e, BLAS1 e, Monad m) => ReadApply a e m where
    unsafeDoApply :: (ReadVector x e m, WriteVector y e m) =>
        a (k,l) e -> x l e -> y k e -> m ()
    unsafeDoApply = unsafeDoSApply 1
      
    unsafeDoApplyAdd :: (ReadVector x e m, WriteVector y e m) =>
        a (k,l) e -> x l e -> e -> y k e -> m ()
    unsafeDoApplyAdd a x beta y =
        unsafeDoSApplyAdd 1 a x beta y

    unsafeDoSApply :: (ReadVector x e m, WriteVector y e m) =>
        e -> a (k,l) e -> x l e -> y k e -> m ()
    unsafeDoSApply alpha a x y =
        unsafeDoSApplyAdd alpha a x 0 y

{-
    unsafeDoApplyMat :: a (m,k) e -> DMatrix t (k,n) e -> IOMatrix (m,n) e -> m ()
    unsafeDoApplyMat a b c =
        unsafeDoSApplyMat 1 a b c
        
    unsafeDoApplyAddMat :: a (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> m ()
    unsafeDoApplyAddMat a b beta c =
        unsafeDoSApplyAddMat 1 a b beta c

    unsafeDoSApplyMat :: e -> a (m,k) e -> DMatrix t (k,n) e -> IOMatrix (m,n) e -> m ()
    unsafeDoSApplyMat alpha a b c =
        unsafeDoSApplyAddMat alpha a b 0 c
-}
    
    unsafeDoApply_ :: (WriteVector y e m) => 
        a (n,n) e -> y n e -> m ()
    unsafeDoApply_ a x =
        unsafeDoSApply_ 1 a x

{-        
    unsafeDoApplyMat_ :: a (m,m) e -> IOMatrix (m,n) e -> m ()
    unsafeDoApplyMat_ a b =
        unsafeDoSApplyMat_ 1 a b
-}
    unsafeDoSApply_ :: (WriteVector y e m) =>
        e -> a (n,n) e -> y n e -> m ()
    unsafeDoSApply_ alpha a x = do
        y <- newVector_ (dim x)
        unsafeDoSApply alpha a x y
        unsafeCopyVector x y

{-
    unsafeDoSApplyMat_ :: e -> a (m,m) e -> IOMatrix (m,n) e -> m ()
    unsafeDoSApplyMat_ alpha a b = do
        c <- newMatrix_ (shape b)
        unsafeDoSApplyMat alpha a b c
        unsafeCopyMatrix b c
-}        
    unsafeDoSApplyAdd :: (ReadVector x e m, WriteVector y e m) =>
        e -> a (k,l) e -> x l e -> e -> y k e -> m ()
    unsafeDoSApplyAdd alpha a x beta y = undefined --do
        --y' <- unsafeGetApply a x
        --scaleBy beta y
        --unsafeAxpy alpha y' y
{-        
    unsafeDoSApplyAddMat :: e -> a (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> m ()
    unsafeDoSApplyAddMat alpha a b beta c = do
        c' <- unsafeGetApplyMat a b
        M.scaleBy beta c
        M.axpy alpha c' c
-}

{-
unsafeGetApply :: (RMatrix a e) =>
    a (k,l) e -> x l e -> IO (DVector r m e)
unsafeGetApply = unsafeGetSApply 1

unsafeGetApplyMat :: (RMatrix a e) =>
    a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
unsafeGetApplyMat = unsafeGetSApplyMat 1

unsafeGetSApply :: (RMatrix a e) =>
    e -> a (k,l) e -> x l e -> IO (DVector r m e)
unsafeGetSApply alpha a x = do
    y <- newZero (numRows a)
    unsafeDoSApply alpha a x y
    return (unsafeCoerce y)

unsafeGetSApplyMat :: (RMatrix a e) =>
    e -> a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
unsafeGetSApplyMat alpha a b = do
    c <- newZero (numRows a, numCols b)
    unsafeDoSApplyMat alpha a b c
    return (unsafeCoerce c)




-- | Apply to a vector
getApply :: (RMatrix a e) =>
    a (k,l) e -> x l e -> IO (DVector r m e)
getApply a x =
    checkMatVecMult (shape a) (dim x) $ do
        unsafeGetApply a x

-- | Apply to a matrix
getApplyMat :: (RMatrix a e) =>
    a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
getApplyMat a b =
    checkMatMatMult (shape a) (shape b) $
        unsafeGetApplyMat a b

-- | Scale and apply to a vector
getSApply :: (RMatrix a e) =>
    e -> a (k,l) e -> x l e -> IO (DVector r m e)
getSApply k a x =
    checkMatVecMult (shape a) (dim x) $ 
        unsafeGetSApply k a x

-- | Scale and apply to a matrix
getSApplyMat :: (RMatrix a e) =>
    e -> a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
getSApplyMat k a b =
    checkMatMatMult (shape a) (shape b) $
        unsafeGetSApplyMat k a b
    
-- | Apply to a vector and store the result in another vector
doApply :: (RMatrix a e) =>
    a (k,l) e -> x l e -> y k e -> m ()
doApply a x y =
    checkMatVecMultAdd (numRows a, numCols a) (dim x) (dim y) $
        unsafeDoApply a x y
    
-- | @y := a x + beta y@
doApplyAdd :: (RMatrix a e) =>
    a (k,l) e -> x l e -> e -> y k e -> m ()
doApplyAdd a x beta y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $
        unsafeDoApplyAdd a x beta y

-- | @y := alpha a x@
doSApply :: (RMatrix a e) =>
    e -> a (k,l) e -> x l e -> y k e -> m ()
doSApply alpha a x y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $
        unsafeDoSApply alpha a x y

-- | @y := alpha a x + beta y@    
doSApplyAdd :: (RMatrix a e) =>
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
doSApplyAdd alpha a x beta y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $
        unsafeDoSApplyAdd alpha a x beta y

-- | Apply to a matrix and store the result in another matrix
doApplyMat :: (RMatrix a e) =>
    a (m,k) e -> DMatrix t (k,n) e -> IOMatrix (m,n) e -> m ()
doApplyMat a b c =
    checkMatMatMultAdd (shape a) (shape b) (shape c) $
        unsafeDoApplyMat a b c
    
-- | @c := a b + beta c@
doApplyAddMat :: (RMatrix a e) =>
    a (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> m ()
doApplyAddMat a b beta c =
    checkMatMatMultAdd (shape a) (shape b) (shape c)
        unsafeDoApplyAddMat a b beta c

-- | Scale and apply to a matrix and store the result in another matrix
doSApplyMat :: (RMatrix a e) => 
    e -> a (m,k) e -> DMatrix t (k,n) e -> IOMatrix (m,n) e -> m ()
doSApplyMat alpha a b c =
    checkMatMatMultAdd (shape a) (shape b) (shape c) $
        unsafeDoSApplyMat alpha a b c
        
-- | @c := alpha a b + beta c@
doSApplyAddMat :: (RMatrix a e) =>
    e -> a (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> m ()
doSApplyAddMat alpha a b beta c =
    checkMatMatMultAdd (shape a) (shape b) (shape c)
        unsafeDoSApplyAddMat alpha a b beta c

-- | @x := a x@    
doApply_ :: (RMatrix a e) =>
    a (n,n) e -> IOVector n e -> m ()
doApply_ a x =
    checkSquare (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoApply_ a x

-- | @ b := a b@
doApplyMat_ :: (RMatrix a e) =>
    a (m,m) e -> IOMatrix (m,n) e -> m ()
doApplyMat_ a b =
    checkSquare (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoApplyMat_ a b

-- | @ x := alpha a x@        
doSApply_ :: (RMatrix a e) =>
    e -> a (n,n) e -> IOVector n e -> m ()
doSApply_ alpha a x =
    checkSquare (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoSApply_ alpha a x

-- | @ b := alpha a b@
doSApplyMat_ :: (RMatrix a e) =>
    e -> a (m,m) e -> IOMatrix (m,n) e -> m ()
doSApplyMat_ alpha a b =
    checkSquare (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoSApplyMat_ alpha a b


instance (BLAS3 e) => RMatrix (DMatrix t) e where
    unsafeDoSApplyAdd    = M.gemv
    unsafeDoSApplyAddMat = M.gemm
-}
