{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.ReadOnly
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.ReadOnly (
    RMatrix(..),
    
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
    
    ) where

import BLAS.Elem ( BLAS1, BLAS3 )
import BLAS.Internal ( checkMatVecMult, checkMatMatMult, checkMatVecMultAdd,
    checkMatMatMultAdd, checkSquare )
import qualified BLAS.Matrix.Base as Base
import Data.Vector.Dense.Internal
import Data.Vector.Dense.Operations
import qualified Data.Vector.Dense.Operations as V
import Data.Matrix.Dense.Internal
import Data.Matrix.Dense.Operations
import qualified Data.Matrix.Dense.Operations as M

import Unsafe.Coerce

-- | Minimal complete definition: (unsafeDoSApplyAdd, unsafeDoSApplyAddMat)
-- or (unsafeDoSApply, unsafeDoSApplyMat)
class (Base.Matrix a, BLAS1 e) => RMatrix a e where
    unsafeDoApply :: a (m,n) e -> DVector t n e -> IOVector m e -> IO ()
    unsafeDoApply a x y =
        unsafeDoApplyAdd a x 0 y
        
    unsafeDoApplyAdd :: a (m,n) e -> DVector t n e -> e -> IOVector m e -> IO ()
    unsafeDoApplyAdd a x beta y =
        unsafeDoSApplyAdd 1 a x beta y

    unsafeDoSApply :: e -> a (m,n) e -> DVector t n e -> IOVector m e -> IO ()
    unsafeDoSApply alpha a x y =
        unsafeDoSApplyAdd alpha a x 0 y

    unsafeDoApplyMat :: a (m,k) e -> DMatrix t (k,n) e -> IOMatrix (m,n) e -> IO ()
    unsafeDoApplyMat a b c =
        unsafeDoApplyAddMat a b 0 c
        
    unsafeDoApplyAddMat :: a (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> IO ()
    unsafeDoApplyAddMat a b beta c =
        unsafeDoSApplyAddMat 1 a b beta c

    unsafeDoSApplyMat :: e -> a (m,k) e -> DMatrix t (k,n) e -> IOMatrix (m,n) e -> IO ()
    unsafeDoSApplyMat alpha a b c =
        unsafeDoSApplyAddMat alpha a b 0 c
    
    unsafeDoApply_ :: a (n,n) e -> IOVector n e -> IO ()
    unsafeDoApply_ a x =
        unsafeDoSApply_ 1 a x
        
    unsafeDoApplyMat_ :: a (m,m) e -> IOMatrix (m,n) e -> IO ()
    unsafeDoApplyMat_ a b =
        unsafeDoSApplyMat_ 1 a b

    unsafeDoSApply_ :: e -> a (n,n) e -> IOVector n e -> IO ()
    unsafeDoSApply_ alpha a x = do
        y <- newVector_ (dim x)
        unsafeDoSApply alpha a x y
        unsafeCopyVector x y

    unsafeDoSApplyMat_ :: e -> a (m,m) e -> IOMatrix (m,n) e -> IO ()
    unsafeDoSApplyMat_ alpha a b = do
        c <- newMatrix_ (shape b)
        unsafeDoSApplyMat alpha a b c
        unsafeCopyMatrix b c
        
    unsafeDoSApplyAdd :: e -> a (m,n) e -> DVector t n e -> e -> IOVector m e -> IO ()
    unsafeDoSApplyAdd alpha a x beta y = do
        y' <- unsafeGetApply a x
        V.scaleBy beta y
        V.unsafeAxpy alpha y' y
        
    unsafeDoSApplyAddMat :: e -> a (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> IO ()
    unsafeDoSApplyAddMat alpha a b beta c = do
        c' <- unsafeGetApplyMat a b
        M.scaleBy beta c
        M.axpy alpha c' c


unsafeGetApply :: (RMatrix a e) =>
    a (m,n) e -> DVector t n e -> IO (DVector r m e)
unsafeGetApply = unsafeGetSApply 1

unsafeGetApplyMat :: (RMatrix a e) =>
    a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
unsafeGetApplyMat = unsafeGetSApplyMat 1

unsafeGetSApply :: (RMatrix a e) =>
    e -> a (m,n) e -> DVector t n e -> IO (DVector r m e)
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
    a (m,n) e -> DVector t n e -> IO (DVector r m e)
getApply a x =
    checkMatVecMult (shape a) (dim x) $ 
        unsafeGetApply a x

-- | Apply to a matrix
getApplyMat :: (RMatrix a e) =>
    a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
getApplyMat a b =
    checkMatMatMult (shape a) (shape b) $
        unsafeGetApplyMat a b

-- | Scale and apply to a vector
getSApply :: (RMatrix a e) =>
    e -> a (m,n) e -> DVector t n e -> IO (DVector r m e)
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
    a (m,n) e -> DVector t n e -> IOVector m e -> IO ()
doApply a x y =
    checkMatVecMultAdd (numRows a, numCols a) (dim x) (dim y) $
        unsafeDoApply a x y
    
-- | @y := a x + beta y@
doApplyAdd :: (RMatrix a e) =>
    a (m,n) e -> DVector t n e -> e -> IOVector m e -> IO ()
doApplyAdd a x beta y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $
        unsafeDoApplyAdd a x beta y

-- | @y := alpha a x@
doSApply :: (RMatrix a e) =>
    e -> a (m,n) e -> DVector t n e -> IOVector m e -> IO ()
doSApply alpha a x y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $
        unsafeDoSApply alpha a x y

-- | @y := alpha a x + beta y@    
doSApplyAdd :: (RMatrix a e) =>
    e -> a (m,n) e -> DVector t n e -> e -> IOVector m e -> IO ()
doSApplyAdd alpha a x beta y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $
        unsafeDoSApplyAdd alpha a x beta y

-- | Apply to a matrix and store the result in another matrix
doApplyMat :: (RMatrix a e) =>
    a (m,k) e -> DMatrix t (k,n) e -> IOMatrix (m,n) e -> IO ()
doApplyMat a b c =
    checkMatMatMultAdd (shape a) (shape b) (shape c) $
        unsafeDoApplyMat a b c
    
-- | @c := a b + beta c@
doApplyAddMat :: (RMatrix a e) =>
    a (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> IO ()
doApplyAddMat a b beta c =
    checkMatMatMultAdd (shape a) (shape b) (shape c)
        unsafeDoApplyAddMat a b beta c

-- | Scale and apply to a matrix and store the result in another matrix
doSApplyMat :: (RMatrix a e) => 
    e -> a (m,k) e -> DMatrix t (k,n) e -> IOMatrix (m,n) e -> IO ()
doSApplyMat alpha a b c =
    checkMatMatMultAdd (shape a) (shape b) (shape c) $
        unsafeDoSApplyMat alpha a b c
        
-- | @c := alpha a b + beta c@
doSApplyAddMat :: (RMatrix a e) =>
    e -> a (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> IO ()
doSApplyAddMat alpha a b beta c =
    checkMatMatMultAdd (shape a) (shape b) (shape c)
        unsafeDoSApplyAddMat alpha a b beta c

-- | @x := a x@    
doApply_ :: (RMatrix a e) =>
    a (n,n) e -> IOVector n e -> IO ()
doApply_ a x =
    checkSquare (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoApply_ a x

-- | @ b := a b@
doApplyMat_ :: (RMatrix a e) =>
    a (m,m) e -> IOMatrix (m,n) e -> IO ()
doApplyMat_ a b =
    checkSquare (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoApplyMat_ a b

-- | @ x := alpha a x@        
doSApply_ :: (RMatrix a e) =>
    e -> a (n,n) e -> IOVector n e -> IO ()
doSApply_ alpha a x =
    checkSquare (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoSApply_ alpha a x

-- | @ b := alpha a b@
doSApplyMat_ :: (RMatrix a e) =>
    e -> a (m,m) e -> IOMatrix (m,n) e -> IO ()
doSApplyMat_ alpha a b =
    checkSquare (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoSApplyMat_ alpha a b


instance (BLAS3 e) => RMatrix (DMatrix t) e where
    unsafeDoSApplyAdd    = M.gemv
    unsafeDoSApplyAddMat = M.gemm
