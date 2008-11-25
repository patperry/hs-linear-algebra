{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Apply.Read
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.MMatrix (
    -- * Getting rows and columns
    getRow,
    getCol,
    getRows,
    getCols,
    getRows',
    getCols',
    
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

    -- * The MMatrix type class
    MMatrix(..),

    -- * Unsafe operations
    unsafeGetApply,
    unsafeDoApply,
    unsafeDoApply_,

    unsafeGetApplyMat,
    unsafeDoApplyMat,
    unsafeDoApplyMat_,

    ) where

import BLAS.Internal( checkSquare, checkMatVecMult, checkMatVecMultAdd,
    checkMatMatMult, checkMatMatMultAdd, checkedRow, checkedCol )
import BLAS.UnsafeIOToM
import BLAS.Matrix.Shaped

import Data.Vector.Dense.Class
import Data.Matrix.Dense.Class.Internal



-- | Get the given row in a matrix.
getRow :: (MMatrix a e m, WriteVector x e m) => a (k,l) e -> Int -> m (x l e)
getRow a = checkedRow (shape a) (unsafeGetRow a)

-- | Get the given column in a matrix.
getCol :: (MMatrix a e m, WriteVector x e m) => a (k,l) e -> Int -> m (x k e)
getCol a = checkedCol (shape a) (unsafeGetCol a)

-- | Get a lazy list the row vectors in the matrix.  See also "getRows'".
getRows :: (MMatrix a e m, WriteVector x e m) => 
    a (k,l) e -> m [x l e]
getRows = unsafeInterleaveM . getRows'

-- | Get a lazy list of the column vectors in the matrix.  See also "getCols'".
getCols :: (MMatrix a e m, WriteVector x e m) => 
    a (k,l) e -> m [x k e]
getCols = unsafeInterleaveM . getCols'

-- | Get a strict list the row vectors in the matrix.  See also "getRows".
getRows' :: (MMatrix a e m, WriteVector x e m) => a (k,l) e -> m [x l e]
getRows' a = mapM (unsafeGetRow a) [0..numRows a - 1]

-- | Get a strict list of the column vectors in the matrix.  See also "getCols".
getCols' :: (MMatrix a e m, WriteVector x e m) => a (k,l) e -> m [x k e]
getCols' a = mapM (unsafeGetCol a) [0..numCols a - 1]

-- | Scale and apply to a vector
getSApply :: (MMatrix a e m, ReadVector x e m, WriteVector y e m) =>
    e -> a (k,l) e -> x l e -> m (y k e)
getSApply k a x =
    checkMatVecMult (shape a) (dim x) $ 
        unsafeGetSApply k a x

-- | Scale and apply to a matrix
getSApplyMat :: (MMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    e -> a (r,s) e -> b (s,t) e -> m (c (r,t) e)
getSApplyMat k a b =
    checkMatMatMult (shape a) (shape b) $
        unsafeGetSApplyMat k a b
    
-- | @y := alpha a x + beta y@    
doSApplyAdd :: (MMatrix a e m, ReadVector x e m, WriteVector y e m) =>
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
doSApplyAdd alpha a x beta y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $
        unsafeDoSApplyAdd alpha a x beta y

-- | @c := alpha a b + beta c@
doSApplyAddMat :: (MMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
doSApplyAddMat alpha a b beta c =
    checkMatMatMultAdd (shape a) (shape b) (shape c)
        unsafeDoSApplyAddMat alpha a b beta c

-- | Apply to a vector
getApply :: (MMatrix a e m, ReadVector x e m, WriteVector y e m) =>
    a (k,l) e -> x l e -> m (y k e)
getApply a x =
    checkMatVecMult (shape a) (dim x) $ do
        unsafeGetApply a x

-- | Apply to a matrix
getApplyMat :: (MMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    a (r,s) e -> b (s,t) e -> m (c (r,t) e)
getApplyMat a b =
    checkMatMatMult (shape a) (shape b) $
        unsafeGetApplyMat a b

-- | @ x := alpha a x@        
doSApply_ :: (MMatrix a e m, WriteVector y e m) =>
    e -> a (n,n) e -> y n e -> m ()
doSApply_ alpha a x =
    checkSquare (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoSApply_ alpha a x

-- | @ b := alpha a b@
doSApplyMat_ :: (MMatrix a e m, WriteMatrix b e m) =>
    e -> a (s,s) e -> b (s,t) e -> m ()
doSApplyMat_ alpha a b =
    checkSquare (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoSApplyMat_ alpha a b

unsafeGetApply :: (MMatrix a e m, ReadVector x e m, WriteVector y e m) =>
    a (k,l) e -> x l e -> m (y k e)
unsafeGetApply = unsafeGetSApply 1

unsafeGetApplyMat :: (MMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    a (r,s) e -> b (s,t) e -> m (c (r,t) e)
unsafeGetApplyMat = unsafeGetSApplyMat 1

-- | Apply to a vector and store the result in another vector
doApply :: (MMatrix a e m, ReadVector x e m, WriteVector y e m) =>
    a (k,l) e -> x l e -> y k e -> m ()
doApply a x y =
    checkMatVecMultAdd (numRows a, numCols a) (dim x) (dim y) $
        unsafeDoApply a x y

-- | Apply to a matrix and store the result in another matrix
doApplyMat :: (MMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    a (r,s) e -> b (s,t) e -> c (r,t) e -> m ()
doApplyMat a b c =
    checkMatMatMultAdd (shape a) (shape b) (shape c) $
        unsafeDoApplyMat a b c
        
unsafeDoApply :: (MMatrix a e m, ReadVector x e m, WriteVector y e m) =>
    a (k,l) e -> x l e -> y k e -> m ()
unsafeDoApply a x y = unsafeDoSApplyAdd 1 a x 0 y

unsafeDoApplyMat :: (MMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    a (r,s) e -> b (s,t) e -> c (r,t) e -> m ()
unsafeDoApplyMat a b c = unsafeDoSApplyAddMat 1 a b 0 c

-- | @x := a x@    
doApply_ :: (MMatrix a e m, WriteVector y e m) =>
    a (n,n) e -> y n e -> m ()
doApply_ a x =
    checkSquare (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoApply_ a x

-- | @ b := a b@
doApplyMat_ :: (MMatrix a e m, WriteMatrix b e m) =>
    a (s,s) e -> b (s,t) e -> m ()
doApplyMat_ a b =
    checkSquare (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoApplyMat_ a b
  
unsafeDoApply_ :: (MMatrix a e m, WriteVector y e m) => 
    a (n,n) e -> y n e -> m ()
unsafeDoApply_ a x =
    unsafeDoSApply_ 1 a x

unsafeDoApplyMat_ :: (MMatrix a e m, WriteMatrix b e m) =>
    a (s,s) e -> b (s,t) e -> m ()
unsafeDoApplyMat_ a b = 
    unsafeDoSApplyMat_ 1 a b
