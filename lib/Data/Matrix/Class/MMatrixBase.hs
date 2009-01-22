{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class.MMatrixBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface for mutable matrices. The type class associates a
-- matrix with a monad type in which operations can be perfomred.  The
-- matrices provide access to rows and columns, and can operate via
-- multiplication on dense vectors and matrices.  
--

module Data.Matrix.Class.MMatrixBase (
    -- * The MMatrix type class
    MMatrix(..),

    -- * Getting rows and columns
    getRow,
    getCol,
    getRows,
    getCols,
    getRows',
    getCols',
    
    -- * Matrix and vector multiplication
    getApplyVector,
    getSApplyVector,
    
    getApplyMatrix,
    getSApplyMatrix,

    -- * In-place multiplication
    doApplyVector,
    doSApplyAddVector,
    doApplyVector_,
    doSApplyVector_,
    
    doApplyMatrix,
    doSApplyAddMatrix,
    doApplyMatrix_,
    doSApplyMatrix_,

    -- * Unsafe operations
    unsafeGetApplyVector,
    unsafeDoApplyVector,
    unsafeDoApplyVector_,

    unsafeGetApplyMatrix,
    unsafeDoApplyMatrix,
    unsafeDoApplyMatrix_,

    ) where

import BLAS.Internal( checkSquare, checkMatVecMult, checkMatVecMultAdd,
    checkMatMatMult, checkMatMatMultAdd, checkedRow, checkedCol )

import Data.Elem.BLAS
import Data.Matrix.Class
import Data.Tensor.Class

import Data.Vector.Dense.Class
import Data.Matrix.Dense.Base


-- | Get the given row in a matrix.
getRow :: (MMatrix a m, WriteVector x m, Elem e) => a (n,p) e -> Int -> m (x p e)
getRow a = checkedRow (shape a) (unsafeGetRow a)
{-# INLINE getRow #-}

-- | Get the given column in a matrix.
getCol :: (MMatrix a m, WriteVector x m, Elem e) => a (n,p) e -> Int -> m (x n e)
getCol a = checkedCol (shape a) (unsafeGetCol a)
{-# INLINE getCol #-}

-- | Get a strict list the row vectors in the matrix.
getRows' :: (MMatrix a m, WriteVector x m, Elem e) => a (n,p) e -> m [x p e]
getRows' a = mapM (unsafeGetRow a) [0..numRows a - 1]
{-# INLINE getRows' #-}

-- | Get a strict list of the column vectors in the matrix.
getCols' :: (MMatrix a m, WriteVector x m, Elem e) => a (n,p) e -> m [x n e]
getCols' a = mapM (unsafeGetCol a) [0..numCols a - 1]
{-# INLINE getCols' #-}

-- | Scale and apply to a vector
getSApplyVector :: (MMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> a (n,p) e -> x p e -> m (y n e)
getSApplyVector k a x =
    checkMatVecMult (shape a) (dim x) $ 
        unsafeGetSApplyVector k a x
{-# INLINE getSApplyVector #-}

-- | Scale and apply to a matrix
getSApplyMatrix :: (MMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
    e -> a (n,p) e -> b (p,q) e -> m (c (n,q) e)
getSApplyMatrix k a b =
    checkMatMatMult (shape a) (shape b) $
        unsafeGetSApplyMatrix k a b
{-# INLINE getSApplyMatrix #-}

-- | @y := alpha a x + beta y@    
doSApplyAddVector :: (MMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> a (n,p) e -> x p e -> e -> y n e -> m ()
doSApplyAddVector alpha a x beta y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $
        unsafeDoSApplyAddVector alpha a x beta y
{-# INLINE doSApplyAddVector #-}

-- | @c := alpha a b + beta c@
doSApplyAddMatrix :: (MMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
    e -> a (n,p) e -> b (p,q) e -> e -> c (n,q) e -> m ()
doSApplyAddMatrix alpha a b beta c =
    checkMatMatMultAdd (shape a) (shape b) (shape c)
        unsafeDoSApplyAddMatrix alpha a b beta c
{-# INLINE doSApplyAddMatrix #-}

-- | ApplyVector to a vector
getApplyVector :: (MMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    a (n,p) e -> x p e -> m (y n e)
getApplyVector a x =
    checkMatVecMult (shape a) (dim x) $ do
        unsafeGetApplyVector a x
{-# INLINE getApplyVector #-}

-- | ApplyVector to a matrix
getApplyMatrix :: (MMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
    a (n,p) e -> b (p,q) e -> m (c (n,q) e)
getApplyMatrix a b =
    checkMatMatMult (shape a) (shape b) $
        unsafeGetApplyMatrix a b
{-# INLINE getApplyMatrix #-}

-- | @ x := alpha a x@        
doSApplyVector_ :: (MMatrix a m, WriteVector y m, BLAS2 e) =>
    e -> a (n,n) e -> y n e -> m ()
doSApplyVector_ alpha a x =
    checkSquare ("doSApplyVector_ " ++ show alpha) (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoSApplyVector_ alpha a x
{-# INLINE doSApplyVector_ #-}

-- | @ b := alpha a b@
doSApplyMatrix_ :: (MMatrix a m, WriteMatrix b m, BLAS3 e) =>
    e -> a (n,n) e -> b (n,p) e -> m ()
doSApplyMatrix_ alpha a b =
    checkSquare ("doSApplyMatrix_ " ++ show alpha) (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoSApplyMatrix_ alpha a b
{-# INLINE doSApplyMatrix_ #-}

unsafeGetApplyVector :: (MMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    a (n,p) e -> x p e -> m (y n e)
unsafeGetApplyVector = unsafeGetSApplyVector 1
{-# INLINE unsafeGetApplyVector #-}

unsafeGetApplyMatrix :: (MMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
    a (n,p) e -> b (p,q) e -> m (c (n,q) e)
unsafeGetApplyMatrix = unsafeGetSApplyMatrix 1
{-# INLINE unsafeGetApplyMatrix #-}

-- | ApplyVector to a vector and store the result in another vector
doApplyVector :: (MMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    a (n,p) e -> x p e -> y n e -> m ()
doApplyVector a x y =
    checkMatVecMultAdd (numRows a, numCols a) (dim x) (dim y) $
        unsafeDoApplyVector a x y
{-# INLINE doApplyVector #-}

-- | ApplyVector to a matrix and store the result in another matrix
doApplyMatrix :: (MMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
    a (n,p) e -> b (p,q) e -> c (n,q) e -> m ()
doApplyMatrix a b c =
    checkMatMatMultAdd (shape a) (shape b) (shape c) $
        unsafeDoApplyMatrix a b c
{-# INLINE doApplyMatrix #-}

unsafeDoApplyVector :: (MMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    a (n,p) e -> x p e -> y n e -> m ()
unsafeDoApplyVector a x y = unsafeDoSApplyAddVector 1 a x 0 y
{-# INLINE unsafeDoApplyVector #-}

unsafeDoApplyMatrix :: (MMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
    a (n,p) e -> b (p,q) e -> c (n,q) e -> m ()
unsafeDoApplyMatrix a b c = unsafeDoSApplyAddMatrix 1 a b 0 c
{-# INLINE unsafeDoApplyMatrix #-}

-- | @x := a x@    
doApplyVector_ :: (MMatrix a m, WriteVector y m, BLAS2 e) =>
    a (n,n) e -> y n e -> m ()
doApplyVector_ a x =
    checkSquare "doApplyVector_" (shape a) $
        checkMatVecMult (shape a) (dim x) $
            unsafeDoApplyVector_ a x
{-# INLINE doApplyVector_ #-}

-- | @ b := a b@
doApplyMatrix_ :: (MMatrix a m, WriteMatrix b m, BLAS3 e) =>
    a (n,n) e -> b (n,p) e -> m ()
doApplyMatrix_ a b =
    checkSquare "doApplyMatrix_" (shape a) $
        checkMatMatMult (shape a) (shape b) $
            unsafeDoApplyMatrix_ a b
{-# INLINE doApplyMatrix_ #-}
  
unsafeDoApplyVector_ :: (MMatrix a m, WriteVector y m, BLAS2 e) => 
    a (n,n) e -> y n e -> m ()
unsafeDoApplyVector_ a x =
    unsafeDoSApplyVector_ 1 a x
{-# INLINE unsafeDoApplyVector_ #-}

unsafeDoApplyMatrix_ :: (MMatrix a m, WriteMatrix b m, BLAS3 e) =>
    a (n,n) e -> b (n,p) e -> m ()
unsafeDoApplyMatrix_ a b = 
    unsafeDoSApplyMatrix_ 1 a b
{-# INLINE unsafeDoApplyMatrix_ #-}
