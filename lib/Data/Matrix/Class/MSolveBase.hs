{-# LANGUAGE MultiParamTypeClasses #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class.MSolveBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface for solving matrix systems in a monad.  The
-- matrices can operate via inverse multiplication on immutable dense
-- vectors and matrices.
--

module Data.Matrix.Class.MSolveBase (
    -- * The MSolve type class
    MSolve(..),

    -- * Solving linear systems
    getSolveVector,
    getSolveMatrix,
    getSSolveVector,
    getSSolveMatrix,
    
    -- * In-place operations
    doSolveVector,
    doSolveMatrix,
    doSSolveVector,
    doSSolveMatrix,
    doSolveVector_,
    doSolveMatrix_,
    doSSolveVector_,
    doSSolveMatrix_,
    
    -- * Unsafe operations
    unsafeGetSolveVector,
    unsafeGetSolveMatrix,
    unsafeGetSSolveVector,
    unsafeGetSSolveMatrix,

    ) where

import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv, checkMatVecSolvTo,
    checkMatMatSolvTo, checkSquare )

import Data.Elem.BLAS
import Data.Matrix.Class
import Data.Tensor.Class
import Data.Vector.Dense.Class
import Data.Matrix.Dense.Base


unsafeGetSolveVector :: (MSolve a m, ReadVector y m, WriteVector x m, BLAS2 e) => 
    a (k,l) e -> y k e -> m (x l e)
unsafeGetSolveVector a y = do
    x  <- newVector_ (numCols a)
    unsafeDoSolveVector a y x
    return x
{-# INLINE unsafeGetSolveVector #-}
    
unsafeGetSSolveVector :: (MSolve a m, ReadVector y m, WriteVector x m, BLAS2 e) => 
    e -> a (k,l) e -> y k e -> m (x l e)
unsafeGetSSolveVector alpha a y = do
    x  <- newVector_ (numCols a)
    unsafeDoSSolveVector alpha a y x
    return x
{-# INLINE unsafeGetSSolveVector #-}
    
unsafeGetSolveMatrix :: (MSolve a m, ReadMatrix c m, WriteMatrix b m, BLAS3 e) => 
    a (r,s) e -> c (r,t) e -> m (b (s,t) e)
unsafeGetSolveMatrix a c = do
    b  <- newMatrix_ (numCols a, numCols c)
    unsafeDoSolveMatrix a c b
    return b
{-# INLINE unsafeGetSolveMatrix #-}

unsafeGetSSolveMatrix :: (MSolve a m, ReadMatrix c m, WriteMatrix b m, BLAS3 e) => 
    e -> a (r,s) e -> c (r,t) e -> m (b (s,t) e)                         
unsafeGetSSolveMatrix alpha a c = do
    b  <- newMatrix_ (numCols a, numCols c)
    unsafeDoSSolveMatrix alpha a c b
    return b
{-# INLINE unsafeGetSSolveMatrix #-}

-- | Return @x@ such that @a x = y@.
getSolveVector :: (MSolve a m, ReadVector y m, WriteVector x m, BLAS2 e) =>
    a (k,l) e -> y k e -> m (x l e)
getSolveVector a y = 
    checkMatVecSolv (shape a) (dim y) $
        unsafeGetSolveVector a y
{-# INLINE getSolveVector #-}

-- | Return @x@ such that @a x = alpha y@.    
getSSolveVector :: (MSolve a m, ReadVector y m, WriteVector x m, BLAS2 e) => 
    e -> a (k,l) e -> y k e -> m (x l e)
getSSolveVector alpha a y = 
    checkMatVecSolv (shape a) (dim y) $
        unsafeGetSSolveVector alpha a y
{-# INLINE getSSolveVector #-}

-- | Return @b@ such that @a b = c@.
getSolveMatrix :: (MSolve a m, ReadMatrix c m, WriteMatrix b m, BLAS3 e) => 
    a (r,s) e -> c (r,t) e -> m (b (s,t) e)                     
getSolveMatrix a c =
    checkMatMatSolv (shape a) (shape c) $
            unsafeGetSolveMatrix a c
{-# INLINE getSolveMatrix #-}
            
-- | Return @b@ such that @a b = alpha c@.
getSSolveMatrix :: (MSolve a m, ReadMatrix c m, WriteMatrix b m, BLAS3 e) => 
    e -> a (r,s) e -> c (r,t) e -> m (b (s,t) e)                 
getSSolveMatrix alpha a b =
    checkMatMatSolv (shape a) (shape b) $
            unsafeGetSSolveMatrix alpha a b
{-# INLINE getSSolveMatrix #-}

-- | Set @x := a^{-1} y@.
doSolveVector :: (MSolve a m, ReadVector y m, WriteVector x m, BLAS2 e) => 
    a (r,s) e -> y r e -> x s e -> m ()                 
doSolveVector a y x =
    checkMatVecSolvTo (shape a) (dim y) (dim x) $
        unsafeDoSolveVector a y x
{-# INLINE doSolveVector #-}
        
-- | Set @b := a^{-1} c@.
doSolveMatrix :: (MSolve a m, ReadMatrix c m, WriteMatrix b m, BLAS3 e) => 
    a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()                
doSolveMatrix a c b =
    checkMatMatSolvTo (shape a) (shape c) (shape b) $
        unsafeDoSolveMatrix a c b
{-# INLINE doSolveMatrix #-}
    
-- | Set @x := a^{-1} (alpha y)@.    
doSSolveVector :: (MSolve a m, ReadVector y m, WriteVector x m, BLAS2 e) => 
    e -> a (k,l) e -> y k e -> x l e -> m ()
doSSolveVector alpha a y x =
    checkMatVecSolvTo (shape a) (dim y) (dim x) $
        unsafeDoSSolveVector alpha a y x
{-# INLINE doSSolveVector #-}

-- | Set @b := a^{-1} (alpha c)@.
doSSolveMatrix :: (MSolve a m, ReadMatrix c m, WriteMatrix b m, BLAS3 e) => 
    e -> a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()                  
doSSolveMatrix alpha a c b =
    checkMatMatSolvTo (shape a) (shape c) (shape b) $
        unsafeDoSSolveMatrix alpha a c b
{-# INLINE doSSolveMatrix #-}

-- | Set @x := a^{-1} x@.
doSolveVector_ :: (MSolve a m, ReadVector y m, WriteVector x m, BLAS2 e) => 
    a (k,k) e -> x k e -> m ()
doSolveVector_ a x =
    checkSquare "doSolveVector_" (shape a) $
        checkMatVecSolv (shape a) (dim x) $
            unsafeDoSolveVector_ a x
{-# INLINE doSolveVector_ #-}

-- | Set @x := a^{-1} (alpha x)@.
doSSolveVector_ :: (MSolve a m, WriteVector x m, BLAS2 e) => 
    e -> a (k,k) e -> x k e -> m ()
doSSolveVector_ alpha a x =
    checkSquare ("doSSolveVector_ " ++ show alpha) (shape a) $
        checkMatVecSolv (shape a) (dim x) $
            unsafeDoSSolveVector_ alpha a x
{-# INLINE doSSolveVector_ #-}

-- | Set @b := a^{-1} b@.
doSolveMatrix_ :: (MSolve a m, WriteMatrix b m, BLAS3 e) => 
    a (k,k) e -> b (k,l) e -> m ()          
doSolveMatrix_ a b =
    checkSquare "doSolveMatrix_" (shape a) $
        checkMatMatSolv (shape a) (shape b) $
            unsafeDoSolveMatrix_ a b
{-# INLINE doSolveMatrix_ #-}

-- | Set @b := a^{-1} (alpha b)@.
doSSolveMatrix_ :: (MSolve a m, WriteMatrix b m, BLAS3 e) =>
    e -> a (k,k) e -> b (k,l) e -> m ()          
doSSolveMatrix_ alpha a b =
    checkSquare ("doSSolveMatrix_ " ++ show alpha) (shape a) $
        checkMatMatSolv (shape a) (shape b) $
            unsafeDoSSolveMatrix_ alpha a b
{-# INLINE doSSolveMatrix_ #-}
