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
    getSolve,
    getSolveMat,
    getSSolve,
    getSSolveMat,
    
    -- * In-place operations
    doSolve,
    doSolveMat,
    doSSolve,
    doSSolveMat,
    doSolve_,
    doSolveMat_,
    doSSolve_,
    doSSolveMat_,
    
    -- * Unsafe operations
    unsafeGetSolve,
    unsafeGetSolveMat,
    unsafeGetSSolve,
    unsafeGetSSolveMat,

    ) where

import Data.Matrix.Class
import Data.Tensor.Class
import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv, checkMatVecSolvTo,
    checkMatMatSolvTo, checkSquare )

import Data.Vector.Dense.Class
import Data.Matrix.Dense.Base


unsafeGetSolve :: (MSolve a e m, ReadVector y e m, WriteVector x e m) => 
    a (k,l) e -> y k e -> m (x l e)
unsafeGetSolve a y = do
    x <- newVector_ (numCols a)
    unsafeDoSolve a y x
    return x
{-# INLINE unsafeGetSolve #-}
    
unsafeGetSSolve :: (MSolve a e m, ReadVector y e m, WriteVector x e m) => 
    e -> a (k,l) e -> y k e -> m (x l e)
unsafeGetSSolve alpha a y = do
    x <- newVector_ (numCols a)
    unsafeDoSSolve alpha a y x
    return x
{-# INLINE unsafeGetSSolve #-}
    
unsafeGetSolveMat :: (MSolve a e m, ReadMatrix c e m, WriteMatrix b e m) => 
    a (r,s) e -> c (r,t) e -> m (b (s,t) e)
unsafeGetSolveMat a c = do
    b <- newMatrix_ (numCols a, numCols c)
    unsafeDoSolveMat a c b
    return b
{-# INLINE unsafeGetSolveMat #-}

unsafeGetSSolveMat :: (MSolve a e m, ReadMatrix c e m, WriteMatrix b e m) => 
    e -> a (r,s) e -> c (r,t) e -> m (b (s,t) e)                         
unsafeGetSSolveMat alpha a c = do
    b <- newMatrix_ (numCols a, numCols c)
    unsafeDoSSolveMat alpha a c b
    return b
{-# INLINE unsafeGetSSolveMat #-}

-- | Return @x@ such that @a x = y@.
getSolve :: (MSolve a e m, ReadVector y e m, WriteVector x e m) =>
    a (k,l) e -> y k e -> m (x l e)
getSolve a y = 
    checkMatVecSolv (shape a) (dim y) $
        unsafeGetSolve a y
{-# INLINE getSolve #-}

-- | Return @x@ such that @a x = alpha y@.    
getSSolve :: (MSolve a e m, ReadVector y e m, WriteVector x e m) => 
    e -> a (k,l) e -> y k e -> m (x l e)
getSSolve alpha a y = 
    checkMatVecSolv (shape a) (dim y) $
        unsafeGetSSolve alpha a y
{-# INLINE getSSolve #-}

-- | Return @b@ such that @a b = c@.
getSolveMat :: (MSolve a e m, ReadMatrix c e m, WriteMatrix b e m) => 
    a (r,s) e -> c (r,t) e -> m (b (s,t) e)                     
getSolveMat a c =
    checkMatMatSolv (shape a) (shape c) $
            unsafeGetSolveMat a c
{-# INLINE getSolveMat #-}
            
-- | Return @b@ such that @a b = alpha c@.
getSSolveMat :: (MSolve a e m, ReadMatrix c e m, WriteMatrix b e m) => 
    e -> a (r,s) e -> c (r,t) e -> m (b (s,t) e)                 
getSSolveMat alpha a b =
    checkMatMatSolv (shape a) (shape b) $
            unsafeGetSSolveMat alpha a b
{-# INLINE getSSolveMat #-}

-- | Set @x := a^{-1} y@.
doSolve :: (MSolve a e m, ReadVector y e m, WriteVector x e m) => 
    a (r,s) e -> y r e -> x s e -> m ()                 
doSolve a y x =
    checkMatVecSolvTo (shape a) (dim y) (dim x) $
        unsafeDoSolve a y x
{-# INLINE doSolve #-}
        
-- | Set @b := a^{-1} c@.
doSolveMat :: (MSolve a e m, ReadMatrix c e m, WriteMatrix b e m) => 
    a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()                
doSolveMat a c b =
    checkMatMatSolvTo (shape a) (shape c) (shape b) $
        unsafeDoSolveMat a c b
{-# INLINE doSolveMat #-}
    
-- | Set @x := a^{-1} (alpha y)@.    
doSSolve :: (MSolve a e m, ReadVector y e m, WriteVector x e m) => 
    e -> a (k,l) e -> y k e -> x l e -> m ()
doSSolve alpha a y x =
    checkMatVecSolvTo (shape a) (dim y) (dim x) $
        unsafeDoSSolve alpha a y x
{-# INLINE doSSolve #-}

-- | Set @b := a^{-1} (alpha c)@.
doSSolveMat :: (MSolve a e m, ReadMatrix c e m, WriteMatrix b e m) => 
    e -> a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()                  
doSSolveMat alpha a c b =
    checkMatMatSolvTo (shape a) (shape c) (shape b) $
        unsafeDoSSolveMat alpha a c b
{-# INLINE doSSolveMat #-}

-- | Set @x := a^{-1} x@.
doSolve_ :: (MSolve a e m, ReadVector y e m, WriteVector x e m) => 
    a (k,k) e -> x k e -> m ()
doSolve_ a x =
    checkSquare "doSolve_" (shape a) $
        checkMatVecSolv (shape a) (dim x) $
            unsafeDoSolve_ a x
{-# INLINE doSolve_ #-}

-- | Set @x := a^{-1} (alpha x)@.
doSSolve_ :: (MSolve a e m, WriteVector x e m) => 
    e -> a (k,k) e -> x k e -> m ()
doSSolve_ alpha a x =
    checkSquare ("doSSolve_ " ++ show alpha) (shape a) $
        checkMatVecSolv (shape a) (dim x) $
            unsafeDoSSolve_ alpha a x
{-# INLINE doSSolve_ #-}

-- | Set @b := a^{-1} b@.
doSolveMat_ :: (MSolve a e m, WriteMatrix b e m) => 
    a (k,k) e -> b (k,l) e -> m ()          
doSolveMat_ a b =
    checkSquare "doSolveMat_" (shape a) $
        checkMatMatSolv (shape a) (shape b) $
            unsafeDoSolveMat_ a b
{-# INLINE doSolveMat_ #-}

-- | Set @b := a^{-1} (alpha b)@.
doSSolveMat_ :: (MSolve a e m, WriteMatrix b e m) =>
    e -> a (k,k) e -> b (k,l) e -> m ()          
doSSolveMat_ alpha a b =
    checkSquare ("doSSolveMat_ " ++ show alpha) (shape a) $
        checkMatMatSolv (shape a) (shape b) $
            unsafeDoSSolveMat_ alpha a b
{-# INLINE doSSolveMat_ #-}
