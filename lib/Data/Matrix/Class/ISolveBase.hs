{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class.ISolveBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface for solving immutable matrix systems.  The
-- matrices can operate via inverse multiplication on immutable dense
-- vectors and matrices.
--

module Data.Matrix.Class.ISolveBase (
    -- * The IMatrix type class
    ISolve(..),
    
    -- * Solving linear systems
    (<\>),
    (<\\>),
    ssolve,
    ssolveMat,
    
    ) where

import Data.Elem.BLAS
import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv )
import Data.Matrix.Class
import Data.Matrix.Class.MSolveBase

import Data.Vector.Dense ( Vector, dim, (*>) )
import Data.Vector.Dense.ST ( runSTVector )
import Data.Matrix.Dense ( Matrix, shape )
import Data.Matrix.Dense.ST ( runSTMatrix )
import Data.Matrix.Tri.Internal

infixr 7 <\>, <\\>

-- | A type class for immutable matrices with inverses.  The member
-- functions of the type class do not perform any checks on the validity
-- of shapes or indices, so in general their safe counterparts should be
-- preferred.
class (MatrixShaped a e, BLAS3 e) => ISolve a e where
    unsafeSolve :: a (m,n) e -> Vector m e -> Vector n e
    unsafeSolve = unsafeSSolve 1
    {-# INLINE unsafeSolve #-}
    
    unsafeSolveMat :: a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
    unsafeSolveMat = unsafeSSolveMat 1
    {-# INLINE unsafeSolveMat #-}

    unsafeSSolve :: e -> a (m,n) e -> Vector m e -> Vector n e
    
    unsafeSSolveMat :: e -> a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e

-- | Solve for a vector.
(<\>) :: (ISolve a e) => a (m,n) e -> Vector m e -> Vector n e
(<\>) a y =
    checkMatVecSolv (shape a) (dim y) $
        unsafeSolve a y
{-# INLINE (<\>) #-}

-- | Solve for a matrix.
(<\\>) :: (ISolve a e) => a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
(<\\>) a b =
    checkMatMatSolv (shape a) (shape b) $
        unsafeSolveMat a b
{-# INLINE (<\\>) #-}

-- | Solve for a vector and scale.
-- @ssolve k a y@ is equal to @a \<\\> (k *> y)@ but is often faster.
ssolve :: (ISolve a e) => e -> a (m,n) e -> Vector m e -> Vector n e
ssolve alpha a y =
    checkMatVecSolv (shape a) (dim y) $
        unsafeSSolve alpha a y
{-# INLINE ssolve #-}

-- | Solve for a matrix and scale.
-- @ssolveMat k a c@ is equal to @a \<\\\\> (k *> c)@ but is often faster.
ssolveMat :: (ISolve a e) => e -> a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
ssolveMat alpha a b =
    checkMatMatSolv (shape a) (shape b) $
        unsafeSSolveMat alpha a b
{-# INLINE ssolveMat #-}

instance (BLAS3 e) => ISolve (Tri Matrix) e where
    unsafeSSolve    alpha a y = runSTVector $ unsafeGetSSolve    alpha a y
    {-# INLINE unsafeSSolve #-}
    unsafeSSolveMat alpha a c = runSTMatrix $ unsafeGetSSolveMat alpha a c
    {-# INLINE unsafeSSolveMat #-}

{-# RULES
"scale.solve/ssolve"       forall k a y. a <\>  (k *> y) = ssolve k a y
"scale.solveMat/ssolveMat" forall k a c. a <\\> (k *> c) = ssolveMat k a c
  #-}
