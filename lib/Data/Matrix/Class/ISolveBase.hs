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
    solveVector,
    solveMatrix,
    ssolveVector,
    ssolveMatrix,
    
    ) where

import Data.Elem.BLAS
import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv )
import Data.Matrix.Class
import Data.Matrix.Class.MSolveBase

import Data.Vector.Dense ( Vector, dim )
import Data.Vector.Dense.ST ( runSTVector )
import Data.Matrix.Dense ( Matrix, shape )
import Data.Matrix.Dense.ST ( runSTMatrix )
import Data.Matrix.TriBase

-- | A type class for immutable matrices with inverses.  The member
-- functions of the type class do not perform any checks on the validity
-- of shapes or indices, so in general their safe counterparts should be
-- preferred.
class (MatrixShaped a) => ISolve a where
    unsafeSolveVector :: (BLAS3 e) 
                   => a (m,n) e -> Vector m e -> Vector n e
    unsafeSolveVector = unsafeSSolveVector 1
    {-# INLINE unsafeSolveVector #-}
    
    unsafeSolveMatrix :: (BLAS3 e) 
                   => a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
    unsafeSolveMatrix = unsafeSSolveMatrix 1
    {-# INLINE unsafeSolveMatrix #-}

    unsafeSSolveVector :: (BLAS3 e)
                    => e -> a (m,n) e -> Vector m e -> Vector n e
    
    unsafeSSolveMatrix :: (BLAS3 e)
                    => e -> a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e

-- | SolveVector for a vector.
solveVector :: (ISolve a, BLAS3 e) => a (m,n) e -> Vector m e -> Vector n e
solveVector a y =
    checkMatVecSolv (shape a) (dim y) $
        unsafeSolveVector a y
{-# INLINE solveVector #-}

-- | SolveVector for a matrix.
solveMatrix :: (ISolve a, BLAS3 e) => a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
solveMatrix a b =
    checkMatMatSolv (shape a) (shape b) $
        unsafeSolveMatrix a b
{-# INLINE solveMatrix #-}

-- | SolveVector for a vector and scale.
-- @ssolveVector k a y@ is equal to @a `solveVector` (k *> y)@ but is often faster.
ssolveVector :: (ISolve a, BLAS3 e) => e -> a (m,n) e -> Vector m e -> Vector n e
ssolveVector alpha a y =
    checkMatVecSolv (shape a) (dim y) $
        unsafeSSolveVector alpha a y
{-# INLINE ssolveVector #-}

-- | SolveVector for a matrix and scale.
-- @ssolveMatrix k a c@ is equal to @a `solveMatrix` (k *> c)@ but is often faster.
ssolveMatrix :: (ISolve a, BLAS3 e) => e -> a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
ssolveMatrix alpha a b =
    checkMatMatSolv (shape a) (shape b) $
        unsafeSSolveMatrix alpha a b
{-# INLINE ssolveMatrix #-}

instance ISolve (Tri Matrix) where
    unsafeSSolveVector alpha a y = runSTVector $ unsafeGetSSolveVector alpha a y
    {-# INLINE unsafeSSolveVector #-}
    unsafeSSolveMatrix alpha a c = runSTMatrix $ unsafeGetSSolveMatrix alpha a c
    {-# INLINE unsafeSSolveMatrix #-}
