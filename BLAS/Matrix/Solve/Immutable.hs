{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Solve.Immutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Solve.Immutable (
    ISolve(..),
    (<\>),
    (<\\>),
    ssolve,
    ssolveMat,
    
    ) where

import BLAS.Elem
import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv )
import BLAS.Matrix.Shaped

import Data.Vector.Dense ( Vector, dim )
import Data.Matrix.Dense ( Matrix, shape )

infixr 7 <\>, <\\>

class (MatrixShaped a e, BLAS1 e) => ISolve a e where
    unsafeSolve :: a (m,n) e -> Vector m e -> Vector n e
    unsafeSolve = unsafeSSolve 1
    
    unsafeSolveMat :: a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
    unsafeSolveMat = unsafeSSolveMat 1

    unsafeSSolve :: e -> a (m,n) e -> Vector m e -> Vector n e
    
    unsafeSSolveMat :: e -> a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e


-- | Solve for a vector
(<\>) :: (ISolve a e) => a (m,n) e -> Vector m e -> Vector n e
(<\>) a y =
    checkMatVecSolv (shape a) (dim y) $
        unsafeSolve a y

-- | Solve for a matrix
(<\\>) :: (ISolve a e) => a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
(<\\>) a b =
    checkMatMatSolv (shape a) (shape b) $
        unsafeSolveMat a b

-- | Solve for a vector and scale
ssolve :: (ISolve a e) => e -> a (m,n) e -> Vector m e -> Vector n e
ssolve alpha a y =
    checkMatVecSolv (shape a) (dim y) $
        unsafeSSolve alpha a y

-- | Solve for a matrix and scale
ssolveMat :: (ISolve a e) => e -> a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
ssolveMat alpha a b =
    checkMatMatSolv (shape a) (shape b) $
        unsafeSSolveMat alpha a b

