{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class.ISolve
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Class.ISolve (
    ISolve(..),
    (<\>),
    (<\\>),
    ssolve,
    ssolveMat,
    
    ) where

import Data.Elem.BLAS
import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv )
import Data.Matrix.Class
import Data.Matrix.Class.MSolve

import Data.Vector.Dense ( Vector, dim )
import Data.Vector.Dense.ST ( runSTVector )
import Data.Matrix.Dense ( Matrix, shape )
import Data.Matrix.Dense.ST ( runSTMatrix )
import Data.Matrix.Tri.Internal

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

instance (BLAS3 e) => ISolve (Tri Matrix) e where
    unsafeSSolve    alpha a y = runSTVector $ unsafeGetSSolve    alpha a y
    unsafeSSolveMat alpha a c = runSTMatrix $ unsafeGetSSolveMat alpha a c

