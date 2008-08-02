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
    ISolve(..)
    ) where

import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv )
import BLAS.Matrix.Immutable
import BLAS.Matrix.Base ( numRows, numCols )

import Data.Vector.Dense ( Vector, dim )
import Data.Matrix.Dense ( Matrix, shape )

infixr 7 <\>, <\\>

class IMatrix a e => ISolve a e where
    -- | Solve for a vector
    (<\>) :: a (m,n) e -> Vector m e -> Vector n e
    (<\>) a y = 
        checkMatVecSolv (numRows a, numCols a) (dim y) $
            unsafeSolve a y
    
    -- | Solve for a matrix
    (<\\>) :: a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
    (<\\>) a b =
        checkMatMatSolv (numRows a, numCols a) (shape b) $
            unsafeSolveMat a b

    unsafeSolve :: a (m,n) e -> Vector m e -> Vector n e
    unsafeSolveMat :: a (m,n) e -> Matrix (m,k) e -> Matrix (n,k) e
    