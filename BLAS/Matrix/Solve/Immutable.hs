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

import BLAS.Matrix.Immutable

import Data.Vector.Dense ( Vector )
import Data.Matrix.Dense ( Matrix )

infixl 7 <\>, <\\>

class IMatrix a e => ISolve a e where
    -- | Solve for a vector
    (<\>) :: a (m,n) e -> Vector m e -> Vector n e
    
    -- | Solve for a matrix
    (<\\>) :: a (m,k) e -> Matrix (m,n) e -> Matrix (k,n) e

