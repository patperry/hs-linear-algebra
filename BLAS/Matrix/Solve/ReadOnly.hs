{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Solve.ReadOnly
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Solve.ReadOnly (
    RSolve(..)
    ) where

import BLAS.Matrix.ReadOnly

import Data.Vector.Dense.IO ( DVector )
import Data.Matrix.Dense.IO ( DMatrix )

class RMatrix a e => RSolve a e where
    -- | Solve for a vector
    getSolve :: a (m,n) e -> DVector t m e -> IO (DVector r n e)
    
    -- | Solve for a matrix
    getSolveMat :: a (m,n) e -> DMatrix t (m,k) e -> IO (DMatrix r (n,k) e)
