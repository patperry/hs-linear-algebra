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

import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv )
import BLAS.Matrix.Base ( numRows, numCols )
import BLAS.Matrix.ReadOnly


import Data.Vector.Dense.IO ( DVector, dim )
import Data.Matrix.Dense.IO ( DMatrix, shape )

class RMatrix a e => RSolve a e where
    -- | Solve for a vector
    getSolve :: a (m,n) e -> DVector t m e -> IO (DVector r n e)
    getSolve a y = 
        checkMatVecSolv (numRows a, numCols a) (dim y) $
            unsafeGetSolve a y
    
    -- | Solve for a matrix
    getSolveMat :: a (m,n) e -> DMatrix t (m,k) e -> IO (DMatrix r (n,k) e)
    getSolveMat a b =
        checkMatMatSolv (numRows a, numCols a) (shape b) $
                unsafeGetSolveMat a b

    unsafeGetSolve :: a (m,n) e -> DVector t m e -> IO (DVector r n e)
    unsafeGetSolveMat :: a (m,n) e -> DMatrix t (m,k) e -> IO (DMatrix r (n,k) e)
