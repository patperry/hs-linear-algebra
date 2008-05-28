{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.ReadOnly
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.ReadOnly (
    RMatrix(..)
    ) where

import BLAS.Elem ( BLAS3 )
import qualified BLAS.Matrix.Base as Base
import Data.Vector.Dense.Internal
import Data.Matrix.Dense.Internal
import qualified Data.Matrix.Dense.Operations as M

class Base.Matrix a => RMatrix a e where
    -- | Apply to a vector
    getApply :: a (m,n) e -> DVector t n e -> IO (DVector r m e)
    
    -- | Apply to a matrix
    getApplyMat :: a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)

instance (BLAS3 e) => RMatrix (DMatrix t) e where
    getApply    = M.getApply
    getApplyMat = M.getApplyMat
