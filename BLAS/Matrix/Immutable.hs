{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Immutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Immutable (
    IMatrix(..)
    ) where

import BLAS.Access
import BLAS.Elem ( BLAS3 )
import qualified BLAS.Matrix.Base as Base
import Data.Vector.Dense
import Data.Matrix.Dense.Internal
import Data.Matrix.Dense.Operations ( apply, applyMat )

infixl 7 <*>, <**>

class Base.Matrix a => IMatrix a e where
    -- | Apply to a vector
    (<*>) :: a (m,n) e -> Vector n e -> Vector m e

    -- | Apply to a matrix
    (<**>) :: a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e

    unsafeApply :: a (m,n) e -> Vector n e -> Vector m e
    
    unsafeApplyMat :: a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e


instance (BLAS3 e) => IMatrix (DMatrix Imm) e where
    (<*>) = apply
    (<**>)  = applyMat
    