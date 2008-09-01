-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Operations
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Operations (
    -- * Matrix operations
    axpyMatrix,
    unsafeAxpyMatrix,
    
    module BLAS.Numeric,
    module BLAS.Matrix.Mutable,
    
    ) where

import BLAS.Elem( BLAS1 )
import BLAS.Internal( checkBinaryOp )
import BLAS.Tensor( shape )
import BLAS.Numeric
import BLAS.Matrix.Mutable

import Data.Matrix.Dense.Class.Internal

axpyMatrix :: (ReadMatrix a x e m, WriteMatrix b y e m, BLAS1 e) =>
    e -> a n e -> b n e -> m ()
axpyMatrix alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyMatrix alpha x y
{-# INLINE axpyMatrix #-}
