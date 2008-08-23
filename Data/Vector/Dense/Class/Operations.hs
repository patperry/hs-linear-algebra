-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Operations
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Operations (
    -- * Vector operations
    axpyVector,
    unsafeAxpyVector,
    
    module BLAS.Numeric,
    ) where

import BLAS.Elem( BLAS1 )
import BLAS.Internal( checkBinaryOp )
import BLAS.Numeric
import BLAS.Tensor( shape )
import Data.Vector.Dense.Class.Internal

axpyVector :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    e -> x n e -> y n e -> m ()
axpyVector alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyVector alpha x y
{-# INLINE axpyVector #-}
