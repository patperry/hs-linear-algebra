-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Special
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Special (
    -- * Special matrices
    newZeroMatrix,
    setZeroMatrix,
    newConstantMatrix,
    setConstantMatrix,
    newIdentityMatrix,
    setIdentityMatrix,
    ) where

import BLAS.Tensor( unsafeWriteElem )
import BLAS.Matrix.Shaped( numRows, numCols )

import Data.Matrix.Dense.Class.Internal


-- | Create a new matrix of the given shape with ones along the diagonal, 
-- and zeros everywhere else.
newIdentityMatrix :: (WriteMatrix a e m) => (Int,Int) -> m (a mn e)
newIdentityMatrix mn = do
    a <- newMatrix_ mn
    setIdentityMatrix a
    return a

-- | Set diagonal elements to one and all other elements to zero.
setIdentityMatrix :: (WriteMatrix a e m) => a mn e -> m ()
setIdentityMatrix a = do
    setZeroMatrix a
    mapM_ (\i -> unsafeWriteElem a (i,i) 1) [0..(mn-1)]
  where
    mn = min (numRows a) (numCols a)
