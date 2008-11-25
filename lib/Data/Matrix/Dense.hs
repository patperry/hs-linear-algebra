-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense (
    module Data.Matrix.Dense.Internal,
    
    -- * Matrix and vector multiplication
    module BLAS.Matrix.IMatrix,
    
    -- * Converting between mutable and immutable matrices
    UnsafeFreezeMatrix(..),
    UnsafeThawMatrix(..),
    freezeMatrix,
    thawMatrix,
    ) where

import Data.Matrix.Dense.Internal hiding ( M )
import qualified Data.Matrix.Dense.Internal as I
import Data.Matrix.Dense.ST
import Data.Matrix.Dense.IO
import BLAS.Matrix.IMatrix

class UnsafeFreezeMatrix a where
    unsafeFreezeMatrix :: a mn e -> Matrix mn e
instance UnsafeFreezeMatrix IOMatrix where
    unsafeFreezeMatrix = I.M
instance UnsafeFreezeMatrix (STMatrix s) where
    unsafeFreezeMatrix = unsafeFreezeMatrix . unsafeSTMatrixToIOMatrix    
    
class UnsafeThawMatrix a where
    unsafeThawMatrix :: Matrix mn e -> a mn e
instance UnsafeThawMatrix IOMatrix where
    unsafeThawMatrix (I.M a) = a
instance UnsafeThawMatrix (STMatrix s) where
    unsafeThawMatrix = unsafeIOMatrixToSTMatrix . unsafeThawMatrix
    
freezeMatrix :: (ReadMatrix a e m, WriteMatrix b e m, UnsafeFreezeMatrix b) =>
    a mn e -> m (Matrix mn e)
freezeMatrix x = do
    x' <- newCopyMatrix x
    return (unsafeFreezeMatrix x')

thawMatrix :: (WriteMatrix a e m) =>
    Matrix mn e -> m (a mn e)
thawMatrix = newCopyMatrix
