-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Read
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Read (
    ReadMatrix,
    
    unsafeIOToM,
    module Data.Matrix.Dense.Class.Base
    ) where

import BLAS.Tensor
import BLAS.Matrix.RowCol

import Data.Matrix.Dense.Class.Base
import Data.Vector.Dense.Class

class (ReadTensor a (Int,Int) e m, BaseMatrix a e, RowColView a x e, ReadVector x e m) => 
    ReadMatrix a x e m | a -> x where
