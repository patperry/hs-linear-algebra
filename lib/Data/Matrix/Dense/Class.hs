-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class (
    -- * The dense matrix type classes
    BaseMatrix_(..),
    BaseMatrix,
    ReadMatrix,
    WriteMatrix,
    
    -- * Matrix shape
    module Data.Tensor.Class,
    module Data.Matrix.Class,
    coerceMatrix,

    module Data.Matrix.Dense.Class.Creating,
    module Data.Matrix.Dense.Class.Elements,
    module Data.Matrix.Dense.Class.Special,
    module Data.Matrix.Dense.Class.Views,
    module Data.Matrix.Dense.Class.Copying,
    module Data.Matrix.Dense.Class.Operations,
    
    -- * Low-level functions
    ldaOfMatrix,
    isHermMatrix,
    withMatrixPtr,
    
    ) where

import Data.Matrix.Dense.Class.Internal( BaseMatrix_(..), ldaOfMatrix, 
    isHermMatrix, BaseMatrix, ReadMatrix, WriteMatrix, coerceMatrix, withMatrixPtr )
import Data.Tensor.Class
import Data.Matrix.Class
import Data.Matrix.Dense.Class.Creating
import Data.Matrix.Dense.Class.Elements
import Data.Matrix.Dense.Class.Special
import Data.Matrix.Dense.Class.Views
import Data.Matrix.Dense.Class.Copying
import Data.Matrix.Dense.Class.Operations
