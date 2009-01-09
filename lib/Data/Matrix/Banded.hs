-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Immutable banded matrices.
--

module Data.Matrix.Banded (
    -- * Banded matrix type
    Banded,

    -- * Overloaded interface for banded matrices
    BaseBanded( numLower, numUpper, bandwidths
              , maybeMatrixStorageFromBanded, maybeBandedFromMatrixStorage, coerceBanded ),

    -- * Overloaded interface for matrices
    module Data.Matrix.Class,
    module Data.Matrix.Class.IMatrix,    
    
    -- * Creating banded matrices
    banded,
    listsBanded,

    -- * Special banded matrices
    zeroBanded,
    constantBanded,

    -- * Conversions between vectors and banded matrices
    bandedFromVector,
    diagBandedFromVector,
    maybeVectorFromBanded,

    -- * Vector views
    diagBanded,

    -- * Overloaded interface for reading banded matrix elements
    module Data.Tensor.Class,
    module Data.Tensor.Class.ITensor,
    ) where

import Data.Matrix.Banded.Base

import Data.Matrix.Class
import Data.Matrix.Class.IMatrix

import Data.Tensor.Class
import Data.Tensor.Class.ITensor
