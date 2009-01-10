-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface to mutable banded matrices.  For matrix types
-- than can be used with this interface, see "Data.Matrix.Banded.IO" and
-- "Data.Matrix.Banded.ST".  Many of these functions can also be used with
-- the immutable type defined in "Data.Matrix.Banded".
--

module Data.Matrix.Banded.Class (
    -- * Banded matrix type classes
    BaseBanded( numLower, numUpper, bandwidths, ldaBanded, isHermBanded
              , transEnumBanded, maybeMatrixStorageFromBanded
              , maybeBandedFromMatrixStorage, coerceBanded ),
    ReadBanded,
    WriteBanded,

    -- * Overloaded interface for matrices
    module Data.Matrix.Class,
    module Data.Matrix.Class.MMatrix,

    -- * Creating banded matrices
    newBanded_,
    newBanded,
    newListsBanded,

    -- * Special banded matrices
    newZeroBanded,
    setZeroBanded,
    newConstantBanded,
    setConstantBanded,

    -- * Copying banded matrices
    newCopyBanded,
    copyBanded,

    -- * Conversions between banded matrices and vectors
    viewVectorAsBanded,
    viewVectorAsDiagBanded,
    maybeViewBandedAsVector,

    -- * Row and column views
    rowViewBanded,
    colViewBanded,
    diagViewBanded,
    
    -- * Getting diagonals
    getDiagBanded,

    -- * Overloaded interface for reading and writing banded matrix elements
    module Data.Tensor.Class,
    module Data.Tensor.Class.MTensor,

    -- * Conversions between mutable and immutable banded matrices
    freezeBanded,
    thawBanded,
    unsafeFreezeBanded,
    unsafeThawBanded,

    -- * Conversions from @IOBanded@s
    unsafeBandedToIOBanded,
    unsafeConvertIOBanded,
    unsafePerformIOWithBanded,
    
    ) where

import Data.Matrix.Banded.Base

import Data.Matrix.Class
import Data.Matrix.Class.MMatrix

import Data.Tensor.Class
import Data.Tensor.Class.MTensor
