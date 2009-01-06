-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Immutable dense vectors.

module Data.Vector.Dense (
    -- * The Vector type
    Vector,

    -- * Overloaded interface for vectors
    BaseVector( dim, conj, coerceVector ),

    -- * Creating new vectors
    vector, 
    listVector,

    -- * Special vectors
    zeroVector,
    constantVector,
    basisVector,

    -- * Overloaded interface for reading vector elements
    module Data.Tensor.Class.ITensor,

    -- * Vector views
    subvector,
    subvectorWithStride,

    -- * Vector properties
    sumAbs,
    norm2,
    whichMaxAbs,
    (<.>),

    ) where

import Data.Vector.Dense.Base
import Data.Tensor.Class.ITensor
