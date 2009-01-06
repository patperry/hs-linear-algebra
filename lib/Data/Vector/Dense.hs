-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense (
    -- * The Vector type
    Vector,

    -- * Overloaded interface for vectors
    BaseVector( dim ),

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
    conj,
    subvector,
    subvectorWithStride,

    -- * Vector properties
    sumAbs,
    norm2,
    whichMaxAbs,
    (<.>),

    -- * Coercing the phantom shape type
    coerceVector,
    ) where

import Data.Vector.Dense.Base
import Data.Tensor.Class.ITensor
