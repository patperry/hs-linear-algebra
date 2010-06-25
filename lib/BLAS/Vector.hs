-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Vector
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Immutable dense vectors.

module BLAS.Vector (
    -- * The Vector type
    Vector,
    dimVector,

    -- * Vector construction
    vector, 
    listVector,
    constantVector,

    -- * Accessing vectors
    atVector,
    indicesVector,
    elemsVector,
    assocsVector,

    -- * Incremental vector updates
    replaceVector,
    accumVector,

    -- * Derived vectors
    mapVector,
    zipWithVector,

    -- * Vector views
    spliceVector,
    splitVectorAt,

    -- * Vector properties
    sumAbsVector,
    norm2Vector,
    whichMaxAbsVector,
    dotVector,

    ) where

import BLAS.Vector.Base
import BLAS.Vector.STBase
