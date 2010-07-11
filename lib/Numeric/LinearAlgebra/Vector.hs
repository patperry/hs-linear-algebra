-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Vector
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Immutable dense vectors.

module Numeric.LinearAlgebra.Vector (
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
    sumVector,
    sumAbsVector,    
    norm2Vector,
    whichMaxAbsVector,
    dotVector,

    -- * Vector math operations
    -- ** Num
    shiftVector,
    addVector,
    addVectorWithScale,
    subVector,
    scaleVector,
    mulVector,
    negateVector,
    conjVector,
    absVector,
    signumVector,

    -- ** Fractional
    divVector,
    recipVector,        

    -- ** Floating
    sqrtVector,
    expVector,
    logVector,
    powVector,
    sinVector,
    cosVector,
    tanVector,
    asinVector,
    acosVector,
    atanVector,
    sinhVector,
    coshVector,
    tanhVector,
    asinhVector,
    acoshVector,
    atanhVector,

    ) where

import Numeric.LinearAlgebra.Vector.Base
import Numeric.LinearAlgebra.Vector.STBase
