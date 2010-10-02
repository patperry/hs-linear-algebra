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
    -- * Immutable vectors
    Vector,
    
    -- * Read-only vectors
    RVector(..),

    -- * Vector construction
    fromAssocs, 
    fromList,
    constant,

    -- * Accessing vectors
    at,
    unsafeAt,
    indices,
    elems,
    assocs,

    -- * Incremental vector updates
    replace,
    accum,

    -- * Derived vectors
    map,
    zipWith,
    concat,

    -- * Vector views
    slice,
    splitAt,
    drop,
    take,

    -- * Vector properties
    sumAbs,    
    norm2,
    whichMaxAbs,
    dot,
    kronecker,

    -- * Vector math operations
    -- ** Num
    shift,
    add,
    addWithScales,
    sub,
    scale,
    mul,
    negate,
    conjugate,
    abs,
    signum,

    -- ** Fractional
    div,
    recip,        

    -- ** Floating
    sqrt,
    exp,
    log,
    pow,
    sin,
    cos,
    tan,
    asin,
    acos,
    atan,
    sinh,
    cosh,
    tanh,
    asinh,
    acosh,
    atanh,
    
    -- * Mutable interface
    module Numeric.LinearAlgebra.Vector.ST

    ) where

import Prelude()
import Numeric.LinearAlgebra.Vector.Base
import Numeric.LinearAlgebra.Vector.STBase( RVector(..) )
import Numeric.LinearAlgebra.Vector.ST hiding ( dim, slice, indices, drop,
    take, splitAt, RVector(..) )
