-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Operations
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Operations (
    -- * Vector operations
    -- ** Unary
    getConjVector,
    getScaledVector,
    getShiftedVector,
    
    -- ** Binary
    getAddVector,
    getSubVector,
    getMulVector,
    getDivVector,
    addVector,
    subVector,
    axpyVector,
    mulVector,
    divVector,

    -- ** Unsafe
    unsafeGetAddVector,
    unsafeGetSubVector,
    unsafeGetMulVector,
    unsafeGetDivVector,
    unsafeAddVector,
    unsafeSubVector,
    unsafeAxpyVector,
    unsafeMulVector,
    unsafeDivVector,
    ) where

import BLAS.Internal( checkBinaryOp )
import Data.Tensor.Class( Shaped(..) )
import Data.Vector.Dense.Class.Internal

---------------------------- Unary Operations -----------------------------



---------------------------- Binary Operations -----------------------------


