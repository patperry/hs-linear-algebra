-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface to mutable dense vectors.  For vector types
-- than can be used with this interface, see "Data.Vector.Dense.IO" and
-- "Data.Vector.Dense.ST".

module Data.Vector.Dense.Class (
    -- * Overloaded dense vector type classes
    BaseVector( dim, conj, coerceVector ),
    ReadVector,
    WriteVector,
    
    -- * Creating new vectors
    newVector_,
    newVector,
    newListVector,
    
    -- * Special vectors
    newZeroVector,
    newConstantVector,
    newBasisVector,
    setZeroVector,
    setConstantVector,
    setBasisVector,

    -- * Copying vectors
    newCopyVector,
    newCopyVector',
    copyVector,
    swapVector,

    -- * Vector views
    subvectorView,
    subvectorViewWithStride,

    -- * Overloaded interface for reading and writing vector elements
    module Data.Tensor.Class.MTensor,

    -- * Vector operations
    -- ** Unary
    getConjVector,
    getScaledVector,
    getShiftedVector,
    doConjVector,
    scaleByVector,
    shiftByVector,
    
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

    -- * Vector Properties
    getSumAbs,
    getNorm2,
    getWhichMaxAbs,
    getDot,

    -- * Conversions between mutable and immutable vectors
    freezeVector,
    thawVector,
    unsafeFreezeVector,
    unsafeThawVector,
    
    ) where

import Data.Vector.Dense.Base
import Data.Tensor.Class.MTensor
