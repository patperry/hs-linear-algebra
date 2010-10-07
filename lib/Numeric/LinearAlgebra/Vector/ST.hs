-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Vector.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Mutable vectors in the ST monad.

module Numeric.LinearAlgebra.Vector.ST (
    -- * Mutable vectors
    STVector,
    IOVector,
    create,
    
    -- * Read-only vectors
    RVector(..),
    
    -- * Creating new vectors
    new_,
    new,
    
    -- * Copying vectors
    newCopy,
    copyTo,
    swap,

    -- * Reading and writing vector elements
    read,
    write,
    modify,
    indices,
    getElems,
    getElems',
    getAssocs,
    getAssocs',
    setElems,
    setAssocs,    
    clear,

    -- * List-like operations
    mapTo,
    zipWithTo,

    -- * Vector linear algebra
    getSumAbs,
    getNorm2,
    getWhichMaxAbs,
    getDot,
    scaleByM_,
    addWithScaleM_,
    kroneckerTo,
    
    -- * Vector views    
    withSliceView,
    withDropView,
    withTakeView,
    withSplitAtView,
    
    withSTSliceView,
    withSTDropView,
    withSTTakeView,
    withSTSplitAtView,
    
    -- * Vector math operations
    -- ** Num
    addTo,
    subTo,
    mulTo,
    negateTo,
    conjugateTo,
    absTo,
    signumTo,

    -- ** Fractional
    divTo,
    recipTo,        

    -- ** Floating
    sqrtTo,
    expTo,
    logTo,
    powTo,
    sinTo,
    cosTo,
    tanTo,
    asinTo,
    acosTo,
    atanTo,
    sinhTo,
    coshTo,
    tanhTo,
    asinhTo,
    acoshTo,
    atanhTo,

    -- * Conversions between mutable and immutable vectors
    freeze,
    unsafeFreeze,
    
    -- * Unsafe operations
    unsafeCopyTo,
    unsafeRead,
    unsafeWrite,
    unsafeModify,
    unsafeMapTo,
    unsafeZipWithTo,
    unsafeAddWithScaleM_,
    
    ) where

import Prelude()
import Numeric.LinearAlgebra.Vector.Base
import Numeric.LinearAlgebra.Vector.STBase
