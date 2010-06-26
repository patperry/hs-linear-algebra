-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Vector.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Mutable vectors in the ST monad.

module BLAS.Vector.ST (
    -- * The @STVector@ data type
    STVector,
    runVector,
    
    -- * Read-only Dense vector type class
    RVector( dimVector, unsafeWithVector ),
    
    -- * Creating new vectors
    newVector_,
    newVector,
    
    -- * Copying vectors
    newCopyVector,
    copyToVector,
    swapVector,

    -- * Vector views
    spliceVector,
    splitVectorAt,

    -- * Reading and writing vector elements
    readVector,
    writeVector,
    getElemsVector,
    getElemsVector',
    getAssocsVector,
    getAssocsVector',
    setElemsVector,
    setAssocsVector,    

    -- * List-like operations
    mapToVector,
    zipWithToVector,

    -- * Vector properties
    getSumVector,    
    getSumAbsVector,
    getNorm2Vector,
    getWhichMaxAbsVector,
    getDotVector,

    -- * Vector math operations
    -- ** Num
    shiftToVector,
    addToVector,
    addToVectorWithScale,
    subToVector,
    scaleToVector,
    mulToVector,
    negateToVector,
    conjToVector,
    absToVector,
    signumToVector,

    -- ** Fractional
    divToVector,
    recipToVector,        

    -- ** Floating
    sqrtToVector,
    expToVector,
    logToVector,
    powToVector,
    sinToVector,
    cosToVector,
    tanToVector,
    asinToVector,
    acosToVector,
    atanToVector,
    sinhToVector,
    coshToVector,
    tanhToVector,
    asinhToVector,
    acoshToVector,
    atanhToVector,

    -- * Conversions between mutable and immutable vectors
    freezeVector,
    unsafeFreezeVector,
    thawVector,
    unsafeThawVector,
    
    -- * Vector views of arrays
    vectorViewArray,

    ) where

import BLAS.Vector.Base
import BLAS.Vector.STBase
