-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Mutable matrices in the ST monad.

module Numeric.LinearAlgebra.Matrix.ST (
    -- * The @STMatrix@ data type
    STMatrix,
    runMatrix,
    
    -- * Read-only Dense matrix type class
    RMatrix( dimMatrix, maybeVectorViewMatrix, unsafeWithMatrix ),
    
    -- * Creating new matrices
    newMatrix_,
    newMatrix,
    
    -- * Copying matrices
    newCopyMatrix,
    copyToMatrix,

    -- * Matrix views
    colMatrix,
    colsMatrix,
    sliceMatrix,
    splitRowsMatrixAt,
    splitColsMatrixAt,
    
    -- * Matrix rows
    getRowMatrix,
    setRowMatrix,   

    -- * Reading and writing matrix elements
    readMatrix,
    writeMatrix,
    updateMatrix,
    indicesMatrix,
    getElemsMatrix,
    getElemsMatrix',
    getAssocsMatrix,
    getAssocsMatrix',
    setElemsMatrix,
    setAssocsMatrix,

    -- * List-like operations
    mapToMatrix,
    zipWithToMatrix,

    -- * Matrix math operations
    shiftToMatrix,
    shiftDiagToMatrix,
    shiftDiagToMatrixWithScale,    
    addToMatrix,
    addToMatrixWithScale,
    subToMatrix,
    scaleToMatrix,
    scaleRowsToMatrix,
    scaleColsToMatrix,
    negateToMatrix,

    -- * Linear algebra
    rank1UpdateToMatrix,

    -- * Conversions between mutable and immutable matrices
    freezeMatrix,
    unsafeFreezeMatrix,
    thawMatrix,
    unsafeThawMatrix,
    
    -- * Matrix views of arrays and vectors
    matrixViewArray,
    matrixViewVector,

    ) where

import Numeric.LinearAlgebra.Matrix.Base
import Numeric.LinearAlgebra.Matrix.STBase
