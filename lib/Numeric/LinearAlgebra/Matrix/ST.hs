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
    addToMatrixWithScales,
    subToMatrix,
    scaleToMatrix,
    scaleRowsToMatrix,
    scaleColsToMatrix,
    negateToMatrix,
    conjToMatrix,

    -- * Linear algebra
    transToMatrix,
    conjTransToMatrix,
    rank1UpdateToMatrix,
    
    -- ** Matrix-Vector multiplication
    mulMatrixToVector,
    mulMatrixToVectorWithScale,
    mulMatrixAddToVectorWithScales,
    
    -- ** Matrix-Matrix multiplication
    mulMatrixToMatrix,
    mulMatrixToMatrixWithScale,
    mulMatrixAddToMatrixWithScales,

    -- * Conversions between mutable and immutable matrices
    freezeMatrix,
    unsafeFreezeMatrix,
    
    -- * Matrix views of arrays and vectors
    matrixViewArray,
    matrixViewVector,

    ) where

import Numeric.LinearAlgebra.Matrix.Base
import Numeric.LinearAlgebra.Matrix.STBase
