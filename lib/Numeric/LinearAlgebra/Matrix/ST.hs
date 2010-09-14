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
    IOMatrix,
    runMatrix,
    
    -- * Read-only Dense matrix type class
    RMatrix( dimMatrix, withColsMatrix, maybeWithVectorViewMatrix, 
        unsafeWithMatrix, unsafeMatrixToForeignPtr,
        unsafeMatrixFromForeignPtr ),
    
    -- * Creating new matrices
    newMatrix_,
    newMatrix,
    
    -- * Copying matrices
    newCopyMatrix,
    copyToMatrix,

    -- * Matrix views
    sliceMatrix,
    takeRowsMatrix,
    dropRowsMatrix,
    splitRowsMatrixAt,
    takeColsMatrix,
    dropColsMatrix,
    splitColsMatrixAt,
    
    -- * Matrix rows and columns
    getRowMatrix,
    setRowMatrix,
    withColMatrix,

    withColMatrixST,    
    withColsMatrixST,

    -- * Matrix diagonals
    getDiagMatrix,
    setDiagMatrix,

    -- * Reading and writing matrix elements
    readMatrix,
    writeMatrix,
    modifyMatrix,
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
    
    -- * Matrix views of vectors
    withMatrixViewVector,
    withMatrixViewColVector,
    withMatrixViewRowVector,

    matrixViewVectorST,
    matrixViewColVectorST,
    matrixViewRowVectorST,
    
    ) where

import Numeric.LinearAlgebra.Matrix.Base
import Numeric.LinearAlgebra.Matrix.STBase
