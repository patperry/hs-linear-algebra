-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Immutable dense matrices.

module Numeric.LinearAlgebra.Matrix (
    -- * The Matrix type
    Matrix,
    dimMatrix,

    -- * Matrix construction
    matrix, 
    listMatrix,
    colListMatrix,
    rowListMatrix,
    constantMatrix,

    -- * Accessing Matrixs
    atMatrix,
    indicesMatrix,
    elemsMatrix,
    assocsMatrix,

    -- * Incremental Matrix updates
    replaceMatrix,
    accumMatrix,

    -- * Derived matrices
    mapMatrix,
    zipWithMatrix,

    -- * Matrix views
    colMatrix,
    colsMatrix,
    sliceMatrix,
    splitRowsMatrixAt,
    splitColsMatrixAt,

    -- * Matrix rows
    rowMatrix,
    rowsMatrix,

    -- * Matrix math operations
    shiftMatrix,
    shiftDiagMatrix,
    shiftDiagMatrixWithScale,    
    addMatrix,
    addMatrixWithScales,
    subMatrix,
    scaleMatrix,
    scaleRowsMatrix,
    scaleColsMatrix,
    negateMatrix,
    conjMatrix,

    -- * Linear algebra
    transMatrix,
    conjTransMatrix,
    rank1UpdateMatrix,
    
    -- ** Matrix-Vector multiplication
    mulMatrixVector,
    mulMatrixVectorWithScale,
    mulMatrixAddVector,
    mulMatrixAddVectorWithScales,
    
    -- ** Matrix-Matrix multiplication
    mulMatrixMatrix,
    mulMatrixMatrixWithScale,
    mulMatrixAddMatrix,
    mulMatrixAddMatrixWithScales,

    ) where

import Numeric.LinearAlgebra.Matrix.Base
import Numeric.LinearAlgebra.Matrix.STBase
