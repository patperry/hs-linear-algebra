-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Immutable dense matrices.

module BLAS.Matrix (
    -- * The Matrix type
    Matrix,
    dimMatrix,

    -- * Matrix construction
    matrix, 
    listMatrix,
    constantMatrix,

    -- * Accessing Matrixs
    atMatrix,
    indicesMatrix,
    elemsMatrix,
    assocsMatrix,

    -- * Incremental Matrix updates
    replaceMatrix,
    accumMatrix,

    -- * Matrix views
    colMatrix,
    colsMatrix,
    spliceMatrix,
    splitRowsMatrixAt,
    splitColsMatrixAt,

    -- * Matrix math operations
    shiftMatrix,
    shiftDiagMatrix,
    shiftDiagMatrixWithScale,    
    addMatrix,
    addMatrixWithScale,
    subMatrix,
    scaleMatrix,
    scaleRowsMatrix,
    scaleColsMatrix,
    negateMatrix,

    ) where

import BLAS.Matrix.Base
import BLAS.Matrix.STBase
