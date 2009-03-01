-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface to mutable dense matrices.  For matrix types
-- than can be used with this interface, see "Data.Matrix.Dense.IO" and
-- "Data.Matrix.Dense.ST".  Many of these functions can also be used with
-- the immutable type defined in "Data.Matrix.Dense".
--

module Data.Matrix.Dense.Class (
    -- * Dense matrix type classes
    BaseMatrix( ldaMatrix, isHermMatrix, transEnumMatrix ),
    ReadMatrix,
    WriteMatrix,
    
    -- * Overloaded interface for matrices
    module Data.Matrix.Class,
    module Data.Matrix.Class.MMatrix,
    
    -- * Creating matrices
    newMatrix_,
    newMatrix,
    newListMatrix,
    newRowsMatrix,
    newColsMatrix,
    newRowMatrix,
    newColMatrix,

    -- * Special matrices
    newZeroMatrix,
    setZeroMatrix,
    newConstantMatrix,
    setConstantMatrix,
    newIdentityMatrix,
    setIdentityMatrix,
    
    -- * Copying matrices
    newCopyMatrix,
    newCopyMatrix',
    copyMatrix,
    swapMatrix,
    
    -- * Swapping rows and columns
    swapRows,
    swapCols,

    -- * Matrix views
    submatrixView,
    splitRowsAt,
    splitColsAt,

    -- * Row and column views
    rowViews,
    colViews,
    rowView,
    colView,
    diagView,

    -- * Conversions between matrices and vectors
    liftMatrix,
    liftMatrix2,
    maybeViewMatrixAsVector,
    maybeViewVectorAsMatrix,
    maybeViewVectorAsRow,
    maybeViewVectorAsCol,  
    
    -- * Getting diagonals
    getDiag,
    
    -- * Overloaded interface for reading and writing matrix elements
    module Data.Tensor.Class,
    module Data.Tensor.Class.MTensor,

    -- * Matrix operations
    -- ** Unary
    getConjMatrix,
    getScaledMatrix,
    getShiftedMatrix,
    doConjMatrix,
    scaleByMatrix,
    shiftByMatrix,
    
    -- ** Binary
    getAddMatrix,
    getSubMatrix,
    getMulMatrix,
    getDivMatrix,
    addMatrix,
    subMatrix,
    axpyMatrix,
    mulMatrix,
    divMatrix,
    
    -- ** Low rank updates
    rank1UpdateMatrix,

    -- * Conversions between mutable and immutable matrices
    freezeMatrix,
    thawMatrix,
    unsafeFreezeMatrix,
    unsafeThawMatrix,

    -- * Conversions from @IOMatrix@s
    unsafeMatrixToIOMatrix,
    unsafeConvertIOMatrix,
    unsafePerformIOWithMatrix,
        
    ) where

import Data.Matrix.Dense.Base

import Data.Matrix.Class
import Data.Matrix.Class.MMatrix

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

