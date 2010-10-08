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
    -- * Mutable matrices
    STMatrix,
    IOMatrix,
    create,
    
    -- * Read-only matrices
    RMatrix(..),
    
    -- * Creating new matrices
    new_,
    new,
    
    -- * Copying matrices
    newCopy,
    copyTo,

    -- * Matrix views
    withSlice,
    withTakeRows,
    withDropRows,
    withSplitRowsAt,
    withTakeCols,
    withDropCols,
    withSplitColsAt,
    withSliceM,
    withTakeRowsM,
    withDropRowsM,
    withSplitRowsAtM,
    withTakeColsM,
    withDropColsM,
    withSplitColsAtM,
    
    -- * Matrix rows and columns
    rowTo,
    unsafeRowTo,
    setRow,
    unsafeSetRow,
    withCol,

    withColM,
    withColsM,

    -- * Matrix diagonals
    getDiag,
    setDiag,

    -- * Reading and writing matrix elements
    read,
    write,
    modify,
    getIndices,
    getElems,
    getElems',
    getAssocs,
    getAssocs',
    setElems,
    setAssocs,

    -- * List-like operations
    mapTo,
    zipWithTo,

    -- * Matrix math operations
    shiftDiagByM_,
    shiftDiagByWithScaleM_,
    addTo,
    subTo,
    scaleByM_,
    addWithScaleM_,
    scaleRowsByM_,
    scaleColsByM_,
    negateTo,
    conjugateTo,

    -- * Linear algebra
    transTo,
    conjTransTo,
    rank1UpdateM_,
    
    -- ** Matrix-Vector multiplication
    mulVectorTo,
    mulVectorWithScaleTo,
    addMulVectorWithScalesM_,
    
    -- ** Matrix-Matrix multiplication
    mulMatrixTo,
    mulMatrixWithScaleTo,
    addMulMatrixWithScalesM_,

    -- * Conversions between mutable and immutable matrices
    freeze,
    
    -- * Vector views of matrices
    maybeWithVectorM,
    
    -- * Matrix views of vectors
    withFromVector,
    withFromCol,
    withFromRow,

    withFromVectorM,
    withFromColM,
    withFromRowM,
    
    -- * Unsafe operations
    unsafeCopyTo,
    unsafeAddWithScaleM_,
    
    ) where

import Prelude()
import Numeric.LinearAlgebra.Matrix.STBase
