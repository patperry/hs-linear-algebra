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
    slice,
    takeRows,
    dropRows,
    splitRowsAt,
    takeCols,
    dropCols,
    splitColsAt,
    
    -- * Matrix rows and columns
    getRow,
    setRow,
    withColView,

    withSTColView,
    withSTColViews,

    -- * Matrix diagonals
    getDiag,
    setDiag,

    -- * Reading and writing matrix elements
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

    -- * List-like operations
    mapTo,
    zipWithTo,

    -- * Matrix math operations
    shiftTo,
    shiftDiagTo,
    shiftDiagToWithScale,    
    addTo,
    addToWithScales,
    subTo,
    scaleTo,
    scaleRowsTo,
    scaleColsTo,
    negateTo,
    conjugateTo,

    -- * Linear algebra
    transTo,
    conjTransTo,
    rank1UpdateTo,
    
    -- ** Matrix-Vector multiplication
    mulToVector,
    mulToVectorWithScale,
    mulAddToVectorWithScales,
    
    -- ** Matrix-Matrix multiplication
    mulToMatrix,
    mulToMatrixWithScale,
    mulAddToMatrixWithScales,

    -- * Conversions between mutable and immutable matrices
    freeze,
    unsafeFreeze,
    
    -- * Vector views of matrices
    maybeWithSTVectorView,
    
    -- * Matrix views of vectors
    withViewFromVector,
    withViewFromCol,
    withViewFromRow,

    withViewFromSTVector,
    withViewFromSTCol,
    withViewFromSTRow,
    
    -- * Unsafe operations
    unsafeCopyTo,
    
    ) where

import Prelude()
import Numeric.LinearAlgebra.Matrix.Base
import Numeric.LinearAlgebra.Matrix.STBase
