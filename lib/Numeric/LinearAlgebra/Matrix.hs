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
    dim,

    -- * Matrix construction
    fromAssocs, 
    fromList,
    fromCols,
    fromRows,
    constant,

    -- * Accessing matrices
    at,
    indices,
    elems,
    assocs,

    -- * Incremental matrix updates
    replace,
    accum,

    -- * Derived matrices
    map,
    zipWith,

    -- * Matrix views
    slice,
    takeRows,
    dropRows,
    splitRowsAt,
    takeCols,
    dropCols,
    splitColsAt,

    -- * Matrix rows and columns
    col,
    cols,
    row,
    rows,
    
    -- * Matrix diagonals
    diag,
    
    -- * Conversions to vectors
    toVector,
    
    -- * Conversions from vectors
    fromVector,
    fromCol,
    fromRow,

    -- * Matrix math operations
    shift,
    shiftDiag,
    shiftDiagWithScale,    
    add,
    addWithScales,
    sub,
    scale,
    scaleRows,
    scaleCols,
    negate,
    conjugate,

    -- * Linear algebra
    trans,
    conjTrans,
    rank1Update,
    
    -- ** Matrix-Vector multiplication
    mulVector,
    mulVectorWithScale,
    mulAddVectorWithScales,
    
    -- ** Matrix-Matrix multiplication
    mulMatrix,
    mulMatrixWithScale,
    mulAddMatrixWithScales,

    ) where

import Prelude()
import Numeric.LinearAlgebra.Matrix.Base
import Numeric.LinearAlgebra.Matrix.STBase
