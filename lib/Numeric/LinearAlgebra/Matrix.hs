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
    -- * Immutable matrices
    Matrix,

    -- * Read-only matrices
    RMatrix(..),

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
    shiftDiag,
    shiftDiagWithScale,    
    add,
    addWithScale,
    sub,
    scaleBy,
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

    -- * Mutable interface
    module Numeric.LinearAlgebra.Matrix.ST,
    
    -- * Hermitian views
    module Numeric.LinearAlgebra.Matrix.Herm,
    
    -- * Cholesky factorizations
    module Numeric.LinearAlgebra.Matrix.Cholesky,
    
    -- * Eigenvalues and eigenvectors
    module Numeric.LinearAlgebra.Matrix.Eigen,

    -- * Basic multivariate statistics
    module Numeric.LinearAlgebra.Matrix.Statistics,

    ) where

import Prelude()
import Numeric.LinearAlgebra.Matrix.Base
import Numeric.LinearAlgebra.Matrix.Herm
import Numeric.LinearAlgebra.Matrix.STBase( indices, slice, takeRows,
    dropRows, splitRowsAt, takeCols, dropCols, splitColsAt, RMatrix(..) )
import Numeric.LinearAlgebra.Matrix.ST hiding ( indices, slice, takeRows,
    dropRows, splitRowsAt, takeCols, dropCols, splitColsAt, RMatrix(..) )
import Numeric.LinearAlgebra.Matrix.Cholesky
import Numeric.LinearAlgebra.Matrix.Eigen
import Numeric.LinearAlgebra.Matrix.Statistics
