-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.IO
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- This modules defines a mutable dense matrix and associate operations.
module Data.Matrix.Dense.IO (
    -- * The mutable dense matrix data type
    DMatrix(..),
    IOMatrix,
    
    module BLAS.Matrix.Base,
    module BLAS.Matrix.ReadOnly,
    module BLAS.Tensor.Base,
    module BLAS.Tensor.Dense.ReadOnly,
    module BLAS.Tensor.ReadOnly,
    module BLAS.Tensor.Mutable,
    
    -- * Creating new matrices
    newMatrix,
    newMatrix_,
    newListMatrix,
    newColsMatrix,
    newRowsMatrix,

    -- * Special matrices
    newIdentity,
    setIdentity,
    
    -- * Views
    -- ** Rows and columns
    row,
    col,
    rows,
    cols,
    
    -- ** Diagonals
    diag,
    
    -- ** Matrix views
    submatrix,

    -- * Operations
    module Data.Matrix.Dense.Operations,
    
    -- ** Lifting scalar and vector operations
    liftV,
    liftV2,

    -- * Converting to and from matrices
    -- ** Vectors
    maybeFromRow,
    maybeFromCol,
    maybeToVector,
    
    -- ** @ForeignPtr@s
    toForeignPtr,
    fromForeignPtr,
    
    -- ** Coercing
    coerceMatrix,
    
    -- * Unsafe operations
    unsafeNewMatrix,
    unsafeWithElemPtr,
    unsafeRow,
    unsafeCol,
    unsafeDiag,
    unsafeSubmatrix,
    unsafeFreeze,
    unsafeThaw,
    
    ) where

import Data.Matrix.Dense.Internal
import Data.Matrix.Dense.Operations hiding ( apply, applyMat,
    sapply, sapplyMat, add, plus, minus, times, divide, getApply, getApplyMat )
    
import BLAS.Matrix.Base hiding ( Matrix )
import BLAS.Matrix.ReadOnly
import BLAS.Tensor.Base
import BLAS.Tensor.Dense.ReadOnly
import BLAS.Tensor.ReadOnly
import BLAS.Tensor.Mutable
