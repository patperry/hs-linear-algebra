-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.IO
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.IO (
    -- * The mutable banded matrix data type
    BMatrix(..),
    IOBanded,
    
    module BLAS.Matrix.Base,
    module BLAS.Matrix.ReadOnly,
    module BLAS.Tensor.Base,
    module BLAS.Tensor.Dense.ReadOnly,
    module BLAS.Tensor.ReadOnly,
    module BLAS.Tensor.Mutable,
    
    -- * Creating new matrices
    newBanded,
    newBanded_,
    newListsBanded,

    -- * Views
    -- ** Rows and columns
    rowView,
    colView,
    getRow,
    getCol,
    
    -- ** Diagonals
    diag,
    
    -- * Operations
    module Data.Matrix.Banded.Operations,
    
    -- * Converting to and from banded matrices
    -- ** @ForeignPtr@s
    toForeignPtr,
    fromForeignPtr,

    -- * Bandwith properties
    bandwidth,
    numLower,
    numUpper,
    
    -- * Coercing
    coerceBanded,
    
    -- * Unsafe operations
    unsafeNewBanded,
    unsafeWithElemPtr,
    unsafeRowView,
    unsafeColView,
    unsafeGetRow,
    unsafeGetCol,
    unsafeDiag,
    unsafeFreeze,
    unsafeThaw,
    
    ) where

import Data.Matrix.Banded.Internal
import Data.Matrix.Banded.Operations hiding ( apply, applyMat,
    sapply, sapplyMat, getApply, getApplyMat )
    
import BLAS.Matrix.Base hiding ( Matrix )
import BLAS.Matrix.ReadOnly
import BLAS.Tensor.Base
import BLAS.Tensor.Dense.ReadOnly
import BLAS.Tensor.ReadOnly
import BLAS.Tensor.Mutable
