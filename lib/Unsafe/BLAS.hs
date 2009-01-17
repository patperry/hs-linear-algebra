-----------------------------------------------------------------------------
-- |
-- Module     : Unsafe.BLAS
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Unsafe BLAS functions.  Most of these functions to not check the
-- shapes of their arguments, and should be used with caution.
--

module Unsafe.BLAS (

    -- * Utility functions
    clearArray,
    inlinePerformIO,

    -- * Vector functions
    IOVector(..),
    unsafeIOVectorToVector,
    unsafePerformIOWithVector,
    unsafeVector,
    unsafeSubvector,
    unsafeSubvectorWithStride,
    unsafeDot,
    
    unsafeNewVector,
    unsafeCopyVector,
    unsafeSwapVector,
    unsafeSubvectorView,
    unsafeSubvectorViewWithStride,
    unsafeGetAddVector,
    unsafeGetSubVector,
    unsafeGetMulVector,
    unsafeGetDivVector,
    unsafeAxpyVector,
    unsafeAddVector,
    unsafeSubVector,
    unsafeMulVector,
    unsafeDivVector,
    unsafeGetDot,

    -- * Matrix functions
    IOMatrix(..),
    unsafeIOMatrixToMatrix,
    unsafePerformIOWithMatrix,
    
    unsafeSubmatrixView,
    unsafeDiagView,
    unsafeRowView,
    unsafeColView,

    unsafeMatrix,
    unsafeSubmatrix,
    unsafeDiag,
    
    unsafeNewMatrix,
    unsafeCopyMatrix,
    unsafeSwapMatrix,
    unsafeSwapRows,
    unsafeSwapCols,
    unsafeGetDiag,
    unsafeGetAddMatrix,
    unsafeGetSubMatrix,
    unsafeGetMulMatrix,
    unsafeGetDivMatrix,
    unsafeAxpyMatrix,
    unsafeAddMatrix,
    unsafeSubMatrix,
    unsafeMulMatrix,
    unsafeDivMatrix,
    unsafeRank1UpdateMatrix,
    
    -- * Banded functions
    IOBanded(..),
    unsafeIOBandedToBanded,
    unsafePerformIOWithBanded,
    
    unsafeDiagViewBanded,
    unsafeRowViewBanded,
    unsafeColViewBanded,
    
    unsafeBanded,
    
    unsafeNewBanded,
    unsafeGetDiagBanded,
    unsafeGetRowBanded,
    unsafeGetColBanded,
    
    unsafeDiagBanded,
    
    -- * Matrix type classes
    IMatrix(..),
    ISolve(..),
    MMatrix(..),
    MSolve(..),
    unsafeGetSSolveVector,
    unsafeGetSSolveMatrix,

    ) where

import BLAS.Internal

import Data.Vector.Dense.Base
import Data.Vector.Dense.IOBase
import Data.Matrix.Dense.Base
import Data.Matrix.Dense.IOBase
import Data.Matrix.Banded.Base
import Data.Matrix.Banded.IOBase

import Data.Matrix.Class.ISolveBase
import Data.Matrix.Class.IMatrixBase
import Data.Matrix.Class.MMatrixBase
import Data.Matrix.Class.MSolveBase

