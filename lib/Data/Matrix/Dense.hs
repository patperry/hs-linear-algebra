{-# LANGUAGE FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Immutable dense matrices.
--

module Data.Matrix.Dense (
    -- * Dense matrix type
    Matrix,

    -- * Overloaded interface for dense matrices
    BaseMatrix( isHermMatrix, coerceMatrix ),

    -- * Overloaded interface for matrices
    module Data.Matrix.Class,
    module Data.Matrix.Class.IMatrix,    
    
    -- * Creating matrices
    matrix, 
    listMatrix,
    rowsMatrix,
    colsMatrix,
    
    -- * Special matrices
    zeroMatrix,
    constantMatrix,
    identityMatrix,

    -- * Conversions between vectors and matrices
    matrixFromRow,
    matrixFromCol,
    matrixFromVector,
    vectorFromMatrix,

    -- * Matrix views
    submatrix,
    splitRowsAt,
    splitColsAt,
    
    -- * Vector views
    diag,

    -- * Overloaded interface for reading matrix elements
    module Data.Tensor.Class,
    module Data.Tensor.Class.ITensor,

    ) where

import Data.Matrix.Dense.Base

import Data.Matrix.Class
import Data.Matrix.Class.IMatrix
    
import Data.Tensor.Class
import Data.Tensor.Class.ITensor
