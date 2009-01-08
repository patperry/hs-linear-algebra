-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class.IMatrix
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface for immutable matrices.  The matrices provide
-- access to rows and columns, and can operate via multiplication on 
-- immutable dense vectors and matrices.
--

module Data.Matrix.Class.IMatrix (
    -- * The IMatrix type class
    IMatrix,

    -- * Rows and columns
    row,
    col,
    rows,
    cols,

    -- * Multiplication
    (<*>),
    (<**>),
    sapply,
    sapplyMat,
    ) where

import Data.Matrix.Class.IMatrixBase
