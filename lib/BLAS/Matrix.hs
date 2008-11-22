-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix (
    module BLAS.Matrix.Shaped,
    module BLAS.Matrix.Immutable,
    module BLAS.Matrix.Mutable,
    module BLAS.Matrix.Solve,
    ) where

import BLAS.Matrix.Shaped
import BLAS.Matrix.Immutable
import BLAS.Matrix.Mutable
import BLAS.Matrix.Solve
