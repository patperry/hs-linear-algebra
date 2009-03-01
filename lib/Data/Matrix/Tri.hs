-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Tri
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Triangular views of matrices.
--

module Data.Matrix.Tri (
    -- * Triangular matrix types
    Tri(..),

    triFromBase,
    triToBase,
    mapTri,

    lower,
    lowerU,
    upper,
    upperU,

    -- * Overloaded interface for solving linear systems
    module Data.Matrix.Class.ISolve,
    module Data.Matrix.Class.MSolve,
    ) where

import Data.Matrix.TriBase
import Data.Matrix.Class.ISolve
import Data.Matrix.Class.MSolve
