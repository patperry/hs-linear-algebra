-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class.MSolve
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface for solving matrix systems in a monad.  The
-- matrices can operate via inverse multiplication on immutable dense
-- vectors and matrices.
--

module Data.Matrix.Class.MSolve (
    -- * The MSolve type class
    MSolve,

    -- * Solving linear systems
    getSolveVector,
    getSolveMatrix,
    getSSolveVector,
    getSSolveMatrix,
    
    -- * In-place operations
    doSolveVector,
    doSolveMatrix,
    doSSolveVector,
    doSSolveMatrix,
    doSolveVector_,
    doSolveMatrix_,
    doSSolveVector_,
    doSSolveMatrix_,

    ) where

import Data.Matrix.Class.MSolveBase
