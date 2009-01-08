-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Class.ISolve
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- An overloaded interface for solving immutable matrix systems.  The
-- matrices can operate via inverse multiplication on immutable dense
-- vectors and matrices.
--

module Data.Matrix.Class.ISolve (
    -- * The ISolve type class
    ISolve,
    
    -- * Solving linear systems
    (<\>),
    (<\\>),
    ssolve,
    ssolveMat,
    
    ) where

import Data.Matrix.Class.ISolveBase
