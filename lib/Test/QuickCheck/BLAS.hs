-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.BLAS
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Test generators for BLAS types.
--

module Test.QuickCheck.BLAS (
    -- * Generating random objects
    -- ** Elements
    elements,
    realElements,
    
    -- ** Vectors
    dim,
    vector,
    
    -- ** Matrices
    shape,
    bandwidths,
    matrix,
    banded,
    ) where

import Test.QuickCheck.BLASBase
