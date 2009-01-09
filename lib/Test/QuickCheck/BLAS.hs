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
    -- * Testable element types
    TestElem(..),
    
    -- * Generating random objects
    -- ** Elements
    elements,
    realElements,
    
    -- ** Vectors
    dim,
    vector,
    
    -- ** Dense matrices
    shape,
    matrix,

    -- ** Banded matrices
    bandwidths,
    banded,
    ) where

import Test.QuickCheck.BLASBase
