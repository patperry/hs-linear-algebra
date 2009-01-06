-----------------------------------------------------------------------------
-- |
-- Module     : Data.Elem.BLAS
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Type classes for elements with BLAS support.
--

module Data.Elem.BLAS (
    -- * Basic element type classes
    module Data.Elem.BLAS.Base,
    -- * Low-level BLAS functions
    module Data.Elem.BLAS.Level1,
    module Data.Elem.BLAS.Level2,
    module Data.Elem.BLAS.Level3,
    module BLAS.Types,
    -- * Re-export of Complex from Data.Complex
    module Data.Complex
    ) where

import Data.Elem.BLAS.Base
import Data.Elem.BLAS.Level1
import Data.Elem.BLAS.Level2
import Data.Elem.BLAS.Level3
import BLAS.Types
import Data.Complex( Complex(..) )
