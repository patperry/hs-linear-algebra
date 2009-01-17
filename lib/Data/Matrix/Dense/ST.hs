-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Mutable dense matrices in the ST monad.
--

module Data.Matrix.Dense.ST (
    -- * The @STMatrix@ data type
    STMatrix,
    runSTMatrix,

    -- * Overloaded mutable dense matrix interface
    module Data.Matrix.Dense.Class,
    ) where

import Data.Matrix.Dense.Class
import Data.Matrix.Dense.Base
