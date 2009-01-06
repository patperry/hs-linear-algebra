-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Mutable vectors in the ST monad.

module Data.Vector.Dense.ST (
    -- * The @STVector@ data type
    runSTVector,
    STVector,
    
    -- * Overloaded mutable dense vector interface
    module Data.Vector.Dense.Class
    ) where

import Data.Vector.Dense.STBase( STVector, runSTVector )
import Data.Vector.Dense.Class
