-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.ST (
    -- * The @STVector@ data type
    STVector,
    runSTVector,
    
    module Data.Vector.Dense.Class
    ) where

import Data.Vector.Dense.STBase( STVector, runSTVector )
import Data.Vector.Dense.Class
