-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.IO
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Mutable vectors in the IO monad.

module Data.Vector.Dense.IO (
    -- * The IOVector data type
    IOVector,
    vectorViewArray,
    vectorViewArrayWithStride,
    withIOVector,
    
    -- * Overloaded mutable dense vector interface
    module Data.Vector.Dense.Class,
    ) where

import Data.Vector.Dense.IOBase
import Data.Vector.Dense.Class
