-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.IO
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.IO (
    -- * The IOVector data type
    IOVector,
    
    -- * Overloaded interface for mutable dense vectors
    module Data.Vector.Dense.Class,
    ) where

import Data.Vector.Dense.IOBase( IOVector )
import Data.Vector.Dense.Class
