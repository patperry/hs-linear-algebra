-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.IO
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Mutable dense matrices in the IO monad.
--

module Data.Matrix.Dense.IO (
    -- * The IOMatrix data type
    IOMatrix,
    
    -- * Overloaded mutable dense matrix interface
    module Data.Matrix.Dense.Class,
    ) where

import Data.Matrix.Dense.IOBase( IOMatrix )
import Data.Matrix.Dense.Class
