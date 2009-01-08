-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.IO
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Mutable banded matrices in the IO monad.
--

module Data.Matrix.Banded.IO (
    -- * The IOBanded data type
    IOBanded,
    withIOBanded,

    -- * Overloaded mutable banded matrix interface    
    module Data.Matrix.Banded.Class,
    ) where

import Data.Matrix.Banded.IOBase
import Data.Matrix.Banded.Class
