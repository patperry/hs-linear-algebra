-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Mutable banded matrices in the ST monad.
--

module Data.Matrix.Banded.ST (
    -- * The @STBanded@ data type
    STBanded,
    runSTBanded,

    -- * Overloaded mutable banded matrix interface
    module Data.Matrix.Banded.Class,
    ) where

import Data.Matrix.Banded.Base
import Data.Matrix.Banded.Class
