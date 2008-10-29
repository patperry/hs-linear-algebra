-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Class
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Class (
    -- * The banded matrix type classes
    BaseBanded(..),
    ReadBanded,
    WriteBanded,
    
    -- * Matrix shape
    module BLAS.Tensor.Base,
    module BLAS.Matrix.Base,
    coerceBanded,

    module Data.Matrix.Banded.Class.Creating,
    module Data.Matrix.Banded.Class.Elements,
    module Data.Matrix.Banded.Class.Special,
    module Data.Matrix.Banded.Class.Views,
    module Data.Matrix.Banded.Class.Copying,
    
    -- * Low-level functions
    ldaOfBanded,
    isHermBanded,
    withBandedPtr,
    
    ) where

import Data.Matrix.Banded.Class.Internal( BaseBanded(..), ldaOfBanded,
    isHermBanded, ReadBanded, WriteBanded, coerceBanded, withBandedPtr )
import BLAS.Tensor.Base
import BLAS.Matrix.Base hiding ( BaseMatrix )
import Data.Matrix.Banded.Class.Creating
import Data.Matrix.Banded.Class.Elements
import Data.Matrix.Banded.Class.Special
import Data.Matrix.Banded.Class.Views
import Data.Matrix.Banded.Class.Copying
