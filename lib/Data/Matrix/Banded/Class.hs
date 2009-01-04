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
    BaseBanded_(..),
    BaseBanded,
    ReadBanded,
    WriteBanded,
        
    -- * Banded matrix shape
    module BLAS.Tensor.Base,
    module Data.Matrix.Class,
    coerceBanded,
    
    -- * Bandwidth
    numLower,
    numUpper,
    bandwidth,

    module Data.Matrix.Banded.Class.Creating,
    module Data.Matrix.Banded.Class.Elements,
    module Data.Matrix.Banded.Class.Special,
    module Data.Matrix.Banded.Class.Views,
    module Data.Matrix.Banded.Class.Copying,
    
    -- * Low-level functions
    ldaOfBanded,
    isHermBanded,
    withBandedPtr,
    withBandedElemPtr,
    
    ) where

import Data.Matrix.Banded.Class.Internal( BaseBanded_(..), 
    BaseBanded, ReadBanded, 
    WriteBanded, numLower, numUpper, bandwidth, ldaOfBanded, isHermBanded,
    coerceBanded, withBandedPtr, withBandedElemPtr )
import BLAS.Tensor.Base
import Data.Matrix.Class
import Data.Matrix.Banded.Class.Creating
import Data.Matrix.Banded.Class.Elements
import Data.Matrix.Banded.Class.Special
import Data.Matrix.Banded.Class.Views
import Data.Matrix.Banded.Class.Copying
