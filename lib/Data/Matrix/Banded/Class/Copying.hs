-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Class.Copying
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Class.Copying (
    -- * Copying Banded matrices
    newCopyBanded,
    copyBanded,
    unsafeCopyBanded,
    ) where

import Data.Elem.BLAS.Level1( BLAS1 )
import qualified Data.Elem.BLAS.Level1 as BLAS
import BLAS.UnsafeIOToM

import Control.Monad( zipWithM_ )
import Data.Ix( range )
import Foreign( advancePtr )

import Data.Matrix.Class
import Data.Matrix.Banded.Class.Internal
import Data.Matrix.Banded.Class.Views
import Data.Vector.Dense.Base( unsafeCopyVector )

copyBanded :: (WriteBanded b e m, ReadBanded a e m) =>
    b mn e -> a mn e -> m ()
copyBanded dst src
    | shapeBanded dst /= shapeBanded src =
        error "Shape mismatch in copyBanded."
    | bandwidth dst /= bandwidth src =
        error "Bandwidth mismatch in copyBanded."
    | otherwise =
        unsafeCopyBanded dst src


    