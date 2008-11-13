{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded (
    module Data.Matrix.Banded.Internal,
    
    -- * Converting between mutable and immutable banded matrices
    UnsafeFreezeBanded(..),
    UnsafeThawBanded(..),
    freezeBanded,
    thawBanded,
    ) where

import BLAS.Elem
import Data.Matrix.Banded.Internal hiding ( B )
import qualified Data.Matrix.Banded.Internal as I
import Data.Matrix.Banded.ST
import Data.Matrix.Banded.IO

class UnsafeFreezeBanded a where
    unsafeFreezeBanded :: a mn e -> Banded mn e
instance UnsafeFreezeBanded IOBanded where
    unsafeFreezeBanded = I.B
instance UnsafeFreezeBanded (STBanded s) where
    unsafeFreezeBanded = unsafeFreezeBanded . unsafeSTBandedToIOBanded    
    
class UnsafeThawBanded a where
    unsafeThawBanded :: Banded mn e -> a mn e
instance UnsafeThawBanded IOBanded where
    unsafeThawBanded (I.B a) = a
instance UnsafeThawBanded (STBanded s) where
    unsafeThawBanded = unsafeIOBandedToSTBanded . unsafeThawBanded
    
freezeBanded :: (ReadBanded a m, WriteBanded b m, UnsafeFreezeBanded b, BLAS1 e) =>
    a mn e -> m (Banded mn e)
freezeBanded x = do
    x' <- newCopyBanded x
    return (unsafeFreezeBanded x')

thawBanded :: (WriteBanded a m, BLAS1 e) =>
    Banded mn e -> m (a mn e)
thawBanded = newCopyBanded
