{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses, UndecidableInstances,
        Rank2Types #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.ST (
    -- * The @STBanded@ data type
    STBanded,
    runSTBanded,

    unsafeIOBandedToSTBanded,
    unsafeSTBandedToIOBanded,

    module Data.Matrix.Banded.Class,
    ) where

import Control.Monad.ST

import Data.Matrix.Banded.Internal( Banded(..) )
import Data.Matrix.Banded.Class
import Data.Matrix.Banded.Class.Internal( STBanded, unsafeIOBandedToSTBanded,
    unsafeSTBandedToIOBanded )

runSTBanded :: (forall s . ST s (STBanded s mn e)) -> Banded mn e
runSTBanded x = runST $ x >>= return . B . unsafeSTBandedToIOBanded
