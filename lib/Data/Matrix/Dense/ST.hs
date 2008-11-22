{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses, UndecidableInstances,
        Rank2Types #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.ST (
    -- * The @STMatrix@ data type
    STMatrix,
    runSTMatrix,

    unsafeIOMatrixToSTMatrix,
    unsafeSTMatrixToIOMatrix,

    module Data.Matrix.Dense.Class,
    ) where

import Control.Monad.ST

import Data.Matrix.Dense.Internal( Matrix(..) )
import Data.Matrix.Dense.Class
import Data.Matrix.Dense.Class.Internal( STMatrix, unsafeIOMatrixToSTMatrix,
    unsafeSTMatrixToIOMatrix )

runSTMatrix :: (forall s . ST s (STMatrix s n e)) -> Matrix n e
runSTMatrix x = runST $ x >>= return . M . unsafeSTMatrixToIOMatrix
