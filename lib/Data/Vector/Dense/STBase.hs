{-# LANGUAGE Rank2Types, FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Vector.Dense.STBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.STBase (
    -- * The @STVector@ data type
    STVector,
    runSTVector,

    unsafeIOVectorToSTVector,
    unsafeSTVectorToIOVector,
    ) where

import Control.Monad.ST
import Data.Vector.Dense.Class.Internal( STVector, unsafeIOVectorToSTVector,
    unsafeSTVectorToIOVector )
import Data.Vector.Dense.Internal ( Vector(..) )

runSTVector :: (forall s . ST s (STVector s n e)) -> Vector n e
runSTVector x = 
    runST $ x >>= return . V . unsafeSTVectorToIOVector
