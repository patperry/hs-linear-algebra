{-# LANGUAGE Rank2Types, FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Vector.Dense.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.ST (
    -- * The @STVector@ data type
    STVector,
    runSTVector,

    unsafeIOVectorToSTVector,
    unsafeSTVectorToIOVector,
    
    module Data.Vector.Dense.Class
    ) where

import Control.Monad.ST
import Data.Vector.Dense.Class.Internal( STVector, unsafeIOVectorToSTVector,
    unsafeSTVectorToIOVector )
import Data.Vector.Dense.Class
import Data.Vector.Dense.Internal ( Vector(..) )

runSTVector :: (forall s . ST s (STVector s n e)) -> Vector n e
runSTVector x = 
    runST $ x >>= return . V . unsafeSTVectorToIOVector
