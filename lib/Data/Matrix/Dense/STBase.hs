-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.STBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.STBase
    where

import Control.Monad.ST

import Data.Matrix.Class
import Data.Matrix.Dense.Base

newtype STMatrix s np e = STMatrix (IOMatrix np e)

instance HasVectorView (STMatrix s) where
    type VectorView (STMatrix s) = STVector s

