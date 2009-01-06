{-# OPTIONS_HADDOCK hide, prune #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Creating
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Creating (
    -- * Creating new vectors
    newVector_,
    newVector,
    newListVector,
    unsafeNewVector,

    ) where

import Foreign

import BLAS.UnsafeIOToM

import Data.Tensor.Class.MTensor( writeElem, unsafeWriteElem )
import Data.Vector.Dense.Class.Internal



