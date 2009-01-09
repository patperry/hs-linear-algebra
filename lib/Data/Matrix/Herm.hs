-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Herm
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Hermitian views of matrices.
--

module Data.Matrix.Herm (
    Herm(..),
    hermFromBase,
    hermToBase,
    mapHerm,
    hermL,
    hermU,
    coerceHerm,
    ) where

import Data.Matrix.HermBase
