{-# LANGUAGE TypeFamilies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.HasVectorView
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.HasVectorView (
    HasVectorView(..)
    ) where

class HasVectorView (a :: * -> * -> *) where
    type VectorView a :: * -> * -> *
