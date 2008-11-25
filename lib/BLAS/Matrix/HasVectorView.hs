{-# LANGUAGE TypeFamilies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.HasVectorView
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.HasVectorView (
    HasVectorView(..)
    ) where

class HasVectorView (a :: * -> * -> *) where
    type VectorView a :: * -> * -> *
