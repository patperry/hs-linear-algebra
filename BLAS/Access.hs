{-# LANGUAGE EmptyDataDecls #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Access
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

-- | Types for access control.
module BLAS.Access (
    Mut,
    Imm
    ) where

-- | Tag for mutable types.
data Mut

-- | Tag for immutable types.
data Imm
