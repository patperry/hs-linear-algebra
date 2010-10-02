-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Packed.ST
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Mutable packed matrices.
--
module Numeric.LinearAlgebra.Packed.ST (
    -- * Mutable packed matrices
    STPacked,
    IOPacked,
    create,

    -- * Read-only packed matrices
    RPacked(..),

    -- * Conversions between mutable and immutable packed matrices
    freeze,
    unsafeFreeze,

    -- * Creating new packed matrices
    new_,
    
    -- * Copying matrices
    newCopy,

    -- * Vector views of packed matrices
    withSTVectorView,
    
    -- * Packed matrix views of vectors
    withViewFromVector,
    withViewFromSTVector,
    
    ) where

import Numeric.LinearAlgebra.Packed.Base

