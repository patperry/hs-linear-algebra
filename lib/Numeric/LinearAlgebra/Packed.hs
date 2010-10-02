-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Packed
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Packed matrices.
--
module Numeric.LinearAlgebra.Packed (
    -- * Immutable packed matrices
    Packed,
    
    -- * Read-only packed matrices
    RPacked(..),
    
    -- * Conversions between vectors and packed matrices
    fromVector,
    toVector,

    -- * Mutable interface
    module Numeric.LinearAlgebra.Packed.ST,

    -- * Hermitian views
    module Numeric.LinearAlgebra.Packed.Herm,

    ) where

import Numeric.LinearAlgebra.Packed.Base
import Numeric.LinearAlgebra.Packed.ST hiding ( RPacked(..) )
import Numeric.LinearAlgebra.Packed.Herm
