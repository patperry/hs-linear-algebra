{-# LANGUAGE TypeFamilies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Types
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Basic type classes and enums.
--

module Numeric.LinearAlgebra.Types (

    -- * Vector math types
    VNum,
    VFractional,
    VFloating,

    -- * BLAS element types
    BLAS1,
    BLAS2,
    BLAS3,
    
    -- * LAPACK element types
    LAPACK,
    
    -- * Matrix types
    HasVectorView(..),

    -- * Enums
    Trans(..),
    Uplo(..),
    Side(..),
    Diag(..),

    -- * Re-export of Complex from Data.Complex
    module Data.Complex,

    -- * Re-export of Storable from Foreign.Storable
    module Foreign.Storable,
        
    ) where

import Foreign.VMath( VNum, VFractional, VFloating )
import Foreign.BLAS( Trans(..), Uplo(..), Side(..), Diag(..), BLAS1, BLAS2, BLAS3 )
import Foreign.LAPACK( LAPACK )
import Data.Complex( Complex(..) )
import Foreign.Storable( Storable() )

-- | Types that can be viewed as vectors.
class HasVectorView (a :: * -> *) where
    type VectorView a :: * -> *

