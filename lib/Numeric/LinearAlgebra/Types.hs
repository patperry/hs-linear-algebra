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

    -- * Matrix views
    Herm(..),
    withHerm,

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

-- | A hermitian view of an underlying matrix.  The view can either be
-- of the upper or lower triangular part of the matrix.  The type arguments
-- are as follows:
--
--     * @m@: the underlyting matrix type.
--
--     * @e@: the element type of the matrix.
--
data Herm m e = Herm Uplo (m e) deriving (Show)

-- | Apply a function to the unerlying 'Uplo' and matrix.
withHerm :: Herm m e -> (Uplo -> m e -> a) -> a
withHerm (Herm u m) f = f u m
