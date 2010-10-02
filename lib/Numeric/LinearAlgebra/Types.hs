{-# OPTIONS_HADDOCK hide #-}
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

    -- * Matrix factorization views
    Chol(..),

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

    -- * Algorithm-specific parameters
    CovMethod(..),

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

-- | A Cholesky decomposition view of a matrix.
data Chol m e = Chol Uplo (m e) deriving (Show)

-- | The method of scaling the sample covariance matrix.
data CovMethod =
      UnbiasedCov -- ^ This is the default behavior. Corresponds to a
                  -- scaling of @n/(n-1)@ in the unweighed case, and
                  -- @1/(1 - \\sum w_i^2)@ in the weighted case, where @w_i@
                  -- is the normalized weight. Note the unweighted and
                  -- weighted cases agree when @w_i = 1/n@.
                  
    | MLCov       -- ^ Returns the centered second moment matrix without
                  -- scaling the result.
    deriving (Eq, Show)

