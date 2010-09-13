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
    
    -- * Enums
    Trans(..),
    Uplo(..),
    Side(..),
    Diag(..),

    -- ** Eigenvalue computations
    EigJob(..),
    EigRange(..),

    -- * Re-export of Complex from Data.Complex
    module Data.Complex,

    -- * Re-export of Storable from Foreign.Storable
    module Foreign.Storable,
        
    ) where

import Foreign.VMath( VNum, VFractional, VFloating )
import Foreign.BLAS( Trans(..), Uplo(..), Side(..), Diag(..), BLAS1, BLAS2, BLAS3 )
import Foreign.LAPACK( EigJob(..), EigRange(..), LAPACK )
import Data.Complex( Complex(..) )
import Foreign.Storable( Storable() )
