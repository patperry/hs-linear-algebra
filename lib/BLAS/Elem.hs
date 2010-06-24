-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Elem
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Type classes for elements with BLAS support.
--

module BLAS.Elem (
    -- * Vector math types
    VNum,
    VFractional,
    VFloating,

    -- * BLAS element types
    BLAS1,
    BLAS2,
    BLAS3,
    
    -- * Re-export of Complex from Data.Complex
    module Data.Complex,
    -- * Re-export of Storable from Foreign.Storable
    module Foreign.Storable,
    ) where

import BLAS.Elem.VMath( VNum, VFractional, VFloating )
import BLAS.Elem.Level1( BLAS1 )
import BLAS.Elem.Level2( BLAS2 )
import BLAS.Elem.Level3( BLAS3 )
import Data.Complex( Complex(..) )
import Foreign.Storable( Storable() )
