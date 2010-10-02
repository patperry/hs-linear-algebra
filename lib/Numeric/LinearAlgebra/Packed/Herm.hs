-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Packed.Herm
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Hermitian views of packed matrices.
--
module Numeric.LinearAlgebra.Packed.Herm (
    -- * Immutable interface

    -- ** Vector multiplication
    hermMulVector,
    hermMulVectorWithScale,
    hermMulAddVectorWithScales,
    
    -- ** Updates
    hermRank1Update,
    hermRank2Update,

    -- * Mutable interface
    hermCreate,    

    -- ** Vector multiplication
    hermMulToVector,
    hermMulToVectorWithScale,
    hermMulAddToVectorWithScales,
    
    -- ** Updates
    hermRank1UpdateTo,
    hermRank2UpdateTo,

    ) where

import Numeric.LinearAlgebra.Packed.Base
