-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Linear algebra types and operations
--

module Numeric.LinearAlgebra (
    -- * Vector types
    Vector,
    RVector,
    STVector,
    IOVector,
    
    -- * Matrix types
    Matrix,
    RMatrix,
    STMatrix,
    IOMatrix,
    
    -- * Packed matrix types
    Packed,
    RPacked,
    STPacked,
    IOPacked,
    
    module Numeric.LinearAlgebra.Types,
    ) where

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Vector( Vector, RVector, STVector, IOVector )
import Numeric.LinearAlgebra.Matrix( Matrix, RMatrix, STMatrix, IOMatrix )
import Numeric.LinearAlgebra.Packed( Packed, RPacked, STPacked, IOPacked )
