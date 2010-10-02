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
    module Numeric.LinearAlgebra.Types,
    module Numeric.LinearAlgebra.Vector,
    module Numeric.LinearAlgebra.Matrix,
    module Numeric.LinearAlgebra.Packed,    
    module Numeric.LinearAlgebra.Factor.Cholesky,
    module Numeric.LinearAlgebra.Statistics,
    ) where

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Vector( Vector, RVector, STVector )
import Numeric.LinearAlgebra.Matrix( Matrix, RMatrix, STMatrix )
import Numeric.LinearAlgebra.Packed( Packed, RPacked, STPacked )
import Numeric.LinearAlgebra.Factor.Cholesky
import Numeric.LinearAlgebra.Statistics
