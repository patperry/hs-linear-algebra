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
    module Numeric.LinearAlgebra.Vector.ST,
    module Numeric.LinearAlgebra.Matrix,
    module Numeric.LinearAlgebra.Matrix.ST,
    module Numeric.LinearAlgebra.Matrix.Herm,
    module Numeric.LinearAlgebra.Matrix.Packed,    
    module Numeric.LinearAlgebra.Factor.Cholesky,
    module Numeric.LinearAlgebra.Factor.Eigen,    
    module Numeric.LinearAlgebra.Statistics,
    ) where

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Vector( Vector )
import Numeric.LinearAlgebra.Vector.ST( STVector, RVector )
import Numeric.LinearAlgebra.Matrix( Matrix )
import Numeric.LinearAlgebra.Matrix.ST( STMatrix, RMatrix )
import Numeric.LinearAlgebra.Matrix.Herm
import Numeric.LinearAlgebra.Matrix.Packed( Packed, RPacked, STPacked )
import Numeric.LinearAlgebra.Factor.Cholesky
import Numeric.LinearAlgebra.Factor.Eigen
import Numeric.LinearAlgebra.Statistics
