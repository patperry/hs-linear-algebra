-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Mutable
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.IO (
    -- * Mutable dense vector data type
    IOVector,
    DVector(..),
    
    module BLAS.Vector,
    module BLAS.Tensor.Base,
    module BLAS.Tensor.Dense.ReadOnly,
    module BLAS.Tensor.ReadOnly,
    module BLAS.Tensor.Mutable,
    
    -- * Creating vectors
    newVector, 
    newVector_,
    newListVector,
    
    -- * Special vectors
    newBasis,
    setBasis,

    -- * Vector views
    subvector,
    subvectorWithStride,

    module Data.Vector.Dense.Operations,
    
    -- * Conversion to and from @ForeignPtr@s
    fromForeignPtr,
    toForeignPtr,
    isConj,
    strideOf,
    
    -- * Converting between mutable and immutable vectors
    unsafeFreeze,
    unsafeThaw,
    
    -- * Casting vectors
    coerceVector,
    
    -- * Unsafe operations
    unsafeNewVector,
    unsafeWithElemPtr,
    unsafeSubvector,
    unsafeSubvectorWithStride,
    
    ) where

import Data.Vector.Dense.Internal
import Data.Vector.Dense.Operations hiding ( axpy, sumAbs, norm2, whichMaxAbs, 
    (<.>), shift, scale, invScale, plus, minus, times, divide )
import BLAS.Vector hiding ( Vector )
import BLAS.Tensor.Base
import BLAS.Tensor.Dense.ReadOnly
import BLAS.Tensor.ReadOnly
import BLAS.Tensor.Mutable
