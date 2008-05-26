{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
module Data.Vector.Dense (
    Vector,
    module BLAS.Vector,
    module BLAS.Tensor.Base,
    module BLAS.Tensor.Scalable,
    module BLAS.Tensor.Immutable,

    -- * Creating vectors
    vector, 
    listVector,

    -- * Special vectors
    basis,

    -- * Augmenting vectors
    subvector,
    subvectorWithStride,

    -- * Norms and dot product
    sumAbs,
    norm2,
    whichMaxAbs,
    (<.>),

    -- * Vector arithmetic
    shift,
    scale,
    invScale,

    -- * Converting to and from lists
    toList,
    fromList,
    
    -- * Casting vectors
    coerceVector,
    
    -- * Unsafe operations
    unsafeVector,
    unsafeSubvector,
    unsafeSubvectorWithStride,
    ) where

import System.IO.Unsafe           

import Data.Vector.Dense.Internal 
import Data.Vector.Dense.Operations

import BLAS.Access
import BLAS.Elem ( BLAS1, BLAS2 )
import BLAS.Vector hiding ( Vector )
import BLAS.Tensor.Base
import BLAS.Tensor.Scalable
import BLAS.Tensor.Immutable


-- | Create a vector with the given dimension and elements.  The elements
-- given in the association list must all have unique indices, otherwise
-- the result is undefined.
vector :: (BLAS1 e) => Int -> [(Int, e)] -> Vector n e
vector n ies = unsafeFreeze $ unsafePerformIO $ newVector n ies
{-# NOINLINE vector #-}

-- | Same as 'vector', but does not range-check the indices.
unsafeVector :: (BLAS1 e) => Int -> [(Int, e)] -> Vector n e
unsafeVector n ies = unsafeFreeze $ unsafePerformIO $ unsafeNewVector n ies
{-# NOINLINE unsafeVector #-}

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.
listVector :: (BLAS1 e) => Int -> [e] -> Vector n e
listVector n es = unsafePerformIO $ newListVector n es
{-# NOINLINE listVector #-}

-- | @basis n i@ creates a vector of dimension @n@ with zeros everywhere but
-- position @i@, where there is a one.
basis :: (BLAS1 e) => Int -> Int -> Vector n e
basis n i = unsafeFreeze $ unsafePerformIO $ newBasis n i
{-# NOINLINE basis #-}

-- | Convert a vector to a list.  Same as @elems@.
toList :: (BLAS1 e) => Vector n e -> [e]
toList = elems

-- | Convert a list to a vector.  @fromList xs = listVector (length xs) xs@.
fromList :: (BLAS1 e) => [e] -> Vector n e
fromList es = listVector (length es) es
{-# INLINE fromList #-}

instance (BLAS1 e) => Scalable (DVector Imm n) e where
    (*>) = scale
        
instance (BLAS2 e) => Num (DVector Imm n e) where
    (+)         = plus
    (-)         = minus
    (*)         = times
    negate      = (*>) (-1)
    abs         = amap abs
    signum      = amap signum
    fromInteger = (constant 1) . fromInteger
    
instance (BLAS2 e) => Fractional (DVector Imm n e) where
    (/)          = divide
    recip        = amap recip
    fromRational = (constant 1) . fromRational 
    
instance (BLAS2 e, Floating e) => Floating (DVector Imm n e) where
    pi    = constant 1 pi
    exp   = amap exp
    sqrt  = amap sqrt
    log   = amap log
    (**)  = azipWith (**)
    sin   = amap sin
    cos   = amap cos
    tan   = amap tan
    asin  = amap asin
    acos  = amap acos
    atan  = amap atan
    sinh  = amap sinh
    cosh  = amap cosh
    tanh  = amap tanh
    asinh = amap asinh
    acosh = amap acosh
    atanh = amap atanh
    