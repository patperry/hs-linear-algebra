{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Internal
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
module Data.Vector.Dense.Internal (
    -- * The Vector type
    Vector(..),

    -- * Vector shape
    dim,
    coerceVector,
    module Data.Tensor.Class,

    -- * Conjugating vectors
    module Data.Elem.Conj,

    -- * Creating new vectors
    vector, 
    listVector,
    unsafeVector,

    -- * Reading vector elements
    module Data.Tensor.Class.ITensor,

    -- * Special vectors
    zeroVector,
    constantVector,
    basisVector,

    -- * Vector views
    subvector,
    subvectorWithStride,

    -- * Vector properties
    sumAbs,
    norm2,
    whichMaxAbs,
    (<.>),

    -- * Low-level vector properties
    stride,
    isConj,
    ) where

import Data.AEq
import Foreign( Storable )
import System.IO.Unsafe

import Data.Elem.Conj
import Data.Tensor.Class
import Data.Tensor.Class.ITensor
import Data.Tensor.Class.MTensor

import Data.Elem.BLAS ( BLAS1 )
import BLAS.Internal ( inlinePerformIO )
import BLAS.UnsafeIOToM

import Data.Vector.Dense.Class.Internal
import Data.Vector.Dense.Class.Creating
import Data.Vector.Dense.Class.Special
import Data.Vector.Dense.Class.Views
import Data.Vector.Dense.Class.Operations
import Data.Vector.Dense.Class.Properties


infixl 7 <.>



-- | Create a vector with the given dimension and elements.  The elements
-- given in the association list must all have unique indices, otherwise
-- the result is undefined.
vector :: (BLAS1 e) => Int -> [(Int, e)] -> Vector n e
vector n ies = unsafeFreezeIOVector $ unsafePerformIO $ newVector n ies
{-# NOINLINE vector #-}

-- | Same as 'vector', but does not range-check the indices.
unsafeVector :: (BLAS1 e) => Int -> [(Int, e)] -> Vector n e
unsafeVector n ies = unsafeFreezeIOVector $ unsafePerformIO $ unsafeNewVector n ies
{-# NOINLINE unsafeVector #-}


-- | @zeroVector n@ creates a vector of dimension @n@ with all values
-- set to zero.
zeroVector :: (BLAS1 e) => Int -> Vector n e
zeroVector n = unsafeFreezeIOVector $ unsafePerformIO $ newZeroVector n
{-# NOINLINE zeroVector #-}

-- | @constantVector n e@ creates a vector of dimension @n@ with all values
-- set to @e@.
constantVector :: (BLAS1 e) => Int -> e -> Vector n e
constantVector n e = unsafeFreezeIOVector $ unsafePerformIO $ newConstantVector n e
{-# NOINLINE constantVector #-}

-- | @basisVector n i@ creates a vector of dimension @n@ with zeros 
-- everywhere but position @i@, where there is a one.
basisVector :: (BLAS1 e) => Int -> Int -> Vector n e
basisVector n i = unsafeFreezeIOVector $ unsafePerformIO $ newBasisVector n i
{-# NOINLINE basisVector #-}

-- | @subvector x o n@ creates a subvector of @x@ starting at index @o@ 
-- and having length @n@.
subvector :: (BLAS1 e) => Vector n e -> Int -> Int -> Vector n' e
subvector = subvectorView
{-# INLINE subvector #-}

-- | @subvectorWithStride s x o n@ creates a subvector of @x@ starting 
-- at index @o@, having length @n@ and stride @s@.
subvectorWithStride :: (BLAS1 e) =>
    Int -> Vector n e -> Int -> Int -> Vector n' e
subvectorWithStride = subvectorViewWithStride
{-# INLINE subvectorWithStride #-}

-- | Compute the sum of absolute values of entries in the vector.
sumAbs :: (BLAS1 e) => Vector n e -> Double
sumAbs = unsafeLiftVector getSumAbs
{-# NOINLINE sumAbs #-}

-- | Compute the 2-norm of a vector.
norm2 :: (BLAS1 e) => Vector n e -> Double
norm2 = unsafeLiftVector getNorm2
{-# NOINLINE norm2 #-}

-- | Get the index and norm of the element with absulte value.  Not valid 
-- if any of the vector entries are @NaN@.  Raises an exception if the 
-- vector has length @0@.
whichMaxAbs :: (BLAS1 e) => Vector n e -> (Int, e)
whichMaxAbs = unsafeLiftVector getWhichMaxAbs
{-# NOINLINE whichMaxAbs #-}

-- | Compute the dot product of two vectors.
(<.>) :: (BLAS1 e) => Vector n e -> Vector n e -> e
(<.>) = unsafeLiftVector2 getDot
{-# NOINLINE (<.>) #-}
