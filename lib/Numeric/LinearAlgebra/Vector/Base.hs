{-# LANGUAGE Rank2Types #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Vector.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Numeric.LinearAlgebra.Vector.Base (
    Vector,
    dim,
    
    fromList,
    zero,
    constant,
    
    at,
    unsafeAt,
    
    indices,
    elems,
    assocs,
    
    replace,
    unsafeReplace,
    accum,
    unsafeAccum,
    
    map,
    zipWith,
    unsafeZipWith,
    concat,
    
    slice,
    splitAt,
    drop,
    take,
    
    unsafeSlice,
    unsafeFromForeignPtr,
    unsafeToForeignPtr,
    unsafeWith,
    ) where

import Prelude hiding ( concat, drop, map, read, splitAt, take, zipWith, 
    negate, signum, abs, div, recip, sqrt, exp, log, sin, cos, tan, asin,
    acos, atan, sinh, cosh, tanh, asinh, acosh, atanh )
import qualified Prelude as P

import Data.AEq( AEq(..) )
import Data.Vector.Storable( Vector )
import qualified Data.Vector.Storable as Vector
import Foreign( ForeignPtr, Ptr, Storable )
import Text.Printf( printf )


-- | The dimension of a vector.  This is equal to the number of
-- elements in the vector.                          
dim :: (Storable e) => Vector e -> Int
dim = Vector.length
{-# INLINE dim #-}

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.  @fromList n es@ is equivalent to 
-- @fromAssocs n (zip [0..(n-1)] es)@.
fromList :: (Storable e) => Int -> [e] -> Vector e
fromList = Vector.fromListN
{-# INLINE fromList #-}

-- | Create a zero vector of the given dimension with all elements initialized
-- to the given zero.
zero :: (Storable e, Num e) => Int -> Vector e
zero n = constant n 0
{-# NOINLINE zero #-}

-- | Create a vector of the given dimension with all elements initialized
-- to the given value.
constant :: (Storable e) => Int -> e -> Vector e
constant n e = Vector.fromListN n $ repeat e
{-# NOINLINE constant #-}

-- | Returns the element of a vector at the specified index.
at :: (Storable e) => Vector e -> Int -> e
at v i
    | i < 0 || i >= n = error $
        printf "at <vector with dim %d> %d: invalid index" n i
    | otherwise =
        unsafeAt v i
  where
    n = dim v
{-# INLINE at #-}

-- | Same as 'at' but does not range-check argument.
unsafeAt :: (Storable e) => Vector e -> Int -> e
unsafeAt v i = Vector.unsafeIndex v i
{-# INLINE unsafeAt #-}

-- | Get the indices of the elements in the vector, @[ 0..n-1 ]@, where
-- @n@ is the dimension of the vector.
indices :: (Storable e) => Vector e -> [Int]
indices x = [ 0..n-1 ] where n = dim x
{-# INLINE indices #-}

-- | Returns a list of the elements of a vector, in the same order as their
-- indices.
elems :: (Storable e) => Vector e -> [e]
elems = Vector.toList
{-# INLINE elems #-}

-- | Returns the contents of a vector as a list of associations.
assocs :: (Storable e) => Vector e -> [(Int,e)]
assocs x = zip (indices x) (elems x)
{-# INLINE assocs #-}

unsafeReplace :: (Storable e) => Vector e -> [(Int,e)] -> Vector e
unsafeReplace = Vector.unsafeUpd

-- | Create a new vector by replacing the values at the specified indices.
replace :: (Storable e) => Vector e -> [(Int,e)] -> Vector e
replace = (Vector.//)

-- | @accum f@ takes a vector and an association list and accumulates
-- pairs from the list into the vector with the accumulating function @f@.
accum :: (Storable e)
      => (e -> e' -> e) 
            -> Vector e
            -> [(Int, e')]
            -> Vector e
accum = Vector.accum
{-# INLINE accum #-}

-- | Same as 'accum' but does not range-check indices.
unsafeAccum :: (Storable e)
            => (e -> e' -> e)
            -> Vector e
            -> [(Int, e')]
            -> Vector e
unsafeAccum = Vector.unsafeAccum
{-# INLINE unsafeAccum #-}

-- | Construct a new vector by applying a function to every element of
-- a vector.
map :: (Storable e, Storable e')
    => (e -> e')
    -> Vector e
    -> Vector e'
map = Vector.map
{-# INLINE map #-}

-- | Construct a new vector by applying a function to every pair of elements
-- of two vectors.  The two vectors must have identical dimensions.
zipWith :: (Storable e, Storable e', Storable f)
        => (e -> e' -> f)
        -> Vector e
        -> Vector e'
        -> Vector f
zipWith f v v'
    | n /= n' = error $
        printf ("zipWith <function> <vector with dim %d>"
                ++ " <vector with dim %d>: lengths differ") n n'
    | otherwise =
        unsafeZipWith f v v'
  where
    n  = dim v
    n' = dim v'
{-# INLINE zipWith #-}

-- | Version of 'zipWith' that does not check if the input vectors
-- have the same lengths.
unsafeZipWith :: (Storable e, Storable e', Storable f)
              => (e -> e' -> f)
              -> Vector e
              -> Vector e'
              -> Vector f
unsafeZipWith = Vector.zipWith
{-# INLINE unsafeZipWith #-}

-- | Create a new vector by concatenating a list of vectors.
concat :: (Storable e) => [Vector e] -> Vector e
concat = Vector.concat
{-# INLINE concat #-}

-- | Unsafe version of 'slice' (no range-checking on indices).
unsafeSlice :: (Storable e)
            => Int
            -> Int
            -> Vector e
            -> Vector e            
unsafeSlice = Vector.unsafeSlice
{-# INLINE unsafeSlice #-}

-- | @slice i n v@ creates a subvector view of @v@ starting at
-- index @i@ and having dimension @n@.
slice :: (Storable e)
      => Int
      -> Int
      -> Vector e
      -> Vector e            
slice i n' v
    | i < 0 || n' < 0 || i + n' > n = error $
        printf "slice %d %d <vector with dim %d>: index out of range"
               i n' n
    | otherwise =
        unsafeSlice i n' v
  where
    n = dim v
{-# INLINE slice #-}

-- | @drop i v@ is equal to @slice i (n-i) v@, where @n@ is
-- the dimension of the vector.
drop :: (Storable e) => Int -> Vector e -> Vector e
drop i v = slice i (dim v - i) v
{-# INLINE drop #-}

-- | @take n v@ is equal to @slice 0 n v@.
take :: (Storable e) => Int -> Vector e -> Vector e
take n v = slice 0 n v
{-# INLINE take #-}

-- | Split a vector into two blocks and returns views into the blocks.  If
-- @(v1, v2) = splitAt i v@, then
-- @v1 = slice 0 i v@ and
-- @v2 = slice i (dim v - i) v@.
splitAt :: (Storable e) => Int -> Vector e -> (Vector e, Vector e)
splitAt i v 
    | i < 0 || i > n = error $
        printf "splitAt %d <vector with dim %d>: invalid index" i n
    | otherwise = let
        v1 = unsafeSlice 0 i     v
        v2 = unsafeSlice i (n-i) v
        in (v1, v2)
  where
    n = dim v

-- | Create a vector from a 'ForeignPtr' with an offset and a dimension. The
-- data may not be modified through the ForeignPtr afterwards.
unsafeFromForeignPtr :: (Storable e)
                     => ForeignPtr e -- ^ pointer
                     -> Int	         -- ^ offset
                     -> Int	         -- ^ dimension
                     -> Vector e
unsafeFromForeignPtr = Vector.unsafeFromForeignPtr
{-# INLINE unsafeFromForeignPtr #-}

-- | Yield the underlying 'ForeignPtr' together with the offset to the data
-- and its length. The data may not be modified through the 'ForeignPtr'.
unsafeToForeignPtr :: (Storable e) => Vector e -> (ForeignPtr e, Int, Int)
unsafeToForeignPtr = Vector.unsafeToForeignPtr
{-# INLINE unsafeToForeignPtr #-}

-- | Pass a pointer to the vector's data to the 'IO' action. The data may
-- not be modified through the 'Ptr'.
unsafeWith :: (Storable e) => Vector e -> (Ptr e -> IO a) -> IO a
unsafeWith = Vector.unsafeWith


instance (Storable e, AEq e) => AEq (Vector e) where
    (===) = compareWith (===)
    {-# INLINE (===) #-}
    (~==) = compareWith (~==)
    {-# INLINE (~==) #-}

compareWith :: (Storable e, Storable e')
            => (e -> e' -> Bool)
            -> Vector e
            -> Vector e'
            -> Bool
compareWith cmp v v' =
    dim v == dim v'
    && and (P.zipWith cmp (elems v) (elems v'))
{-# INLINE compareWith #-}

