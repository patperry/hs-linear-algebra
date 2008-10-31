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
    module BLAS.Tensor.Base,

    -- * Conjugating vectors
    module BLAS.Conj,

    -- * Creating new vectors
    vector, 
    listVector,
    unsafeVector,

    -- * Reading vector elements
    module BLAS.Tensor.Immutable,

    -- * Special vectors
    zeroVector,
    constantVector,
    basisVector,

    -- * Vector views
    subvector,
    subvectorWithStride,
    unsafeSubvector,
    unsafeSubvectorWithStride,

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
import System.IO.Unsafe

import BLAS.Conj
import BLAS.Tensor.Base
import BLAS.Tensor.Immutable

import BLAS.Elem ( Elem, BLAS1 )
import BLAS.Internal ( inlinePerformIO )
import BLAS.UnsafeIOToM

import Data.Vector.Dense.IO
import Data.Vector.Dense.Class


infixl 7 <.>


newtype Vector n e = V (IOVector n e)


unsafeFreezeIOVector :: IOVector n e -> Vector n e
unsafeFreezeIOVector = V
unsafeThawIOVector :: Vector n e -> IOVector n e
unsafeThawIOVector (V x) = x


liftVector :: (IOVector n e -> a) -> Vector n e -> a
liftVector f (V x) = f x
{-# INLINE liftVector #-}

liftVector2 :: 
    (IOVector n e -> IOVector n e -> a) -> 
        Vector n e -> Vector n e -> a
liftVector2 f x = liftVector (liftVector f x)
{-# INLINE liftVector2 #-}

unsafeLiftVector :: (IOVector n e -> IO a) -> Vector n e -> a
unsafeLiftVector f = unsafePerformIO . liftVector f
{-# NOINLINE unsafeLiftVector #-}

unsafeLiftVector2 :: 
    (IOVector n e -> IOVector n e -> IO a) -> 
        Vector n e -> Vector n e -> a
unsafeLiftVector2 f x y = unsafePerformIO $ liftVector2 f x y
{-# NOINLINE unsafeLiftVector2 #-}

inlineLiftVector :: (IOVector n e -> IO a) -> Vector n e -> a
inlineLiftVector f = inlinePerformIO . liftVector f
{-# INLINE inlineLiftVector #-}

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

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.  @listVector n es@ is equivalent to 
-- @vector n (zip [0..(n-1)] es)@, except that the result is undefined 
-- if @length es@ is less than @n@.
listVector :: (BLAS1 e) => Int -> [e] -> Vector n e
listVector n es = unsafeFreezeIOVector $ unsafePerformIO $ newListVector n es
{-# NOINLINE listVector #-}

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


instance (Elem e) => BaseTensor Vector Int e where
    shape  = liftVector shape
    bounds = liftVector bounds

instance (BLAS1 e) => ITensor Vector Int e where
    (//)          = replaceHelp writeElem
    unsafeReplace = replaceHelp unsafeWriteElem
    
    unsafeAt x i  = inlineLiftVector (flip unsafeReadElem i) x
    {-# INLINE unsafeAt #-}
    
    size          = inlineLiftVector getSize
    elems         = inlineLiftVector getElems
    indices       = inlineLiftVector getIndices
    assocs        = inlineLiftVector getAssocs

    tmap f x      = listVector (dim x) (map f $ elems x)

    (*>) k x = unsafeFreezeIOVector $ unsafeLiftVector (getScaledVector k) x
    {-# NOINLINE (*>) #-}

    shift k x = unsafeFreezeIOVector $ unsafeLiftVector (getShiftedVector k) x
    {-# NOINLINE shift #-}


replaceHelp :: (BLAS1 e) => 
    (IOVector n e -> Int -> e -> IO ()) ->
        Vector n e -> [(Int, e)] -> Vector n e
replaceHelp set x ies =
    unsafePerformIO $ do
        y  <- newCopyVector (unsafeThawIOVector x)
        mapM_ (uncurry $ set y) ies
        return (unsafeFreezeIOVector y)
{-# NOINLINE replaceHelp #-}

instance (BLAS1 e, Monad m) => ReadTensor Vector Int e m where
    getSize        = return . size
    getAssocs      = return . assocs
    getIndices     = return . indices
    getElems       = return . elems
    getAssocs'     = return . assocs
    getIndices'    = return . indices
    getElems'      = return . elems
    unsafeReadElem x i = return $ unsafeAt x i
    


instance (Elem e) => BaseVector Vector e where
    vectorViewArray f o n s c = V $ vectorViewArray f o n s c
    arrayFromVector           = liftVector arrayFromVector

instance (BLAS1 e, UnsafeIOToM m) => ReadVector Vector e m where
    
instance (BLAS1 e) => Num (Vector n e) where
    (+) x y     = unsafeFreezeIOVector $ unsafeLiftVector2 getAddVector x y
    (-) x y     = unsafeFreezeIOVector $ unsafeLiftVector2 getSubVector x y
    (*) x y     = unsafeFreezeIOVector $ unsafeLiftVector2 getMulVector x y
    negate      = ((-1) *>)
    abs         = tmap abs
    signum      = tmap signum
    fromInteger = (constantVector 1) . fromInteger
    
instance (BLAS1 e) => Fractional (Vector n e) where
    (/) x y      = unsafeFreezeIOVector $ unsafeLiftVector2 getDivVector x y
    recip        = tmap recip
    fromRational = (constantVector 1) . fromRational 
    
instance (BLAS1 e, Floating e) => Floating (Vector n e) where
    pi    = constantVector 1 pi
    exp   = tmap exp
    sqrt  = tmap sqrt
    log   = tmap log
    (**)  = tzipWith (**)
    sin   = tmap sin
    cos   = tmap cos
    tan   = tmap tan
    asin  = tmap asin
    acos  = tmap acos
    atan  = tmap atan
    sinh  = tmap sinh
    cosh  = tmap cosh
    tanh  = tmap tanh
    asinh = tmap asinh
    acosh = tmap acosh
    atanh = tmap atanh

tzipWith :: (BLAS1 e) =>
    (e -> e -> e) -> Vector n e -> Vector n e -> Vector n e
tzipWith f x y
    | dim y /= n =
        error ("tzipWith: vector lengths differ; first has length `" ++
                show n ++ "' and second has length `" ++
                show (dim y) ++ "'")
    | otherwise =
        listVector n (zipWith f (elems x) (elems y))
  where
    n = dim x


instance (BLAS1 e, Show e) => Show (Vector n e) where
    show x
        | isConj x  = "conj (" ++ show (conj x) ++ ")"
        | otherwise = "listVector " ++ show (dim x) ++ " " ++ show (elems x)

instance (BLAS1 e, Eq e) => Eq (Vector n e) where
    (==) = compareHelp (==)

instance (BLAS1 e, AEq e) => AEq (Vector n e) where
    (===) = compareHelp (===)
    (~==) = compareHelp (~==)

compareHelp :: (BLAS1 e) => 
    (e -> e -> Bool) -> 
        Vector n e -> Vector n e -> Bool
compareHelp cmp x y
    | isConj x && isConj y =
        compareHelp cmp (conj x) (conj y)
    | otherwise =
        (dim x == dim y) && (and $ zipWith cmp (elems x) (elems y))

            