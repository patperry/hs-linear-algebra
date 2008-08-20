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

    -- * Creating vectors
    vector, 
    listVector,

    -- * Special vectors
    zeroVector,
    constantVector,
    basisVector,

    -- * Vector Operations
    sumAbs,
    norm2,
    whichMaxAbs,
    (<.>),

    -- * Converting between mutable and immutable vectors
    UnsafeFreezeVector(..),
    UnsafeThawVector(..),
    freezeVector,
    thawVector,
    
    -- * Unsafe operations
    unsafeVector,

    module Data.Vector.Dense.Class.Read,
    module BLAS.Tensor.Immutable,
    module BLAS.Numeric.Immutable,
    
    ) where

import Data.AEq
import System.IO.Unsafe

import BLAS.Elem ( Elem, BLAS1 )
import BLAS.Internal ( UnsafeIOToM(..), inlinePerformIO )
import BLAS.Tensor.Immutable
import BLAS.Numeric.Immutable


import Data.Vector.Dense.IO hiding ( IOVector )
import qualified Data.Vector.Dense.IO as IO
import Data.Vector.Dense.Class.Read hiding ( conjVector )
import Data.Vector.Dense.Class.Base


infixl 7 <.>


newtype Vector n e = V (IO.IOVector n e)


unsafeFreezeIOVector :: IO.IOVector n e -> Vector n e
unsafeFreezeIOVector = V

unsafeThawIOVector :: Vector n e -> IO.IOVector n e
unsafeThawIOVector (V x) = x

class UnsafeFreezeVector x where
    unsafeFreezeVector :: x n e -> Vector n e
    
class UnsafeThawVector x where
    unsafeThawVector :: Vector n e -> x n e
    
instance UnsafeFreezeVector IO.IOVector where
    unsafeFreezeVector = unsafeFreezeIOVector
instance UnsafeThawVector IO.IOVector where
    unsafeThawVector = unsafeThawIOVector
    
freezeVector :: (CopyTensor x y i e m, UnsafeFreezeVector y) =>
    x n e -> m (Vector n e)
freezeVector x = do
    x' <- newCopy x
    return (unsafeFreezeVector x')

thawVector :: (CopyTensor Vector y Int e m) =>
    Vector n e -> m (y n e)
thawVector = newCopy


liftVector :: (IO.IOVector n e -> a) -> Vector n e -> a
liftVector f (V x) = f x
{-# INLINE liftVector #-}

liftVector2 :: 
    (IO.IOVector n e -> IO.IOVector n e -> a) -> 
        Vector n e -> Vector n e -> a
liftVector2 f x = liftVector (liftVector f x)
{-# INLINE liftVector2 #-}

unsafeLiftVector :: (IO.IOVector n e -> IO a) -> Vector n e -> a
unsafeLiftVector f = unsafePerformIO . liftVector f
{-# NOINLINE unsafeLiftVector #-}

unsafeLiftVector2 :: 
    (IO.IOVector n e -> IO.IOVector n e -> IO a) -> 
        Vector n e -> Vector n e -> a
unsafeLiftVector2 f x y = unsafePerformIO $ liftVector2 f x y
{-# NOINLINE unsafeLiftVector2 #-}

inlineLiftVector :: (IO.IOVector n e -> IO a) -> Vector n e -> a
inlineLiftVector f = inlinePerformIO . liftVector f
{-# INLINE inlineLiftVector #-}

-- | Create a vector with the given dimension and elements.  The elements
-- given in the association list must all have unique indices, otherwise
-- the result is undefined.
vector :: (Elem e) => Int -> [(Int, e)] -> Vector n e
vector n ies = unsafeFreezeIOVector $ unsafePerformIO $ newVector n ies
{-# NOINLINE vector #-}

-- | Same as 'vector', but does not range-check the indices.
unsafeVector :: (Elem e) => Int -> [(Int, e)] -> Vector n e
unsafeVector n ies = unsafeFreezeIOVector $ unsafePerformIO $ unsafeNewVector n ies
{-# NOINLINE unsafeVector #-}

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.  @listVector n es@ is equivalent to 
-- @vector n (zip [0..(n-1)] es)@, except that the result is undefined 
-- if @length es@ is less than @n@.
listVector :: (Elem e) => Int -> [e] -> Vector n e
listVector n es = unsafeFreezeIOVector $ unsafePerformIO $ newListVector n es
{-# NOINLINE listVector #-}

-- | @zeroVector n@ creates a vector of dimension @n@ with all values
-- set to zero.
zeroVector :: (Elem e) => Int -> Vector n e
zeroVector n = unsafeFreezeIOVector $ unsafePerformIO $ newZeroVector n
{-# NOINLINE zeroVector #-}

-- | @constantVector n e@ creates a vector of dimension @n@ with all values
-- set to @e@.
constantVector :: (Elem e) => Int -> e -> Vector n e
constantVector n e = unsafeFreezeIOVector $ unsafePerformIO $ newConstantVector n e
{-# NOINLINE constantVector #-}

-- | @basisVector n i@ creates a vector of dimension @n@ with zeros 
-- everywhere but position @i@, where there is a one.
basisVector :: (Elem e) => Int -> Int -> Vector n e
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


instance BaseTensor Vector Int e where
    shape  = liftVector shape
    bounds = liftVector bounds

instance (BLAS1 e) => ITensor Vector Int e where
    constant      = constantVector
    zero          = zeroVector

    (//)          = replaceHelp writeElem
    unsafeReplace = replaceHelp unsafeWriteElem
    
    unsafeAt x i  = inlineLiftVector (flip unsafeReadElem i) x
    {-# INLINE unsafeAt #-}
    
    size          = inlineLiftVector getSize
    elems         = inlineLiftVector getElems
    indices       = inlineLiftVector getIndices
    assocs        = inlineLiftVector getAssocs

    tmap f x      = listVector (dim x) (map f $ elems x)

replaceHelp :: (BLAS1 e) => 
    (IO.IOVector n e -> Int -> e -> IO ()) ->
        Vector n e -> [(Int, e)] -> Vector n e
replaceHelp set x ies =
    unsafePerformIO $ do
        y  <- newCopy (unsafeThawIOVector x)
        mapM_ (uncurry $ set y) ies
        return (unsafeFreezeVector y)
{-# NOINLINE replaceHelp #-}

instance (BLAS1 e, Monad m) => ReadTensor Vector Int e m where
    getSize        = return . size

    getAssocs      = return . assocs
    getIndices     = return . indices
    getElems       = return . elems

    getAssocs'     = getAssocs
    getIndices'    = getIndices
    getElems'      = getElems
    
    unsafeReadElem x i = return (unsafeAt x i)

instance (BLAS1 e) => INumeric Vector Int e where
    (*>) k x = unsafeFreezeIOVector $ unsafeLiftVector (getScaled k) x
    {-# NOINLINE (*>) #-}

    shift k x = unsafeFreezeIOVector $ unsafeLiftVector (getShifted k) x
    {-# NOINLINE shift #-}

instance (BLAS1 e, Monad m) => ReadNumeric Vector Int e m where

instance (Elem e) => BaseVector Vector e where
    stride                    = liftVector stride
    isConj                    = liftVector isConj
    conjVector (V x)          = V $ conjVector x
    unsafeSubvectorWithStride s (V x) o n = V $ unsafeSubvectorWithStride s x o n
    vectorViewArray f o n s c = V $ vectorViewArray f o n s c
    withVectorPtr             = liftVector withVectorPtr

instance (BLAS1 e, UnsafeIOToM m) => ReadVector Vector e m where
    
instance (BLAS1 e) => Num (Vector n e) where
    (+) x y     = unsafeFreezeIOVector $ unsafeLiftVector2 getAdd x y
    (-) x y     = unsafeFreezeIOVector $ unsafeLiftVector2 getSub x y
    (*) x y     = unsafeFreezeIOVector $ unsafeLiftVector2 getMul x y
    negate x    = unsafeFreezeIOVector $ unsafeLiftVector (getScaled (-1)) x
    abs         = tmap abs
    signum      = tmap signum
    fromInteger = (constant 1) . fromInteger
    
instance (BLAS1 e) => Fractional (Vector n e) where
    (/) x y      = unsafeFreezeIOVector $ unsafeLiftVector2 getDiv x y
    recip        = tmap recip
    fromRational = (constant 1) . fromRational 
    
instance (BLAS1 e, Floating e) => Floating (Vector n e) where
    pi    = constant 1 pi
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

            