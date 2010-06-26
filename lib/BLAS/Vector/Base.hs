{-# LANGUAGE DeriveDataTypeable, GeneralizedNewtypeDeriving, Rank2Types #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Vector.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module BLAS.Vector.Base
    where

import Control.Monad
import Control.Monad.ST
import Data.AEq( AEq(..) )
import Data.Typeable
import Foreign( advancePtr, peek, peekElemOff, touchForeignPtr )
import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import BLAS.Internal( inlinePerformIO )



import BLAS.Elem
-- import qualified Data.Elem.BLAS.Level1 as BLAS

import BLAS.Vector.STBase


-- | Immutable dense vectors. The type arguments are as follows:
--
--     * @e@: the element type of the vector.
--
newtype Vector e = Vector { unVector :: STVector RealWorld e }
    deriving (RVector, Typeable)

-- | A safe way to create and work with a mutable vector before returning 
-- an immutable vector for later perusal. This function avoids copying
-- the vector before returning it - it uses 'unsafeFreezeVector' internally,
-- but this wrapper is a safe interface to that function. 
runVector :: (forall s . ST s (STVector s e)) -> Vector e
runVector mx = runST $ mx >>= unsafeFreezeVector
{-# INLINE runVector #-}


-- | Converts a mutable vector to an immutable one by taking a complete
-- copy of it.
freezeVector :: (Storable e) => STVector s e -> ST s (Vector e)
freezeVector = fmap Vector . unsafeCoerce . newCopyVector
{-# INLINE freezeVector #-}

-- | Converts a mutable vector into an immutable vector. This simply casts
-- the vector from one type to the other without copying the vector.
--
-- Note that because the vector is possibly not copied, any subsequent
-- modifications made to the mutable version of the vector may be shared with
-- the immutable version. It is safe to use, therefore, if the mutable
-- version is never modified after the freeze operation.
unsafeFreezeVector :: STVector s e -> ST s (Vector e)
unsafeFreezeVector = return . Vector . unsafeCoerce
{-# INLINE unsafeFreezeVector #-}

-- | Converts an immutable vector to a mutable one by taking a complete
-- copy of it.
thawVector :: (Storable e) => Vector e -> ST s (STVector s e)
thawVector = newCopyVector
{-# INLINE thawVector #-}

-- | Converts an immutable vector into a mutable vector. This simply casts
-- the vector from one type to the other without copying the vector.
--
-- Note that because the vector is possibly not copied, any subsequent
-- modifications made to the mutable version of the vector may be shared with
-- the immutable version. It is only safe to use, therefore, if the immutable
-- vector is never referenced again in this thread, and there is no
-- possibility that it can be also referenced in another thread.
unsafeThawVector :: Vector e -> ST s (STVector s e)
unsafeThawVector = return . unsafeCoerce . unVector
{-# INLINE unsafeThawVector #-}


-- | Create a vector with the given dimension and elements.  The elements
-- given in the association list must all have unique indices, otherwise
-- the result is undefined.
--
-- Not every index within the bounds of the array need appear in the
-- association list, but the values associated with indices that do not
-- appear will be undefined.
vector :: (Storable e) => Int -> [(Int, e)] -> Vector e
vector n ies = runVector $ do
    v <- newVector_ n
    setAssocsVector v ies
    return v
{-# INLINE vector #-}

-- | Same as 'vector', but does not range-check the indices.
unsafeVector :: (Storable e) => Int -> [(Int, e)] -> Vector e
unsafeVector n ies = runVector $ do
    v <- newVector_ n
    unsafeSetAssocsVector v ies
    return v
{-# INLINE unsafeVector #-}

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.  @listVector n es@ is equivalent to 
-- @vector n (zip [0..(n-1)] es)@.
listVector :: (Storable e) => Int -> [e] -> Vector e
listVector n es = runVector $ do
    v <- newVector_ n
    setElemsVector v es
    return v
{-# INLINE listVector #-}

-- | Create a vector of the given dimension with all elements initialized
-- to the given value
constantVector :: (Storable e) => Int -> e -> Vector e
constantVector n e = runVector $ newVector n e
{-# INLINE constantVector #-}

-- | Returns the element of a vector at the specified index.
atVector :: (Storable e) => Vector e -> Int -> e
atVector v i
    | i < 0 || i >= n = error $
        printf "indexVector <vector with dimension %d> %d: invalid index" n i
    | otherwise =
        unsafeAtVector v i
  where
    n = dimVector v
{-# INLINE atVector #-}

unsafeAtVector :: (Storable e) => Vector e -> Int -> e
unsafeAtVector (Vector (STVector p _ f)) i = inlinePerformIO $ do
    e  <- peekElemOff p i
    touchForeignPtr f
    return $! e
{-# INLINE unsafeAtVector #-}

-- | Returns a list of the elements of a vector, in the same order as their
-- indices.
elemsVector :: (Storable e) => Vector e -> [e]
elemsVector (Vector (STVector p n f)) =
    let end = p `advancePtr` n
        go p' | p' == end = inlinePerformIO $ do
                                touchForeignPtr f
                                return []
              | otherwise = let e  = inlinePerformIO (peek p')
                                es = go (p' `advancePtr` 1)
                            in e `seq` (e:es)
    in go p
{-# SPECIALIZE INLINE elemsVector :: Vector Double -> [Double] #-}
{-# SPECIALIZE INLINE elemsVector :: Vector (Complex Double) -> [Complex Double] #-}

-- | Returns the contents of a vector as a list of associations.
assocsVector :: (Storable e) => Vector e -> [(Int,e)]
assocsVector x = zip (indicesVector x) (elemsVector x)
{-# INLINE assocsVector #-}

unsafeReplaceVector :: (Storable e) => Vector e -> [(Int,e)] -> Vector e
unsafeReplaceVector v ies = runVector $ do
    v' <- newCopyVector v
    unsafeSetAssocsVector v' ies
    return v'

replaceVector :: (Storable e) => Vector e -> [(Int,e)] -> Vector e
replaceVector v ies = runVector $ do
    v' <- newCopyVector v
    setAssocsVector v' ies
    return v'

-- | @accumVector f@ takes a vector and an association list and accumulates
-- pairs from the list into the vector with the accumulating function @f@.
accumVector :: (Storable e)
            => (e -> e' -> e) 
            -> Vector e
            -> [(Int, e')]
            -> Vector e
accumVector f v ies = runVector $ do
    v' <- newCopyVector v
    forM_ ies $ \(i,new) -> do
        old <- readVector v' i
        unsafeWriteVector v' i (f old new) -- index checked on prev. line
    return v'

-- | Same as 'accumVector' but does not range-check indices.
unsafeAccumVector :: (Storable e)
                  => (e -> e' -> e)
                  -> Vector e
                  -> [(Int, e')]
                  -> Vector e
unsafeAccumVector f v ies = runVector $ do
    v' <- newCopyVector v
    forM_ ies $ \(i,new) -> do
        old <- unsafeReadVector v' i
        unsafeWriteVector v' i (f old new)
    return v'   

-- | Construct a new vector by applying a function to every element of
-- a vector.
mapVector :: (Storable e, Storable e')
          => (e -> e')
          -> Vector e
          -> Vector e'
mapVector f v = runVector $ do
    v' <- newVector_ (dimVector v)
    unsafeMapToVector f v v'
    return v'

-- | Construct a new vector by applying a function to every pair of elements
-- of two vectors.  The two vectors must have identical dimensions.
zipWithVector :: (Storable e, Storable e', Storable f)
              => (e -> e' -> f)
              -> Vector e
              -> Vector e'
              -> Vector f
zipWithVector f v v'
    | n /= n' = error $
        printf ("zipWithVector <function> <vector with dimension %d>"
                ++ " <vector with dimension %d>: lengths differ") n n'
    | otherwise =
        unsafeZipWithVector f v v'
  where
    n  = dimVector v
    n' = dimVector v'
{-# INLINE zipWithVector #-}

-- | Version of 'zipWithVector' that does not check if the input vectors
-- have the same lengths.
unsafeZipWithVector :: (Storable e, Storable e', Storable f)
                    => (e -> e' -> f)
                    -> Vector e
                    -> Vector e'
                    -> Vector f
unsafeZipWithVector f v v' =
    listVector (dimVector v) $ zipWith f (elemsVector v) (elemsVector v')

-- | Compute the sum of absolute values of entries in the vector.
sumAbsVector :: (BLAS1 e) => Vector e -> Double
sumAbsVector v = runST $ getSumAbsVector v
{-# INLINE sumAbsVector #-}

-- | Compute the 2-norm (Euclidean norm) of a vector.
norm2Vector :: (BLAS1 e) => Vector e -> Double
norm2Vector v = runST $ getNorm2Vector v
{-# INLINE norm2Vector #-}

-- | Get the index and norm of the element with absulte value.  Not valid 
-- if any of the vector entries are @NaN@.  Raises an exception if the 
-- vector has length @0@.
whichMaxAbsVector :: (BLAS1 e) => Vector e -> (Int, e)
whichMaxAbsVector v = runST $ getWhichMaxAbsVector v
{-# INLINE whichMaxAbsVector #-}

-- | Compute the dot product of two vectors.
dotVector :: (BLAS1 e) => Vector e -> Vector e -> e
dotVector v v' = runST $ getDotVector v v'
{-# INLINE dotVector #-}

unsafeDotVector :: (BLAS1 e) => Vector e -> Vector e -> e
unsafeDotVector v v' = runST $ unsafeGetDotVector v v'
{-# INLINE unsafeDotVector #-}


instance (Storable e, Show e) => Show (Vector e) where
    show x = "listVector " ++ show (dimVector x) ++ " " ++ show (elemsVector x)
    {-# INLINE show #-}

instance (Storable e, Eq e) => Eq (Vector e) where
    (==) = compareVectorWith (==)
    {-# INLINE (==) #-}

instance (Storable e, AEq e) => AEq (Vector e) where
    (===) = compareVectorWith (===)
    {-# INLINE (===) #-}
    (~==) = compareVectorWith (~==)
    {-# INLINE (~==) #-}

compareVectorWith :: (Storable e, Storable e')
                  => (e -> e' -> Bool)
                  -> Vector e
                  -> Vector e'
                  -> Bool
compareVectorWith cmp v v' =
    dimVector v == dimVector v'
    && and (zipWith cmp (elemsVector v) (elemsVector v'))
{-# INLINE compareVectorWith #-}


instance (VNum e) => Num (Vector e) where
    (+) = resultVector2 addToVector
    {-# INLINE (+) #-}
    (-) = resultVector2 subToVector
    {-# INLINE (-) #-}
    (*) = resultVector2 mulToVector
    {-# INLINE (*) #-}
    negate = resultVector negateToVector
    {-# INLINE negate #-}
    abs = resultVector absToVector
    {-# INLINE abs #-}
    signum = resultVector signumToVector
    {-# INLINE signum #-}
    fromInteger = singletonVector . fromInteger
    {-# INLINE fromInteger #-}

instance (VFractional e) => Fractional (Vector e) where
    (/) = resultVector2 divToVector
    {-# INLINE (/) #-}
    recip = resultVector recipToVector
    {-# INLINE recip #-}
    fromRational = singletonVector . fromRational
    {-# INLINE fromRational #-}

instance (VFloating e) => Floating (Vector e) where
    pi = singletonVector pi
    {-# INLINE pi #-}
    exp = resultVector expToVector
    {-# INLINE exp #-}
    sqrt = resultVector sqrtToVector
    {-# INLINE sqrt #-}
    log = resultVector logToVector
    {-# INLINE log #-}
    (**) = resultVector2 powToVector
    {-# INLINE (**) #-}
    sin = resultVector sinToVector
    {-# INLINE sin #-}
    cos = resultVector cosToVector
    {-# INLINE cos #-}
    tan = resultVector tanToVector
    {-# INLINE tan #-}
    asin = resultVector asinToVector
    {-# INLINE asin #-}
    acos = resultVector acosToVector
    {-# INLINE acos #-}
    atan = resultVector atanToVector
    {-# INLINE atan #-}
    sinh = resultVector sinhToVector
    {-# INLINE sinh #-}
    cosh = resultVector coshToVector
    {-# INLINE cosh #-}
    tanh = resultVector tanhToVector
    {-# INLINE tanh #-}
    asinh = resultVector asinhToVector
    {-# INLINE asinh #-}
    acosh = resultVector acoshToVector
    {-# INLINE acosh #-}
    atanh = resultVector atanhToVector
    {-# INLINE atanh #-}

singletonVector :: (Storable e) => e -> Vector e
singletonVector e = listVector 1 [e]
{-# INLINE singletonVector #-}

resultVector :: (Storable f)
             => (forall s . Vector e -> STVector s f -> ST s a)
             -> Vector e
             -> Vector f
resultVector f v = runVector $ newResultVector f v
{-# INLINE resultVector #-}

resultVector2 :: (Storable g)
              => (forall s . Vector e -> Vector f -> STVector s g -> ST s a)
              -> Vector e
              -> Vector f
              -> Vector g
resultVector2 f v1 v2 = runVector $ newResultVector2 f v1 v2
{-# INLINE resultVector2 #-}
