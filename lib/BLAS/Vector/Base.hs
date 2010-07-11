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


import BLAS.Vector.STBase

infixr 8 `powVector`
infixl 7 `divVector`
infixl 7 `mulVector`, `scaleVector`
infixl 6 `addVector`, `shiftVector`, `subVector`


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
        printf "atVector <vector with dim %d> %d: invalid index" n i
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

-- | Create a new vector by replacing the values at the specified indices.
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
{-# INLINE mapVector #-}

-- | Construct a new vector by applying a function to every pair of elements
-- of two vectors.  The two vectors must have identical dimensions.
zipWithVector :: (Storable e, Storable e', Storable f)
              => (e -> e' -> f)
              -> Vector e
              -> Vector e'
              -> Vector f
zipWithVector f v v'
    | n /= n' = error $
        printf ("zipWithVector <function> <vector with dim %d>"
                ++ " <vector with dim %d>: lengths differ") n n'
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
{-# INLINE unsafeZipWithVector #-}

-- | Compute the sum of the entries in the vector.
sumVector :: (VNum e) => Vector e -> e
sumVector v = runST $ getSumVector v
{-# INLINE sumVector #-}

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

-- | @shiftVector k x@ returns @k + x@.
shiftVector :: (VNum e) => e -> Vector e -> Vector e
shiftVector k = resultVector $ shiftToVector k

-- | @addVector x y@ returns @x + y@.
addVector :: (VNum e) => Vector e -> Vector e -> Vector e
addVector = resultVector2 addToVector

-- | @addVectorWithScale a x b y@ returns @a*x + b*y@.
addVectorWithScale :: (VNum e) => e -> Vector e -> e -> Vector e -> Vector e
addVectorWithScale a x b y =
    (resultVector2 $ \x' y' -> addToVectorWithScale a x' b y') x y

-- | @subVector x y@ returns @x - y@.
subVector :: (VNum e) => Vector e -> Vector e -> Vector e
subVector = resultVector2 subToVector

-- | @scaleVector k x@ returns @k * x@.
scaleVector :: (VNum e) => e -> Vector e -> Vector e
scaleVector k = resultVector $ scaleToVector k

-- | @mulVector x y@ returns @x + y@.
mulVector :: (VNum e) => Vector e -> Vector e -> Vector e
mulVector = resultVector2 mulToVector

-- | @negateVector x@ returns @-x@.
negateVector :: (VNum e) => Vector e -> Vector e
negateVector = resultVector negateToVector

-- | @conjVector x@ returns @conj(x)@.
conjVector :: (VNum e) => Vector e -> Vector e
conjVector = resultVector conjToVector

-- | @absVector x@ returns @abs(x)@.
absVector :: (VNum e) => Vector e -> Vector e
absVector = resultVector absToVector

-- | @signumVector x@ returns @signum(x)@.
signumVector :: (VNum e) => Vector e -> Vector e
signumVector = resultVector signumToVector

-- | @divVector x y@ returns @x / y@.
divVector :: (VFractional e) => Vector e -> Vector e -> Vector e
divVector = resultVector2 divToVector

-- | @recipVector x y@ returns @1 / x@.
recipVector :: (VFractional e) => Vector e -> Vector e
recipVector = resultVector recipToVector

-- | @sqrtVector x@ returns @sqrt(x)@.
sqrtVector :: (VFloating e) => Vector e -> Vector e
sqrtVector = resultVector sqrtToVector

-- | @expVector x@ returns @exp(x)@.
expVector :: (VFloating e) => Vector e -> Vector e
expVector = resultVector expToVector

-- | @logVector x@ returns @log(x)@.
logVector :: (VFloating e) => Vector e -> Vector e
logVector = resultVector logToVector

-- | @powVector x y@ returns @x ** y@.
powVector :: (VFloating e) => Vector e -> Vector e -> Vector e
powVector = resultVector2 powToVector

-- | @sinVector x@ returns @sin(x)@.
sinVector :: (VFloating e) => Vector e -> Vector e
sinVector = resultVector sinToVector

-- | @cosVector x@ returns @cos(x)@.
cosVector :: (VFloating e) => Vector e -> Vector e
cosVector = resultVector cosToVector

-- | @tanVector x@ returns @tan(x)@.
tanVector :: (VFloating e) => Vector e -> Vector e
tanVector = resultVector tanToVector

-- | @asinVector x@ returns @asin(x)@.
asinVector :: (VFloating e) => Vector e -> Vector e
asinVector = resultVector asinToVector

-- | @acosVector x@ returns @acos(x)@.
acosVector :: (VFloating e) => Vector e -> Vector e
acosVector = resultVector acosToVector

-- | @atanVector x@ returns @atan(x)@.
atanVector :: (VFloating e) => Vector e -> Vector e
atanVector = resultVector atanToVector

-- | @sinhVector x@ returns @sinh(x)@.
sinhVector :: (VFloating e) => Vector e -> Vector e
sinhVector = resultVector sinhToVector

-- | @coshVector x@ returns @cosh(x)@.
coshVector :: (VFloating e) => Vector e -> Vector e
coshVector = resultVector coshToVector

-- | @tanhVector x@ returns @tanh(x)@.
tanhVector :: (VFloating e) => Vector e -> Vector e
tanhVector = resultVector tanhToVector

-- | @asinhVector x@ returns @asinh(x)@.
asinhVector :: (VFloating e) => Vector e -> Vector e
asinhVector = resultVector asinhToVector

-- | @acoshVector x@ returns @acosh(x)@.
acoshVector :: (VFloating e) => Vector e -> Vector e
acoshVector = resultVector acoshToVector

-- | @atanhVector x@ returns @atanh(x)@.
atanhVector :: (VFloating e) => Vector e -> Vector e
atanhVector = resultVector atanhToVector


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
