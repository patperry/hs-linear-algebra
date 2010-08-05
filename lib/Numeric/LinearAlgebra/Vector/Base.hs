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
    dimVector,
    
    vector,
    unsafeVector,
    listVector,
    constantVector,
    
    runVector,
    freezeVector,
    thawVector,
    
    atVector,
    unsafeAtVector,
    
    indicesVector,
    elemsVector,
    assocsVector,
    
    replaceVector,
    unsafeReplaceVector,
    accumVector,
    unsafeAccumVector,
    
    mapVector,
    zipWithVector,
    unsafeZipWithVector,
    
    sliceVector,
    splitVectorAt,
    dropVector,
    takeVector,
    
    sumVector,
    sumAbsVector,    
    norm2Vector,
    whichMaxAbsVector,
    dotVector,
    unsafeDotVector,
    kroneckerVector,
    
    shiftVector,
    addVector,
    addVectorWithScale,
    subVector,
    scaleVector,
    mulVector,
    negateVector,
    conjVector,
    absVector,
    signumVector,
    
    divVector,
    recipVector,        

    sqrtVector,
    expVector,
    logVector,
    powVector,
    sinVector,
    cosVector,
    tanVector,
    asinVector,
    acosVector,
    atanVector,
    sinhVector,
    coshVector,
    tanhVector,
    asinhVector,
    acoshVector,
    atanhVector,
    ) where

import Control.Monad
import Control.Monad.ST
import Data.AEq( AEq(..) )

import Data.Vector.Storable( Vector )
import qualified Data.Vector.Storable as Vector
import qualified Data.Vector.Storable.Mutable as STVector

import Text.Printf( printf )

import Numeric.LinearAlgebra.Elem
import Numeric.LinearAlgebra.Vector.STBase

infixr 8 `powVector`
infixl 7 `divVector`
infixl 7 `mulVector`, `scaleVector`, `kroneckerVector`
infixl 6 `addVector`, `shiftVector`, `subVector`

instance RVector Vector where
    dimVector = Vector.length
    {-# INLINE dimVector #-}
    unsafeSliceVector = Vector.unsafeSlice
    {-# INLINE unsafeSliceVector #-}
    unsafeWithVector = Vector.unsafeWith
    {-# INLINE unsafeWithVector #-}


-- | A safe way to create and work with a mutable vector before returning 
-- an immutable vector for later perusal. This function avoids copying
-- the vector before returning it. 
runVector :: (Storable e) => (forall s . ST s (STVector s e)) -> Vector e
runVector = Vector.create
{-# INLINE runVector #-}

-- | Converts a mutable vector to an immutable one by taking a complete
-- copy of it.
freezeVector :: (Storable e) => STVector s e -> ST s (Vector e)
freezeVector mv = do
    mv' <- newVector_ (dimVector mv)
    unsafeCopyToVector mv mv'
    let (f,o,n) = STVector.unsafeToForeignPtr mv'
        v' = Vector.unsafeFromForeignPtr f o n
    return v'
{-# INLINE freezeVector #-}

-- | Converts an immutable vector to a mutable one by taking a complete
-- copy of it.
thawVector :: (Storable e) => Vector e -> ST s (STVector s e)
thawVector = newCopyVector
{-# INLINE thawVector #-}

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
{-# NOINLINE vector #-}

-- | Same as 'vector', but does not range-check the indices.
unsafeVector :: (Storable e) => Int -> [(Int, e)] -> Vector e
unsafeVector n ies = runVector $ do
    v <- newVector_ n
    unsafeSetAssocsVector v ies
    return v
{-# NOINLINE unsafeVector #-}

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.  @listVector n es@ is equivalent to 
-- @vector n (zip [0..(n-1)] es)@.
listVector :: (Storable e) => Int -> [e] -> Vector e
listVector n es = runVector $ do
    v <- newVector_ n
    setElemsVector v es
    return v
{-# NOINLINE listVector #-}

-- | Create a vector of the given dimension with all elements initialized
-- to the given value
constantVector :: (Storable e) => Int -> e -> Vector e
constantVector n e = runVector $ newVector n e
{-# NOINLINE constantVector #-}

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
unsafeAtVector v i = (Vector.unsafeIndex) v i
{-# INLINE unsafeAtVector #-}

-- | Returns a list of the elements of a vector, in the same order as their
-- indices.
elemsVector :: (Storable e) => Vector e -> [e]
elemsVector = Vector.toList
{-# INLINE elemsVector #-}

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
mapVector = Vector.map
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
        Vector.zipWith f v v'
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
unsafeZipWithVector = Vector.zipWith
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

-- | Compute the kronecker product of two vectors.
kroneckerVector :: (VNum e) => Vector e -> Vector e -> Vector e
kroneckerVector x y = runVector $ do
    z <- newVector_ (dimVector x * dimVector y)
    kroneckerToVector x y z
    return z

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
shiftVector k = mapVector (k+)

-- | @addVector x y@ returns @x + y@.
addVector :: (VNum e) => Vector e -> Vector e -> Vector e
addVector = zipWithVector (+)

-- | @addVectorWithScale a x b y@ returns @a*x + b*y@.
addVectorWithScale :: (VNum e) => e -> Vector e -> e -> Vector e -> Vector e
addVectorWithScale a x b y = zipWithVector (\e f -> a * e + b * f) x y

-- | @subVector x y@ returns @x - y@.
subVector :: (VNum e) => Vector e -> Vector e -> Vector e
subVector = zipWithVector (-)

-- | @scaleVector k x@ returns @k * x@.
scaleVector :: (VNum e) => e -> Vector e -> Vector e
scaleVector k = mapVector (k*)

-- | @mulVector x y@ returns @x + y@.
mulVector :: (VNum e) => Vector e -> Vector e -> Vector e
mulVector = zipWithVector (*)

-- | @negateVector x@ returns @-x@.
negateVector :: (VNum e) => Vector e -> Vector e
negateVector = mapVector negate

-- | @conjVector x@ returns @conj(x)@.
conjVector :: (VNum e) => Vector e -> Vector e
conjVector = resultVector conjToVector -- TODO: replace with mapVector

-- | @absVector x@ returns @abs(x)@.
absVector :: (VNum e) => Vector e -> Vector e
absVector = mapVector abs

-- | @signumVector x@ returns @signum(x)@.
signumVector :: (VNum e) => Vector e -> Vector e
signumVector = mapVector signum

-- | @divVector x y@ returns @x / y@.
divVector :: (VFractional e) => Vector e -> Vector e -> Vector e
divVector = zipWithVector (/)

-- | @recipVector x y@ returns @1 / x@.
recipVector :: (VFractional e) => Vector e -> Vector e
recipVector = mapVector recip

-- | @sqrtVector x@ returns @sqrt(x)@.
sqrtVector :: (VFloating e) => Vector e -> Vector e
sqrtVector = mapVector sqrt

-- | @expVector x@ returns @exp(x)@.
expVector :: (VFloating e) => Vector e -> Vector e
expVector = mapVector exp

-- | @logVector x@ returns @log(x)@.
logVector :: (VFloating e) => Vector e -> Vector e
logVector = mapVector log

-- | @powVector x y@ returns @x ** y@.
powVector :: (VFloating e) => Vector e -> Vector e -> Vector e
powVector = zipWithVector (**)

-- | @sinVector x@ returns @sin(x)@.
sinVector :: (VFloating e) => Vector e -> Vector e
sinVector = mapVector sin

-- | @cosVector x@ returns @cos(x)@.
cosVector :: (VFloating e) => Vector e -> Vector e
cosVector = mapVector cos

-- | @tanVector x@ returns @tan(x)@.
tanVector :: (VFloating e) => Vector e -> Vector e
tanVector = mapVector tan

-- | @asinVector x@ returns @asin(x)@.
asinVector :: (VFloating e) => Vector e -> Vector e
asinVector = mapVector asin

-- | @acosVector x@ returns @acos(x)@.
acosVector :: (VFloating e) => Vector e -> Vector e
acosVector = mapVector acos

-- | @atanVector x@ returns @atan(x)@.
atanVector :: (VFloating e) => Vector e -> Vector e
atanVector = mapVector atan

-- | @sinhVector x@ returns @sinh(x)@.
sinhVector :: (VFloating e) => Vector e -> Vector e
sinhVector = mapVector sinh

-- | @coshVector x@ returns @cosh(x)@.
coshVector :: (VFloating e) => Vector e -> Vector e
coshVector = mapVector cosh

-- | @tanhVector x@ returns @tanh(x)@.
tanhVector :: (VFloating e) => Vector e -> Vector e
tanhVector = mapVector tanh

-- | @asinhVector x@ returns @asinh(x)@.
asinhVector :: (VFloating e) => Vector e -> Vector e
asinhVector = mapVector asinh

-- | @acoshVector x@ returns @acosh(x)@.
acoshVector :: (VFloating e) => Vector e -> Vector e
acoshVector = mapVector acosh

-- | @atanhVector x@ returns @atanh(x)@.
atanhVector :: (VFloating e) => Vector e -> Vector e
atanhVector = mapVector atanh


resultVector :: (Storable e, Storable f)
             => (forall s . Vector e -> STVector s f -> ST s a)
             -> Vector e
             -> Vector f
resultVector f v = runVector $ newResultVector f v
{-# INLINE resultVector #-}
