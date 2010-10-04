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
    unVector,
    unSTVector,
    dim,
    
    fromAssocs,
    unsafeFromAssocs,
    fromList,
    constant,
    
    create,
    freeze,
    unsafeFreeze,
    
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
    concat,
    
    slice,
    splitAt,
    drop,
    take,
    
    sumAbs,    
    norm2,
    whichMaxAbs,
    dot,
    unsafeDot,
    kronecker,

    add,
    sub,
    scale,
    mul,
    negate,
    conjugate,
    abs,
    signum,
    
    div,
    recip,        

    sqrt,
    exp,
    log,
    pow,
    sin,
    cos,
    tan,
    asin,
    acos,
    atan,
    sinh,
    cosh,
    tanh,
    asinh,
    acosh,
    atanh,
    ) where

import Prelude hiding ( concat, drop, map, read, splitAt, take, zipWith, 
    negate, signum, abs, div, recip, sqrt, exp, log, sin, cos, tan, asin,
    acos, atan, sinh, cosh, tanh, asinh, acosh, atanh )
import qualified Prelude as P

import Control.Monad
import Control.Monad.ST
import Data.AEq( AEq(..) )
import Data.List( mapAccumL )

import Data.Vector.Storable( Vector )
import qualified Data.Vector.Storable as Vector
import qualified Data.Vector.Storable.Mutable as STVector

import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Vector.STBase

infixr 8 `pow`
infixl 7 `div`
infixl 7 `mul`, `scale`, `kronecker`
infixl 6 `add`, `sub`

unVector :: Vector e -> STVector RealWorld e
unVector = unsafeCoerce
{-# INLINE unVector #-}

unSTVector :: STVector s e -> Vector e
unSTVector = unsafeCoerce
{-# INLINE unSTVector #-}

instance RVector Vector where
    dim = Vector.length
    {-# INLINE dim #-}
    unsafeSlice = Vector.unsafeSlice
    {-# INLINE unsafeSlice #-}
    unsafeWith = Vector.unsafeWith
    {-# INLINE unsafeWith #-}
    unsafeFromForeignPtr = Vector.unsafeFromForeignPtr
    {-# INLINE unsafeFromForeignPtr #-}
    unsafeToForeignPtr = Vector.unsafeToForeignPtr
    {-# INLINE unsafeToForeignPtr #-}


-- | A safe way to create and work with a mutable vector before returning 
-- an immutable vector for later perusal. This function avoids copying
-- the vector before returning it - it uses 'unsafeFreeze' internally,
-- but this wrapper is a safe interface to that function.
create :: (Storable e) => (forall s . ST s (STVector s e)) -> Vector e
create = Vector.create
{-# INLINE create #-}

-- | Converts a mutable vector to an immutable one by taking a complete
-- copy of it.
freeze :: (Storable e) => STVector s e -> ST s (Vector e)
freeze mv = do
    mv' <- newCopy mv
    unsafeFreeze mv'
{-# INLINE freeze #-}

-- | Converts a mutable vector into an immutable vector. This simply casts
-- the vector from one type to the other without copying the vector.
-- Note that because the vector is possibly not copied, any subsequent
-- modifications made to the mutable version of the vector may be shared with
-- the immutable version. It is safe to use, therefore, if the mutable
-- version is never modified after the freeze operation.
unsafeFreeze :: (Storable e) => STVector s e -> ST s (Vector e)
unsafeFreeze mv =
    let (f,o,n) = STVector.unsafeToForeignPtr mv
        v = Vector.unsafeFromForeignPtr f o n
    in return v
{-# INLINE unsafeFreeze #-}

-- | Create a vector with the given dimension and elements.  The elements
-- given in the association list must all have unique indices, otherwise
-- the result is undefined.
--
-- Not every index within the bounds of the array need appear in the
-- association list, but the values associated with indices that do not
-- appear will be undefined.
fromAssocs :: (Storable e) => Int -> [(Int, e)] -> Vector e
fromAssocs n ies = create $ do
    v <- new_ n
    setAssocs v ies
    return v
{-# NOINLINE fromAssocs #-}

-- | Same as 'fromAssocs', but does not range-check the indices.
unsafeFromAssocs :: (Storable e) => Int -> [(Int, e)] -> Vector e
unsafeFromAssocs n ies = create $ do
    v <- new_ n
    unsafeSetAssocs v ies
    return v
{-# NOINLINE unsafeFromAssocs #-}

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.  @fromList n es@ is equivalent to 
-- @vector n (zip [0..(n-1)] es)@.
fromList :: (Storable e) => Int -> [e] -> Vector e
fromList n es = create $ do
    v <- new_ n
    setElems v es
    return v
{-# NOINLINE fromList #-}

-- | Create a vector of the given dimension with all elements initialized
-- to the given value
constant :: (Storable e) => Int -> e -> Vector e
constant n e = create $ new n e
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
unsafeReplace v ies = create $ do
    v' <- newCopy v
    unsafeSetAssocs v' ies
    return v'

-- | Create a new vector by replacing the values at the specified indices.
replace :: (Storable e) => Vector e -> [(Int,e)] -> Vector e
replace v ies = create $ do
    v' <- newCopy v
    setAssocs v' ies
    return v'

-- | @accum f@ takes a vector and an association list and accumulates
-- pairs from the list into the vector with the accumulating function @f@.
accum :: (Storable e)
            => (e -> e' -> e) 
            -> Vector e
            -> [(Int, e')]
            -> Vector e
accum f v ies = create $ do
    v' <- newCopy v
    forM_ ies $ \(i,e') -> do
        e <- read v' i
        unsafeWrite v' i (f e e') -- index checked on prev. line
    return v'

-- | Same as 'accum' but does not range-check indices.
unsafeAccum :: (Storable e)
                  => (e -> e' -> e)
                  -> Vector e
                  -> [(Int, e')]
                  -> Vector e
unsafeAccum f v ies = create $ do
    v' <- newCopy v
    forM_ ies $ \(i,e') -> do
        e <- unsafeRead v' i
        unsafeWrite v' i (f e e')
    return v'   

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
unsafeZipWith f v v' =
    fromList (dim v) $ P.zipWith f (elems v) (elems v')
{-# INLINE unsafeZipWith #-}

-- | Create a new vector by concatenating a list of vectors.
concat :: (Storable e) => [Vector e] -> Vector e
concat xs = let
    (n_tot,onxs) = mapAccumL (\o x -> let n = dim x 
                                      in o `seq` (o + n, (o,n,x))) 0 xs
    in create $ do
        y <- new_ n_tot
        forM_ onxs $ \(o,n,x) ->
            unsafeCopyTo (unsafeSlice o n y) x
        return y

-- | Compute the sum of absolute values of entries in the vector.
sumAbs :: (BLAS1 e) => Vector e -> Double
sumAbs v = runST $ getSumAbs v
{-# INLINE sumAbs #-}

-- | Compute the 2-norm (Euclidean norm) of a vector.
norm2 :: (BLAS1 e) => Vector e -> Double
norm2 v = runST $ getNorm2 v
{-# INLINE norm2 #-}

-- | Get the index and norm of the element with absulte value.  Not valid 
-- if any of the vector entries are @NaN@.  Raises an exception if the 
-- vector has length @0@.
whichMaxAbs :: (BLAS1 e) => Vector e -> (Int, e)
whichMaxAbs v = runST $ getWhichMaxAbs v
{-# INLINE whichMaxAbs #-}

-- | Compute the dot product of two vectors.
dot :: (BLAS1 e) => Vector e -> Vector e -> e
dot v v' = runST $ getDot v v'
{-# INLINE dot #-}

unsafeDot :: (BLAS1 e) => Vector e -> Vector e -> e
unsafeDot v v' = runST $ unsafeGetDot v v'
{-# INLINE unsafeDot #-}

-- | Compute the kronecker product of two vectors.
kronecker :: (BLAS2 e) => Vector e -> Vector e -> Vector e
kronecker x y = create $ do
    z <- new_ (dim x * dim y)
    kroneckerTo z x y
    return z


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

-- | @add x y@ returns @x + y@.
add :: (VNum e) => Vector e -> Vector e -> Vector e
add = result2 addTo

-- | @sub x y@ returns @x - y@.
sub :: (VNum e) => Vector e -> Vector e -> Vector e
sub = result2 subTo

-- | @scale x k@ returns @k * x@.
scale :: (BLAS1 e) => Vector e -> e -> Vector e
scale x k = create $ do
    x' <- newCopy x
    scaleM x' k
    return x'

-- | @mul x y@ returns @x * y@.
mul :: (VNum e) => Vector e -> Vector e -> Vector e
mul = result2 mulTo

-- | @negate x@ returns @-x@.
negate :: (VNum e) => Vector e -> Vector e
negate = result negateTo

-- | @conjugate x@ returns @conjugate(x)@.
conjugate :: (VNum e) => Vector e -> Vector e
conjugate = result conjugateTo

-- | @abs x@ returns @abs(x)@.
abs :: (VNum e) => Vector e -> Vector e
abs = result absTo

-- | @signum x@ returns @signum(x)@.
signum :: (VNum e) => Vector e -> Vector e
signum = result signumTo

-- | @div x y@ returns @x / y@.
div :: (VFractional e) => Vector e -> Vector e -> Vector e
div = result2 divTo

-- | @recip x y@ returns @1 / x@.
recip :: (VFractional e) => Vector e -> Vector e
recip = result recipTo

-- | @sqrt x@ returns @sqrt(x)@.
sqrt :: (VFloating e) => Vector e -> Vector e
sqrt = result sqrtTo

-- | @exp x@ returns @exp(x)@.
exp :: (VFloating e) => Vector e -> Vector e
exp = result expTo

-- | @log x@ returns @log(x)@.
log :: (VFloating e) => Vector e -> Vector e
log = result logTo

-- | @pow x y@ returns @x ** y@.
pow :: (VFloating e) => Vector e -> Vector e -> Vector e
pow = result2 powTo

-- | @sin x@ returns @sin(x)@.
sin :: (VFloating e) => Vector e -> Vector e
sin = result sinTo

-- | @cos x@ returns @cos(x)@.
cos :: (VFloating e) => Vector e -> Vector e
cos = result cosTo

-- | @tan x@ returns @tan(x)@.
tan :: (VFloating e) => Vector e -> Vector e
tan = result tanTo

-- | @asin x@ returns @asin(x)@.
asin :: (VFloating e) => Vector e -> Vector e
asin = result asinTo

-- | @acos x@ returns @acos(x)@.
acos :: (VFloating e) => Vector e -> Vector e
acos = result acosTo

-- | @atan x@ returns @atan(x)@.
atan :: (VFloating e) => Vector e -> Vector e
atan = result atanTo

-- | @sinh x@ returns @sinh(x)@.
sinh :: (VFloating e) => Vector e -> Vector e
sinh = result sinhTo

-- | @cosh x@ returns @cosh(x)@.
cosh :: (VFloating e) => Vector e -> Vector e
cosh = result coshTo

-- | @tanh x@ returns @tanh(x)@.
tanh :: (VFloating e) => Vector e -> Vector e
tanh = result tanhTo

-- | @asinh x@ returns @asinh(x)@.
asinh :: (VFloating e) => Vector e -> Vector e
asinh = result asinhTo

-- | @acosh x@ returns @acosh(x)@.
acosh :: (VFloating e) => Vector e -> Vector e
acosh = result acoshTo

-- | @atanh x@ returns @atanh(x)@.
atanh :: (VFloating e) => Vector e -> Vector e
atanh = result atanhTo


result :: (Storable e, Storable f)
       => (forall s . STVector s f -> Vector e -> ST s a)
       -> Vector e
       -> Vector f
result f v = create $ newResult f v
{-# INLINE result #-}

result2 :: (Storable e, Storable f, Storable g)
        => (forall s . STVector s g -> Vector e -> Vector f -> ST s a)
        -> Vector e
        -> Vector f
        -> Vector g
result2 f v1 v2 = create $ newResult2 f v1 v2
{-# INLINE result2 #-}
