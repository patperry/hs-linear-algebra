{-# LANGUAGE Rank2Types #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Vector
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Immutable dense vectors.

module Numeric.LinearAlgebra.Vector (
    -- * Immutable vectors
    Vector,
    dim,
    
    -- * Vector construction
    fromList,
    zero,
    constant,

    -- * Accessing vectors
    at,
    unsafeAt,
    indices,
    elems,
    assocs,

    -- * Incremental vector updates
    replace,
    unsafeReplace,
    accum,
    unsafeAccum,

    -- * Derived vectors
    map,
    zipWith,
    unsafeZipWith,
    concat,

    -- * Vector views
    slice,
    unsafeSlice,
    splitAt,
    drop,
    take,

    -- * Vector properties
    sumAbs,    
    norm2,
    whichMaxAbs,
    dot,
    unsafeDot,
    kronecker,

    -- * Vector math operations
    -- ** Num
    add,
    addWithScale,
    sub,
    scaleBy,
    mul,
    negate,
    conjugate,
    abs,
    signum,

    -- ** Fractional
    div,
    recip,        

    -- ** Floating
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
    
    -- * Conversions between foreign pointers
    unsafeFromForeignPtr,
    unsafeToForeignPtr,
    
    -- * Mutable interface
    module Numeric.LinearAlgebra.Vector.ST,
    
    -- * Basic multivariate statistics
    module Numeric.LinearAlgebra.Vector.Statistics,

    ) where

import Prelude( Int, Double, ($), (*), return )
import Control.Monad.ST( runST, ST )
import Foreign( Storable )
import Foreign.BLAS( BLAS1, BLAS2 )
import Foreign.VMath( VNum, VFractional, VFloating )

import Numeric.LinearAlgebra.Vector.Base
import Numeric.LinearAlgebra.Vector.ST
import Numeric.LinearAlgebra.Vector.Statistics

infixr 8 `pow`
infixl 7 `div`
infixl 7 `mul`, `scaleBy`, `kronecker`
infixl 6 `add`, `sub`


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

-- | @add x y@ returns @x + y@.
add :: (VNum e) => Vector e -> Vector e -> Vector e
add = result2 addTo

-- | @sub x y@ returns @x - y@.
sub :: (VNum e) => Vector e -> Vector e -> Vector e
sub = result2 subTo

-- | @scaleBy k x@ returns @k * x@.
scaleBy :: (BLAS1 e) => e -> Vector e -> Vector e
scaleBy k x = create $ do
    x' <- newCopy x
    scaleByM_ k x'
    return x'

-- | @addWithScale alpha x y@ return @alpha * x + y@.
addWithScale :: (BLAS1 e) => e -> Vector e -> Vector e -> Vector e
addWithScale alpha x y = create $ do
    y' <- newCopy y
    addWithScaleM_ alpha x y'
    return y'

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

newResult :: (RVector v, Storable e, Storable f)
          => (STVector s f -> v e -> ST s a)
          -> v e
          -> ST s (STVector s f)
newResult f v = do
    n <- getDim v
    z <- new_ n
    _ <- f z v
    return z
{-# INLINE newResult #-}

newResult2 :: (RVector v1, RVector v2, Storable e, Storable f, Storable g)
           => (STVector s g -> v1 e -> v2 f -> ST s a)
           -> v1 e
           -> v2 f
           -> ST s (STVector s g)
newResult2 f v1 v2 = do
    n <- getDim v1
    z <- new_ n
    _ <- f z v1 v2
    return z
{-# INLINE newResult2 #-}
