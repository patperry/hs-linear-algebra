{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Operations
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Operations (
    -- * Copy and swap
    copyVector,
    swapVectors,

    -- * Vector norms and inner products
    -- ** Pure
    sumAbs,
    norm2,
    whichMaxAbs,
    (<.>),
    
    getSumAbs,
    getNorm2,
    getWhichMaxAbs,
    getDot,

    -- * Vector arithmetic
    -- ** Pure
    shift,
    scale,
    invScale,
    add,
    plus,
    minus,
    times,
    divide,
    
    -- ** Impure
    getShifted,
    getScaled,
    getInvScaled,
    getSum,
    getDiff,
    getProduct,
    getRatio,
    
    -- * In-place arithmetic
    doConj,
    scaleBy,
    shiftBy,
    invScaleBy,
    (+=),
    (-=),
    (*=),
    (//=),
    
    -- * Unsafe operations
    unsafeCopyVector,
    unsafeSwapVectors,
    unsafeGetDot,
    
    -- * BLAS calls
    axpy,
    
    ) where
  
import Control.Monad ( forM_ )
import Data.Vector.Dense.Internal
import BLAS.Tensor
import BLAS.Elem.Base ( Elem )
import qualified BLAS.Elem.Base as E

import Foreign ( Ptr )
import System.IO.Unsafe
import Unsafe.Coerce

import BLAS.Internal  ( inlinePerformIO, checkVecVecOp )
import BLAS.C hiding ( copy, swap, iamax, conj, axpy, acxpy )
import qualified BLAS.C as BLAS
import qualified BLAS.C.Types as T

infixl 7 <.>, `times`, `divide`, `scale`, `invScale`
infixl 6 `plus`, `minus`, `shift`
infixl 1 +=, -=, *=, //=


-- | @copyVector dst src@ replaces the elements in @dst@ with the values from @src@.
-- This may result in a loss of precision.
copyVector :: (BLAS1 e) => IOVector n e -> DVector t n e -> IO ()
copyVector x y = checkVecVecOp "copyVector" (dim x) (dim y) >> unsafeCopyVector x y

-- | Same as 'copyVector' but does not check the dimensions of the arguments.
unsafeCopyVector :: (BLAS1 e) => IOVector n e -> DVector t n e -> IO ()
unsafeCopyVector x y 
    | isConj x && isConj y =
        unsafeCopyVector (conj x) (conj y)
    | isConj x || isConj y =
        forM_ [0..(dim x - 1)] $ \i -> do
            unsafeReadElem y i >>= unsafeWriteElem x i
    | otherwise =
        call2 BLAS.copy y x

-- | @swapVectors x y@ replaces the elements in @x@ with the values from @y@, and 
-- replaces the elements in @y@ with the values from @x@.  This may result in 
-- a loss of precision.
swapVectors :: (BLAS1 e) => IOVector n e -> IOVector n e -> IO ()
swapVectors x y = checkVecVecOp "swapVectors" (dim x) (dim y) >> unsafeSwapVectors x y

-- | Same as 'swap' but does not check the dimensions of the arguments.
unsafeSwapVectors :: (BLAS1 e) => IOVector n e -> IOVector n e -> IO ()
unsafeSwapVectors x y
    | isConj x && isConj y =
        unsafeSwapVectors (conj x) (conj y)
    | isConj x || isConj y =
        forM_ [0..(dim x - 1)] $ \i -> do
            tmp <- unsafeReadElem x i
            unsafeReadElem y i >>= unsafeWriteElem x i
            unsafeWriteElem y i tmp
    | otherwise =
        call2 BLAS.swap x y

-- | Gets the sum of the absolute values of the vector entries.
getSumAbs :: (BLAS1 e) => DVector t n e -> IO Double
getSumAbs = call asum
    
-- | Gets the 2-norm of a vector.
getNorm2 :: (BLAS1 e) => DVector t n e -> IO Double
getNorm2 = call nrm2

-- | Gets the index and norm of the element with maximum magnitude.  This is undefined
-- if any of the elements are @NaN@.  It will throw an exception if the 
-- dimension of the vector is 0.
getWhichMaxAbs :: (BLAS1 e) => DVector t n e -> IO (Int, e)
getWhichMaxAbs x =
    case (dim x) of
        0 -> ioError $ userError $ "getWhichMaxAbs of an empty vector"
        _ -> do
            i <- call BLAS.iamax x
            e <- unsafeReadElem x i
            return (i,e)

-- | Computes the dot product of two vectors.
getDot :: (BLAS1 e) => DVector s n e -> DVector t n e -> IO e
getDot x y = checkVecVecOp "dot" (dim x) (dim y) >> unsafeGetDot x y
{-# INLINE getDot #-}

unsafeGetDot :: (BLAS1 e) => DVector s n e -> DVector t n e -> IO e
unsafeGetDot x y =
    case (isConj x, isConj y) of
        (False, False) -> call2 dotc x y
        (True , False) -> call2 dotu x y
        (False, True ) -> call2 dotu x y >>= return . E.conj
        (True , True)  -> call2 dotc x y >>= return . E.conj
{-# INLINE unsafeGetDot #-}

unsafeGetDotDouble :: DVector s n Double -> DVector t n Double -> IO Double
unsafeGetDotDouble x y = call2 dotc x y
{-# INLINE unsafeGetDotDouble #-}

{-# RULES "unsafeGetDot/Double" unsafeGetDot = unsafeGetDotDouble #-}


-- | Create a new vector by adding a value to every element in a vector.
getShifted :: (BLAS1 e) => e -> DVector t n e -> IO (DVector r n e)
getShifted k x = do
    y <- newCopy x
    shiftBy k (unsafeThaw y)
    return (unsafeCoerce y)

-- | Create a new vector by scaling every element by a value.  @scale'k x@
-- is equal to @newCopy' (scale k x)@.
getScaled :: (BLAS1 e) => e -> DVector t n e -> IO (DVector r n e)
getScaled k x = do
    y <- newCopy x
    scaleBy k (unsafeThaw y)
    return (unsafeCoerce y)

-- | Create a new vector by dividing every element by a value.
getInvScaled :: (BLAS1 e) => e -> DVector t n e -> IO (DVector r n e)
getInvScaled k x = do
    y <- newCopy x
    invScaleBy k (unsafeThaw y)
    return (unsafeCoerce y)

-- | Computes the sum of two vectors.
getSum :: (BLAS1 e) => e -> DVector s n e -> e -> DVector t n e -> IO (DVector r n e)
getSum alpha x beta y = checkVecVecOp "getSum" (dim x) (dim y) >> unsafeGetSum alpha x beta y

unsafeGetSum :: (BLAS1 e) => e -> DVector s n e -> e -> DVector t n e -> IO (DVector r n e)
unsafeGetSum 1 x beta y
    | beta /= 1 = unsafeGetSum beta y 1 x
unsafeGetSum alpha x beta y
    | isConj x = do
        s <- unsafeGetSum (E.conj alpha) (conj x) (E.conj beta) (conj y)
        return (conj s)
    | otherwise = do
        s <- newCopy y
        scaleBy beta (unsafeThaw s)
        axpy alpha x (unsafeThaw s)
        return (unsafeCoerce s)
            
-- | Computes the difference of two vectors.
getDiff :: (BLAS1 e) => DVector s n e -> DVector t n e -> IO (DVector r n e)
getDiff x y = checkVecVecOp "getDiff" (dim x) (dim y) >> unsafeGetSum 1 x (-1) y

-- | Computes the elementwise product of two vectors.
getProduct :: (BLAS2 e) => DVector s n e -> DVector t n e -> IO (DVector r n e)
getProduct = binaryOp "getProduct" (*=)

-- | Computes the elementwise ratio of two vectors.
getRatio :: (BLAS2 e) => DVector s n e -> DVector t n e -> IO (DVector r n e)
getRatio = binaryOp "getRatio" (//=)

-- | Replace every element with its complex conjugate.  See also 'conj'.
doConj  :: (BLAS1 e) => IOVector n e -> IO ()
doConj = call BLAS.conj
            
-- | Add a value to every element in a vector.
shiftBy :: (BLAS1 e) => e -> IOVector n e -> IO ()
shiftBy alpha x | isConj x  = shiftBy (E.conj alpha) (conj x)
                | otherwise = modifyWith (alpha+) x

-- | Scale every element in the vector by the given value, and return a view
-- to the scaled vector.  See also 'scale'.
scaleBy :: (BLAS1 e) => e -> IOVector n e -> IO ()
scaleBy 1 _ = return ()
scaleBy k x | isConj x  = scaleBy (E.conj k) (conj x)
            | otherwise = call (flip scal k) x

-- | Divide every element by a value.
invScaleBy :: (BLAS1 e) => e -> IOVector n e -> IO ()
invScaleBy k x | isConj x  = invScaleBy (E.conj k) (conj x)
               | otherwise = modifyWith (/k) x

-- | @y += x@ replaces @y@ by @y + x@.
(+=) :: (BLAS1 e) => IOVector n e -> DVector t n e -> IO ()
(+=) y x = checkVecVecOp "(+=)" (dim y) (dim x) >> axpy 1 x y

axpy :: (BLAS1 e) => e -> DVector t n e -> IOVector n e -> IO ()
axpy alpha x y
    | isConj y =
        axpy (E.conj alpha) (conj x) (conj y)
    | isConj x =
        call2 (flip BLAS.acxpy alpha) x y
    | otherwise =
        call2 (flip BLAS.axpy alpha) x y

-- | @y -= x@ replaces @y@ by @y - x@.
(-=) :: (BLAS1 e) => IOVector n e -> DVector t n e -> IO ()
(-=) y x = checkVecVecOp "(-=)" (dim y) (dim x) >> axpy (-1) x y

-- | @y *= x@ replaces @y@ by @x * y@, the elementwise product.
(*=) :: (BLAS2 e) => IOVector n e -> DVector t n e -> IO ()
(*=) y x = checkVecVecOp "(*=)" (dim y) (dim x) >> timesEquals y x

timesEquals :: (BLAS2 e) => IOVector n e -> DVector t n e -> IO ()
timesEquals y x
    | isConj y =
        timesEquals (conj y) (conj x)
    | isConj x =
        call2 (flip (tbmv T.colMajor T.upper T.conjTrans T.nonUnit) 0) x y    
    | otherwise =
        call2 (flip (tbmv T.colMajor T.upper T.noTrans T.nonUnit) 0) x y

-- | @y //= x@ replaces @y@ by @y / x@, the elementwise ratio.
(//=) :: (BLAS2 e) => IOVector n e -> DVector t n e -> IO ()
(//=) y x = checkVecVecOp "(//=)" (dim y) (dim x) >> divideEquals y x

divideEquals :: (BLAS2 e) => IOVector n e -> DVector t n e -> IO ()
divideEquals y x
    | isConj y =
        divideEquals (conj y) (conj x)
    | isConj x =
        call2 (flip (tbsv T.colMajor T.upper T.conjTrans T.nonUnit) 0) x y
    | otherwise =
        call2 (flip (tbsv T.colMajor T.upper T.noTrans T.nonUnit) 0) x y
               
call :: (Elem e) => (Int -> Ptr e -> Int -> IO a) -> DVector t n e -> IO a
call f x =
    let n    = dim x
        incX = strideOf x
    in unsafeWithElemPtr x 0 $ \pX -> f n pX incX
{-# INLINE call #-}

call2 :: (Elem e) => 
       (Int -> Ptr e -> Int -> Ptr e -> Int -> IO a) 
    -> DVector s n e -> DVector t m e -> IO a
call2 f x y =
    let n    = dim x
        incX = strideOf x
        incY = strideOf y
    in unsafeWithElemPtr x 0 $ \pX ->
           unsafeWithElemPtr y 0 $ \pY ->
               f n pX incX pY incY
{-# INLINE call2 #-}

binaryOp :: (BLAS1 e) => String -> (IOVector n e -> DVector t n e -> IO ()) 
    -> DVector s n e -> DVector t n e -> IO (DVector r n e)
binaryOp name f x y =
    checkVecVecOp name (dim x) (dim y) >> do
        x' <- newCopy x >>= return . unsafeThaw
        f x' y
        return $! (unsafeCoerce x')



-- | Compute the sum of absolute values of entries in the vector.
sumAbs :: (BLAS1 e) => Vector n e -> Double
sumAbs x = inlinePerformIO $ getSumAbs x
{-# NOINLINE sumAbs #-}

-- | Compute the 2-norm of a vector.
norm2 :: (BLAS1 e) => Vector n e -> Double
norm2 x = inlinePerformIO $ getNorm2 x
{-# NOINLINE norm2 #-}

-- | Get the index and norm of the element with absulte value.  Not valid 
-- if any of the vector entries are @NaN@.  Raises an exception if the 
-- vector has length @0@.
whichMaxAbs :: (BLAS1 e) => Vector n e -> (Int, e)
whichMaxAbs x = inlinePerformIO $ getWhichMaxAbs x
{-# NOINLINE whichMaxAbs #-}

-- | Compute the dot product of two vectors.
(<.>) :: (BLAS1 e) => Vector n e -> Vector n e -> e
(<.>) x y = inlinePerformIO $ getDot x y
{-# NOINLINE (<.>) #-}

-- | Add a value to every element in a vector.
shift :: (BLAS1 e) => e -> Vector n e -> Vector n e
shift k x = unsafePerformIO $ getShifted k x
{-# NOINLINE shift #-}

-- | Scale every element by the given value.
scale :: (BLAS1 e) => e -> Vector n e -> Vector n e
scale k x = unsafePerformIO $ getScaled k x
{-# NOINLINE scale #-}

-- | Divide every element by a value.
invScale :: (BLAS1 e) => e -> Vector n e -> Vector n e
invScale k x = unsafePerformIO $ getInvScaled k x
{-# NOINLINE invScale #-}

add :: (BLAS1 e) => e -> Vector n e -> e -> Vector n e -> Vector n e
add alpha x beta y = unsafePerformIO $ getSum alpha x beta y
{-# NOINLINE add #-}

-- | Sum of two vectors.
plus :: (BLAS1 e) => Vector n e -> Vector n e -> Vector n e
plus x y = add 1 x 1 y

-- | Difference of vectors.
minus :: (BLAS1 e) => Vector n e -> Vector n e -> Vector n e
minus x y = unsafePerformIO $ getDiff x y
{-# NOINLINE minus #-}

-- | Elementwise vector product.
times :: (BLAS2 e) => Vector n e -> Vector n e -> Vector n e
times x y = unsafePerformIO $ getProduct x y
{-# NOINLINE times #-}

-- | Elementwise vector ratio.
divide :: (BLAS2 e) => Vector n e -> Vector n e -> Vector n e
divide x y =  unsafePerformIO $ getRatio x y
{-# NOINLINE divide #-}


{-# RULES
"scale/plus"   forall k l x y. plus (scale k x) (scale l y) = add k x l y
"scale1/plus"  forall k x y.   plus (scale k x) y = add k x 1 y
"scale2/plus"  forall k x y.   plus x (scale k y) = add 1 x k y

"scale/minus"  forall k l x y. minus (scale k x) (scale l y) = add k x (-l) y
"scale1/minus" forall k x y.   minus (scale k x) y = add k x (-1) y
"scale2/minus" forall k x y.   minus x (scale k y) = add 1 x (-k) y
  #-}
  