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
    (*>),
    invScale,
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

infixl 7 <.>, `times`, `divide`, *>, `invScale`
infixl 6 `plus`, `minus`, `shift`
infixl 1 +=, -=, *=, //=


-- | @copyVector dst src@ replaces the elements in @dst@ with the values from @src@.
-- This may result in a loss of precision.
copyVector :: (BLAS1 e) => IOVector n e -> DVector t n e -> IO ()
copyVector x y = checkVecVecOp "copyVector" (dim x) (dim y) >> unsafeCopyVector x y

-- | Same as 'copyVector' but does not check the dimensions of the arguments.
unsafeCopyVector :: (BLAS1 e) => IOVector n e -> DVector t n e -> IO ()
unsafeCopyVector (C x) (C y) =
    unsafeCopyVector x y
unsafeCopyVector x@(DV _ _ _ _) y@(DV _ _ _ _) =
    call2 BLAS.copy y x
unsafeCopyVector x y = do
    forM_ [0..(dim x - 1)] $ \i -> do
        unsafeReadElem y i >>= unsafeWriteElem x i


-- | @swapVectors x y@ replaces the elements in @x@ with the values from @y@, and 
-- replaces the elements in @y@ with the values from @x@.  This may result in 
-- a loss of precision.
swapVectors :: (BLAS1 e) => IOVector n e -> IOVector n e -> IO ()
swapVectors x y = checkVecVecOp "swapVectors" (dim x) (dim y) >> unsafeSwapVectors x y

-- | Same as 'swap' but does not check the dimensions of the arguments.
unsafeSwapVectors :: (BLAS1 e) => IOVector n e -> IOVector n e -> IO ()
unsafeSwapVectors (C x) (C y) =
    unsafeSwapVectors x y
unsafeSwapVectors x@(DV _ _ _ _) y@(DV _ _ _ _) =
    call2 BLAS.swap x y
unsafeSwapVectors x y = do
    forM_ [0..(dim x - 1)] $ \i -> do
        tmp <- unsafeReadElem x i
        unsafeReadElem y i >>= unsafeWriteElem x i
        unsafeWriteElem y i tmp

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

unsafeGetDot :: (BLAS1 e) => DVector s n e -> DVector t n e -> IO e
unsafeGetDot x@(DV _ _ _ _) y@(DV _ _ _ _) =
    call2 dotc x y
unsafeGetDot (C x@(DV _ _ _ _)) (y@(DV _ _ _ _)) =
    call2 dotu x y
unsafeGetDot (x@(DV _ _ _ _)) (C y@(DV _ _ _ _)) =
    call2 dotu x y >>= return . E.conj
unsafeGetDot x@(DV _ _ _ _) (C (C y)) = 
    unsafeGetDot x y
unsafeGetDot (C x) y = 
    unsafeGetDot x (conj y) >>= return . E.conj

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
unsafeGetSum alpha (C x) beta y = do
    s <- unsafeGetSum (E.conj alpha) x (E.conj beta) (conj y)
    return (conj s)
unsafeGetSum alpha x@(DV _ _ _ _) beta y = do
    s <- newCopy x
    scaleBy alpha (unsafeThaw s)
    axpy beta y (unsafeThaw s)
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
shiftBy alpha (C x) = shiftBy (E.conj alpha) x
shiftBy alpha x     = modifyWith (alpha+) x

-- | Scale every element in the vector by the given value, and return a view
-- to the scaled vector.  See also 'scale'.
scaleBy :: (BLAS1 e) => e -> IOVector n e -> IO ()
scaleBy 1 _ = return ()
scaleBy k (C x) = do
    scaleBy (E.conj k) x
scaleBy k x@(DV _ _ _ _) =
    call (flip scal k) x

-- | Divide every element by a value.
invScaleBy :: (BLAS1 e) => e -> IOVector n e -> IO ()
invScaleBy k (C x) = invScaleBy (E.conj k) x
invScaleBy k x     = modifyWith (/k) x

-- | @y += x@ replaces @y@ by @y + x@.
(+=) :: (BLAS1 e) => IOVector n e -> DVector t n e -> IO ()
(+=) y x = checkVecVecOp "(+=)" (dim y) (dim x) >> axpy 1 x y

axpy :: (BLAS1 e) => e -> DVector t n e -> IOVector n e -> IO ()
axpy alpha x@(DV _ _ _ _) y@(DV _ _ _ _) =
    call2 (flip BLAS.axpy alpha) x y
axpy alpha (C x@(DV _ _ _ _))  y@(DV _ _ _ _) =
    call2 (flip BLAS.acxpy alpha) x y
axpy alpha (C (C x)) y = 
    axpy alpha x y
axpy alpha x (C y) =
    axpy (E.conj alpha) (conj x) y

-- | @y -= x@ replaces @y@ by @y - x@.
(-=) :: (BLAS1 e) => IOVector n e -> DVector t n e -> IO ()
(-=) y x = checkVecVecOp "(-=)" (dim y) (dim x) >> axpy (-1) x y

-- | @y *= x@ replaces @y@ by @x * y@, the elementwise product.
(*=) :: (BLAS2 e) => IOVector n e -> DVector t n e -> IO ()
(*=) y x = checkVecVecOp "(*=)" (dim y) (dim x) >> timesEquals y x

timesEquals :: (BLAS2 e) => IOVector n e -> DVector t n e -> IO ()
timesEquals y@(DV _ _ _ _) x@(DV _ _ _ _) =
    call2 (flip (tbmv T.colMajor T.upper T.noTrans T.nonUnit) 0) x y
timesEquals y@(DV _ _ _ _) (C x@(DV _ _ _ _)) =
    call2 (flip (tbmv T.colMajor T.upper T.conjTrans T.nonUnit) 0) x y    
timesEquals y@(DV _ _ _ _) (C (C x)) = 
    timesEquals y x
timesEquals (C y) x =
    timesEquals y (conj x)

-- | @y //= x@ replaces @y@ by @y / x@, the elementwise ratio.
(//=) :: (BLAS2 e) => IOVector n e -> DVector t n e -> IO ()
(//=) y x = checkVecVecOp "(//=)" (dim y) (dim x) >> divideEquals y x

divideEquals :: (BLAS2 e) => IOVector n e -> DVector t n e -> IO ()
divideEquals y@(DV _ _ _ _) x@(DV _ _ _ _) =
    call2 (flip (tbsv T.colMajor T.upper T.noTrans T.nonUnit) 0) x y
divideEquals y@(DV _ _ _ _) (C x@(DV _ _ _ _)) =
    call2 (flip (tbsv T.colMajor T.upper T.conjTrans T.nonUnit) 0) x y
divideEquals y@(DV _ _ _ _) (C (C x)) = 
    divideEquals y x
divideEquals (C y) x =
    divideEquals y (conj x)
               
call :: (Elem e) => (Int -> Ptr e -> Int -> IO a) -> DVector t n e -> IO a
call f x =
    let n    = dim x
        incX = strideOf x
    in unsafeWithElemPtr x 0 $ \pX -> f n pX incX

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
(*>) :: (BLAS1 e) => e -> Vector n e -> Vector n e
(*>) k x = unsafePerformIO $ getScaled k x
{-# NOINLINE (*>) #-}

-- | Divide every element by a value.
invScale :: (BLAS1 e) => e -> Vector n e -> Vector n e
invScale k x = unsafePerformIO $ getInvScaled k x
{-# NOINLINE invScale #-}

-- | Sum of two vectors.
plus :: (BLAS1 e) => Vector n e -> Vector n e -> Vector n e
plus x y = unsafePerformIO $ getSum 1 x 1 y
{-# NOINLINE plus #-}

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
