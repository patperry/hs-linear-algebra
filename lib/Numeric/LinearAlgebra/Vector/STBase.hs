{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Vector.STBase
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : asinerimental
--

module Numeric.LinearAlgebra.Vector.STBase (
    STVector,
    IOVector,
    RVector(..),
    
    slice,
    drop,
    take,
    splitAt,
    
    new_,
    new,
    newCopy,
    newResult,
    newResult2,
    
    clear,
    copyTo,
    unsafeCopyTo,
    swap,
    unsafeSwap,
    
    indices,
    getElems,
    getElems',
    getAssocs,
    getAssocs',
    setElems,
    setAssocs,
    unsafeSetAssocs,
    
    read,
    unsafeRead,
    write,
    unsafeWrite,
    modify,
    unsafeModify,
    unsafeSwapElems,

    mapTo,
    unsafeMapTo,
    zipWithTo,
    unsafeZipWithTo,
    
    getSumAbs,
    getNorm2,
    getWhichMaxAbs,
    getDot,
    unsafeGetDot,
    kroneckerTo,
    
    shiftTo,
    addTo,
    addToWithScales,
    subTo,
    scaleTo,
    mulTo,
    negateTo,
    conjugateTo,
    absTo,
    signumTo,
    divTo,
    recipTo,        
    sqrtTo,
    expTo,
    logTo,
    powTo,
    sinTo,
    cosTo,
    tanTo,
    asinTo,
    acosTo,
    atanTo,
    sinhTo,
    coshTo,
    tanhTo,
    asinhTo,
    acoshTo,
    atanhTo,
    ) where

import Prelude hiding ( drop, read, splitAt, take )

import Control.Monad
import Control.Monad.ST
import Foreign hiding ( new )
import System.IO.Unsafe( unsafeInterleaveIO )
import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import Data.Vector.Storable.Mutable( STVector, IOVector, MVector )
import qualified Data.Vector.Storable.Mutable as STVector

import Numeric.LinearAlgebra.Internal( clearArray )
import Numeric.LinearAlgebra.Types
import qualified Foreign.BLAS as BLAS
import qualified Foreign.VMath as VMath


-- | Read-only vectors
class RVector v where
    -- | Get the dimension of the vector.  This is equal to the number of
    -- elements in the vector.                          
    dim :: (Storable e) => v e -> Int

    unsafeSlice :: (Storable e) => Int -> Int -> v e -> v e

    -- | Execute an 'IO' action with a pointer to the first element in the
    -- vector.
    unsafeWith :: (Storable e) => v e -> (Ptr e -> IO a) -> IO a

    -- | Convert a vector to a @ForeignPtr@, and offset, and a length
    unsafeToForeignPtr :: (Storable e)
                       => v e -> (ForeignPtr e, Int, Int)
                             
    -- | Cast a @ForeignPtr@ to a vector.
    unsafeFromForeignPtr :: (Storable e)
                         => ForeignPtr e -- ^ the pointer
                         -> Int          -- ^ the offset
                         -> Int          -- ^ the dimension
                         -> v e


instance RVector (MVector s) where
    dim = STVector.length
    {-# INLINE dim #-}

    unsafeSlice = STVector.unsafeSlice
    {-# INLINE unsafeSlice #-}

    unsafeWith v f =
        STVector.unsafeWith (cast v) f
      where
        cast :: a s e -> a RealWorld e
        cast = unsafeCoerce
    {-# INLINE unsafeWith #-}

    unsafeToForeignPtr = STVector.unsafeToForeignPtr
    {-# INLINE unsafeToForeignPtr #-}
    
    unsafeFromForeignPtr = STVector.unsafeFromForeignPtr
    {-# INLINE unsafeFromForeignPtr #-}


-- | @slice i n v@ creates a subvector view of @v@ starting at
-- index @i@ and having dimension @n@.
slice :: (RVector v, Storable e)
      => Int
      -> Int
      -> v e
      -> v e            
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
drop :: (RVector v, Storable e) => Int -> v e -> v e
drop i v = slice i (dim v - i) v
{-# INLINE drop #-}

-- | @take n v@ is equal to @slice 0 n v@.
take :: (RVector v, Storable e) => Int -> v e -> v e
take n = slice 0 n
{-# INLINE take #-}

-- | Creates a new vector of the given length.  The elements will be
-- uninitialized.
new_ :: (Storable e) => Int -> ST s (STVector s e)
new_ n
    | n < 0 =  error $
        printf "new_ %d: invalid dimension" n
    | otherwise = unsafeIOToST $ do
        f <- mallocForeignPtrArray n
        return $ unsafeFromForeignPtr f 0 n

-- | Create a vector with every element initialized to the same value.
new :: (Storable e) => Int -> e -> ST s (STVector s e)
new n e = do
    x <- new_ n
    setElems x $ replicate n e
    return x

-- | Creates a new vector by copying another one.    
newCopy :: (RVector v, Storable e) => v e -> ST s (STVector s e)
newCopy x = do
    y <- new_ (dim x)
    unsafeCopyTo y x
    return y

-- | @copyTo src dst@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyTo :: (RVector v, Storable e) => STVector s e -> v e -> ST s ()
copyTo = checkOp2 "copyTo" unsafeCopyTo
{-# INLINE copyTo #-}

unsafeCopyTo :: (RVector v, Storable e) => STVector s e -> v e -> ST s ()
unsafeCopyTo dst src = unsafeIOToST $
    unsafeWith dst $ \pdst ->
    unsafeWith src $ \psrc ->
        when (pdst /= psrc) $ copyArray pdst psrc n
  where
    n = dim dst
{-# INLINE unsafeCopyTo #-}

-- | Swap the values stored in two vectors.
swap :: (BLAS1 e) => STVector s e -> STVector s e -> ST s ()
swap = checkOp2 "swap" unsafeSwap
{-# INLINE swap #-}

unsafeSwap :: (BLAS1 e) => STVector s e -> STVector s e -> ST s ()
unsafeSwap = strideCall2 BLAS.swap
{-# INLINE unsafeSwap #-}

-- | Split a vector into two blocks and returns views into the blocks.  If
-- @(x1, x2) = splitAt k x@, then
-- @x1 = slice 0 k x@ and
-- @x2 = slice k (dim x - k) x@.
splitAt :: (RVector v, Storable e) => Int -> v e -> (v e, v e)
splitAt k x
    | k < 0 || k > n = error $
        printf "splitAt %d <vector with dim %d>: invalid index" k n
    | otherwise = let
        x1 = unsafeSlice 0 k     x
        x2 = unsafeSlice k (n-k) x
    in (x1,x2)
  where
    n = dim x
{-# INLINE splitAt #-}

-- | Get the indices of the elements in the vector, @[ 0..n-1 ]@, where
-- @n@ is the dimension of the vector.
indices :: (RVector v, Storable e) => v e -> [Int]
indices x = [ 0..n-1 ] where n = dim x
{-# INLINE indices #-}

-- | Lazily get the elements of the vector.
getElems :: (RVector v, Storable e) => v e -> ST s [e]
getElems v = unsafeIOToST $ unsafeWith v $ \p -> let
    n   = dim v
    end = p `advancePtr` n
    go p' | p' == end = do
              touchVector v
              return []
          | otherwise = unsafeInterleaveIO $ do
              e   <- peek p'
              es  <- go (p' `advancePtr` 1)
              return $ e `seq` (e:es) in
    go p
  where
    touchVector v' = unsafeWith v' $ const (return ())
{-# SPECIALIZE INLINE getElems :: STVector s Double -> ST s [Double] #-}
{-# SPECIALIZE INLINE getElems :: STVector s (Complex Double) -> ST s [Complex Double] #-}

-- | Get the elements of the vector.
getElems' :: (RVector v, Storable e) => v e -> ST s [e]
getElems' v = unsafeIOToST $ unsafeWith v $ \p -> let
    n   = dim v
    end = p `advancePtr` (-1)
    go p' es | p' == end =
                 return es
             | otherwise = do
                 e <- peek p'
                 go (p' `advancePtr` (-1)) (e:es) in
    go (p `advancePtr` (n-1)) []
{-# SPECIALIZE INLINE getElems' :: STVector s Double -> ST s [Double] #-}
{-# SPECIALIZE INLINE getElems' :: STVector s (Complex Double) -> ST s [Complex Double] #-}

-- | Lazily get the association list of the vector.
getAssocs :: (RVector v, Storable e) => v e -> ST s [(Int,e)]
getAssocs x = liftM (zip (indices x)) (getElems x)
{-# INLINE getAssocs #-}

-- | Get the association list of the vector.
getAssocs' :: (RVector v, Storable e) => v e -> ST s [(Int,e)]
getAssocs' x = liftM (zip (indices x)) (getElems' x)
{-# INLINE getAssocs' #-}

-- | Set all of the values of the vector from the elements in the list.
setElems  :: (Storable e) => STVector s e -> [e] -> ST s ()
setElems x es = let
    n = dim x
    go [] i _ | i < n = error $ 
                    printf ("setElems <vector with dim %d>"
                            ++ " <list with length %d>:"
                            ++ " not enough elements") n i
    go (_:_) i _ | i == n = error $ 
                    printf ("setElems <vector with dim %d>"
                            ++ " <list with length at least %d>:"
                            ++ " too many elements") n (i+1)
    go []     _ _ = return ()
    go (f:fs) i p = do
        poke p f
        go fs (i+1) (p `advancePtr` 1)
    
    in unsafeIOToST $ unsafeWith x $ go es 0

-- | Set the given values in the vector.  If an index is repeated twice,
-- the value is implementation-defined.
setAssocs :: (Storable e) => STVector s e -> [(Int,e)] -> ST s ()
setAssocs x ies =
    let n = dim x
        go p ((i,e):ies') = do
            when (i < 0 || i >= n) $ error $
                printf ("setAssocs <vector with dim %d>"
                        ++ " [ ..., (%d,_), ... ]: invalid index") n i
            pokeElemOff p i e
            go p ies'
        go _ [] = return ()
    in unsafeIOToST $ unsafeWith x $ \p -> go p ies

unsafeSetAssocs :: (Storable e) => STVector s e -> [(Int,e)] -> ST s ()
unsafeSetAssocs x ies =
    let go p ((i,e):ies') = do
            pokeElemOff p i e
            go p ies'
        go _ [] = return ()
    in unsafeIOToST $ unsafeWith x $ \p -> go p ies

-- | Get the element stored at the given index.
read :: (RVector v, Storable e) => v e -> Int -> ST s e
read x i
    | i < 0 || i >= n = error $
        printf ("read <vector with dim %d> %d:"
                ++ " invalid index") n i
    | otherwise =
        unsafeRead x i
  where
    n = dim x
{-# SPECIALIZE INLINE read :: STVector s Double -> Int -> ST s (Double) #-}
{-# SPECIALIZE INLINE read :: STVector s (Complex Double) -> Int -> ST s (Complex Double) #-}

unsafeRead :: (RVector v, Storable e) => v e -> Int -> ST s e
unsafeRead x i =
    unsafeIOToST $ unsafeWith x $ \p -> peekElemOff p i
{-# SPECIALIZE INLINE unsafeRead :: STVector s Double -> Int -> ST s (Double) #-}
{-# SPECIALIZE INLINE unsafeRead :: STVector s (Complex Double) -> Int -> ST s (Complex Double) #-}

-- | Set the element stored at the given index.
write :: (Storable e) => STVector s e -> Int -> e -> ST s ()
write x i e
    | i < 0 || i >= n = error $
        printf ("write <vector with dim %d> %d:"
                ++ " invalid index") n i
    | otherwise =
        unsafeWrite x i e
  where
    n = dim x
{-# SPECIALIZE INLINE write :: STVector s Double -> Int -> Double -> ST s () #-}
{-# SPECIALIZE INLINE write :: STVector s (Complex Double) -> Int -> Complex Double -> ST s () #-}

unsafeWrite :: (Storable e) => STVector s e -> Int -> e -> ST s ()
unsafeWrite x i e =
    unsafeIOToST $ unsafeWith x $ \p -> pokeElemOff p i e
{-# SPECIALIZE INLINE unsafeWrite :: STVector s Double -> Int -> Double -> ST s () #-}
{-# SPECIALIZE INLINE unsafeWrite :: STVector s (Complex Double) -> Int -> Complex Double -> ST s () #-}

-- | Modify the element stored at the given index.
modify :: (Storable e) => STVector s e -> Int -> (e -> e) -> ST s ()
modify x i f
    | i < 0 || i >= n = error $
        printf ("modify <vector with dim %d> %d:"
                ++ " invalid index") n i
    | otherwise =
        unsafeModify x i f
  where
    n = dim x
{-# SPECIALIZE INLINE modify :: STVector s Double -> Int -> (Double -> Double) -> ST s () #-}
{-# SPECIALIZE INLINE modify :: STVector s (Complex Double) -> Int -> (Complex Double -> Complex Double) -> ST s () #-}

unsafeModify :: (Storable e) => STVector s e -> Int -> (e -> e) -> ST s ()
unsafeModify x i f =
    unsafeIOToST $ unsafeWith x $ \p -> do
        e <- peekElemOff p i
        pokeElemOff p i $ f e
{-# SPECIALIZE INLINE unsafeModify :: STVector s Double -> Int -> (Double -> Double) -> ST s () #-}
{-# SPECIALIZE INLINE unsafeModify :: STVector s (Complex Double) -> Int -> (Complex Double -> Complex Double) -> ST s () #-}

unsafeSwapElems :: (Storable e) => STVector s e -> Int -> Int -> ST s ()
unsafeSwapElems x i1 i2 = unsafeIOToST $ unsafeWith x $ \p ->
    let p1 = p `advancePtr` i1
        p2 = p `advancePtr` i2
    in do
        e1  <- peek p1
        e2  <- peek p2
        poke p2 e1
        poke p1 e2
{-# SPECIALIZE INLINE unsafeSwapElems :: STVector s Double -> Int -> Int -> ST s () #-}
{-# SPECIALIZE INLINE unsafeSwapElems :: STVector s (Complex Double) -> Int -> Int -> ST s () #-}

-- | @mapTo f x z@ replaces @z@ elementwise with @f(x)@.
mapTo :: (RVector v, Storable e, Storable f)
      => (e -> f)
      -> v e
      -> STVector s f
      -> ST s ()
mapTo f = checkOp2 "mapTo _" $ unsafeMapTo f
{-# INLINE mapTo #-}
                            
unsafeMapTo :: (RVector v, Storable e, Storable f)
            => (e -> f)
            -> v e
            -> STVector s f
            -> ST s ()
unsafeMapTo f v mv =
    let go psrc pdst end
            | pdst == end =
                return ()
            | otherwise = do
                e <- peek psrc
                poke pdst (f e)
                go (psrc `advancePtr` 1) (pdst `advancePtr` 1) end
    in unsafeIOToST $
           unsafeWith v $ \psrc -> 
           unsafeWith mv $ \pdst ->
               go psrc pdst (pdst `advancePtr` dim mv)
  where

{-# INLINE unsafeMapTo #-}

-- | @zipWithTo f x y z@ replaces @z@ elementwise with @f(x, y)@.
zipWithTo :: (RVector v1, RVector v2, Storable e1, Storable e2, Storable f)
          => (e1 -> e2 -> f)
          -> v1 e1
          -> v2 e2
          -> STVector s f
          -> ST s ()
zipWithTo f = checkOp3 "zipWithTo _" $
    unsafeZipWithTo f
{-# INLINE zipWithTo #-}

unsafeZipWithTo :: (RVector v1, RVector v2, Storable e1, Storable e2, Storable f)
                => (e1 -> e2 -> f)
                -> v1 e1
                -> v2 e2
                -> STVector s f
                -> ST s ()
unsafeZipWithTo f v1 v2 mv =
    let go psrc1 psrc2 pdst end
            | pdst == end = 
                return ()
            | otherwise = do
                e1 <- peek psrc1
                e2 <- peek psrc2
                poke pdst (f e1 e2)
                go (psrc1 `advancePtr` 1) (psrc2 `advancePtr` 1)
                   (pdst `advancePtr` 1) end
    in unsafeIOToST $
           unsafeWith v1 $ \psrc1 ->
           unsafeWith v2 $ \psrc2 -> 
           unsafeWith mv $ \pdst ->
               go psrc1 psrc2 pdst (pdst `advancePtr` dim mv)
{-# INLINE unsafeZipWithTo #-}


-- | Set every element in the vector to a default value.  For
-- standard numeric types (including 'Double', 'Complex Double', and 'Int'),
-- the default value is '0'.
clear :: (Storable e) => STVector s e -> ST s ()
clear x = unsafeIOToST $ unsafeWith x $ \p -> clearArray p n
  where
    n = dim x

-- | @negateTo x z@ replaces @z@ with @negate(x)@.
negateTo :: (RVector v, VNum e) => v e -> STVector s e -> ST s ()
negateTo = checkOp2 "negateTo" $
    call2 VMath.vNeg
{-# INLINE negateTo #-}

-- | @absTo x z@ replaces @z@ with @abs(x)@.
absTo :: (RVector v, VNum e) => v e -> STVector s e -> ST s ()
absTo = checkOp2 "absTo" $
    call2 VMath.vAbs
{-# INLINE absTo #-}

-- | @signumTo x z@ replaces @z@ with @signum(x)@.
signumTo :: (RVector v, VNum e) => v e -> STVector s e -> ST s ()
signumTo = checkOp2 "signumTo" $
    call2 VMath.vSgn
{-# INLINE signumTo #-}


-- | @conjugateTo x z@ replaces @z@ with @conjugate(x)@.
conjugateTo :: (RVector v, VNum e) => v e -> STVector s e -> ST s ()
conjugateTo = checkOp2 "conjugateTo" $
    call2 VMath.vConj
{-# INLINE conjugateTo #-}

-- | @shiftTo alpha x z@ replaces @z@ with @alpha + x@.
shiftTo :: (RVector v, VNum e) => e -> v e -> STVector s e -> ST s ()
shiftTo alpha = checkOp2 "shiftTo _" $
    call2 (flip VMath.vShift alpha)
{-# INLINE shiftTo #-}

-- | @scaleTo alpha x z@ replaces @z@ with @alpha * x@.
scaleTo :: (RVector v, VNum e) => e -> v e -> STVector s e -> ST s ()
scaleTo alpha = checkOp2 "scaleTo _" $
    call2 (flip VMath.vScale alpha)
{-# INLINE scaleTo #-}

-- | @addToWithScales alpha x beta y z@ replaces @z@ with
-- @alpha * x + beta * y@.
addToWithScales :: (RVector v1, RVector v2, VNum e)
                => e -> v1 e -> e -> v2 e -> STVector s e -> ST s ()
addToWithScales alpha x beta y z =
    checkOp3 "addToWithScales"
        (call3 (\n v1 v2 v3 -> VMath.vAxpby n alpha v1 beta v2 v3))
        x y z
{-# INLINE addToWithScales #-}

-- | @addTo x y z@ replaces @z@ with @x+y@.
addTo :: (RVector v1, RVector v2, VNum e)
      => v1 e -> v2 e -> STVector s e -> ST s ()
addTo = checkOp3 "addTo" $ call3 VMath.vAdd
{-# INLINE addTo #-}

-- | @subTo x y z@ replaces @z@ with @x-y@.
subTo :: (RVector v1, RVector v2, VNum e)
      => v1 e -> v2 e -> STVector s e -> ST s ()
subTo = checkOp3 "subTo" $ call3 VMath.vSub
{-# INLINE subTo #-}

-- | @mulTo x y z@ replaces @z@ with @x*y@.
mulTo :: (RVector v1, RVector v2, VNum e)
      => v1 e -> v2 e -> STVector s e -> ST s ()
mulTo = checkOp3 "mulTo" $ call3 VMath.vMul
{-# INLINE mulTo #-}
 
-- | @divTo x y z@ replaces @z@ with @x/y@.
divTo :: (RVector v1, RVector v2, VFractional e)
      => v1 e -> v2 e -> STVector s e -> ST s ()
divTo = checkOp3 "divTo" $ call3 VMath.vDiv
{-# INLINE divTo #-}

-- | @recipTo x z@ replaces @z@ with @1/x@.
recipTo :: (RVector v, VFractional e)
        => v e -> STVector s e -> ST s ()
recipTo = checkOp2 "recipTo" $ call2 VMath.vInv
{-# INLINE recipTo #-}

-- | @sqrtTo x z@ replaces @z@ with @sqrt(x)@.
sqrtTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
sqrtTo = checkOp2 "sqrtTo" $
    call2 VMath.vSqrt
{-# INLINE sqrtTo #-}

-- | @expTo x z@ replaces @z@ with @exp(x)@.
expTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
expTo = checkOp2 "expTo" $
    call2 VMath.vExp
{-# INLINE expTo #-}

-- | @logTo x z@ replaces @z@ with @log(x)@.
logTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
logTo = checkOp2 "logTo" $
    call2 VMath.vLog
{-# INLINE logTo #-}

-- | @powTo x y z@ replaces @z@ with @x ** y@.
powTo :: (RVector v1, RVector v2, VFloating e)
            => v1 e -> v2 e -> STVector s e -> ST s ()
powTo = checkOp3 "powTo" $
    call3 VMath.vPow
{-# INLINE powTo #-}

-- | @sinTo x z@ replaces @z@ with @sin(x)@.
sinTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
sinTo = checkOp2 "sinTo" $
    call2 VMath.vSin
{-# INLINE sinTo #-}

-- | @cosTo x z@ replaces @z@ with @cos(x)@.
cosTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
cosTo = checkOp2 "cosTo" $
    call2 VMath.vCos
{-# INLINE cosTo #-}

-- | @tanTo x z@ replaces @z@ with @tan(x)@.
tanTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
tanTo = checkOp2 "tanTo" $
    call2 VMath.vTan
{-# INLINE tanTo #-}

-- | @asinTo x z@ replaces @z@ with @asin(x)@.
asinTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
asinTo = checkOp2 "asinTo" $
    call2 VMath.vASin
{-# INLINE asinTo #-}

-- | @acosTo x z@ replaces @z@ with @acos(x)@.
acosTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
acosTo = checkOp2 "acosTo" $
    call2 VMath.vACos
{-# INLINE acosTo #-}

-- | @atanTo x z@ replaces @z@ with @atan(x)@.
atanTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
atanTo = checkOp2 "atanTo" $
    call2 VMath.vATan
{-# INLINE atanTo #-}

-- | @sinhTo x z@ replaces @z@ with @sinh(x)@.
sinhTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
sinhTo = checkOp2 "sinhTo" $
    call2 VMath.vSinh
{-# INLINE sinhTo #-}

-- | @coshTo x z@ replaces @z@ with @cosh(x)@.
coshTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
coshTo = checkOp2 "coshTo" $
    call2 VMath.vCosh
{-# INLINE coshTo #-}

-- | @tanhTo x z@ replaces @z@ with @tanh(x)@.
tanhTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
tanhTo = checkOp2 "tanhTo" $
    call2 VMath.vTanh
{-# INLINE tanhTo #-}

-- | @asinhTo x z@ replaces @z@ with @asinh(x)@.
asinhTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
asinhTo = checkOp2 "asinhTo" $
    call2 VMath.vASinh
{-# INLINE asinhTo #-}

-- | @acoshTo x z@ replaces @z@ with @acosh(x)@.
acoshTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
acoshTo = checkOp2 "acoshTo" $
    call2 VMath.vACosh
{-# INLINE acoshTo #-}

-- | @atanhTo x z@ replaces @z@ with @atanh(x)@.
atanhTo :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
atanhTo = checkOp2 "atanhTo" $
    call2 VMath.vATanh
{-# INLINE atanhTo #-}


-- | Gets the sum of the absolute values of the vector entries.
getSumAbs :: (RVector v, BLAS1 e) => v e -> ST s Double
getSumAbs = strideCall BLAS.asum
{-# INLINE getSumAbs #-}

-- | Gets the 2-norm of a vector.
getNorm2 :: (RVector v, BLAS1 e) => v e -> ST s Double
getNorm2 = strideCall BLAS.nrm2
{-# INLINE getNorm2 #-}

-- | Gets the index and norm of the element with maximum magnitude.  This is 
-- undefined if any of the elements are @NaN@.  It will throw an exception if 
-- the dimension of the vector is 0.
getWhichMaxAbs :: (RVector v, BLAS1 e) => v e -> ST s (Int, e)
getWhichMaxAbs x =
    case (dim x) of
        0 -> error $ "getWhichMaxAbs <vector with dim 0>: empty vector"
        _ -> do
            i <- strideCall BLAS.iamax x
            e <- unsafeRead x i
            return (i,e)
{-# INLINE getWhichMaxAbs #-}

-- | Computes the dot product of two vectors.
getDot :: (RVector v, RVector v', BLAS1 e)
             => v e -> v' e -> ST s e
getDot = checkOp2 "getDot" unsafeGetDot
{-# INLINE getDot #-}

unsafeGetDot :: (RVector x, RVector y, BLAS1 e)
                   => x e -> y e -> ST s e
unsafeGetDot x y = (strideCall2 BLAS.dotc) y x
{-# INLINE unsafeGetDot #-}

-- | @kroneckerTo x y z@ sets @z := x \otimes y@.
kroneckerTo :: (RVector v1, RVector v2, VNum e)
            => v1 e -> v2 e -> STVector s e -> ST s ()
kroneckerTo x y z
    | dim z /= m * n = error $
        printf ("kroneckerTo"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (dim z)
    | otherwise = do
        ies <- getAssocs x
        forM_ ies $ \(i,e) ->
            scaleTo e y (unsafeSlice (i*n) n z)
  where
    m = dim x
    n = dim y

{-
call :: (RVector x, Storable e)
           => (Int -> Ptr e -> IO a) 
           ->  x e -> ST s a
call f x = 
    let n    = dim x
    in unsafeIOToST $
           unsafeWith x $ \pX ->
               f n pX
{-# INLINE call #-}
-}

call2 :: (RVector x, RVector y, Storable e, Storable f)
      => (Int -> Ptr e -> Ptr f -> IO a) 
      -> x e -> y f -> ST s a
call2 f x y =
    let n    = dim x
    in unsafeIOToST $
           unsafeWith x $ \pX ->
           unsafeWith y $ \pY ->
               f n pX pY
{-# INLINE call2 #-}    

call3 :: (RVector x, RVector y, RVector z, Storable e, Storable f, Storable g)
      => (Int -> Ptr e -> Ptr f -> Ptr g -> IO a) 
      -> x e -> y f -> z g -> ST s a
call3 f x y z =
    let n    = dim x
    in unsafeIOToST $
           unsafeWith x $ \pX ->
           unsafeWith y $ \pY ->
           unsafeWith z $ \pZ ->           
               f n pX pY pZ
{-# INLINE call3 #-}   

strideCall :: (RVector x, Storable e)
           => (Int -> Ptr e -> Int -> IO a) 
           ->  x e -> ST s a
strideCall f x = 
    let n    = dim x
        incX = 1
    in unsafeIOToST $
           unsafeWith x $ \pX ->
               f n pX incX
{-# INLINE strideCall #-}

strideCall2 :: (RVector x, RVector y, Storable e, Storable f)
            => (Int -> Ptr e -> Int -> Ptr f -> Int -> IO a) 
            -> x e -> y f -> ST s a
strideCall2 f x y =
    let n    = dim x
        incX = 1
        incY = 1
    in unsafeIOToST $
           unsafeWith x $ \pX ->
           unsafeWith y $ \pY ->
               f n pX incX pY incY
{-# INLINE strideCall2 #-}    

checkOp2 :: (RVector x, RVector y, Storable e, Storable f)
         => String
         -> (x e -> y f -> a)
         -> x e
         -> y f
         -> a
checkOp2 str f x y
    | n1 /= n2 = error $
        printf ("%s <vector with dim %d> <vector with dim %d>:"
                ++ " dimension mismatch") str n1 n2
    | otherwise =
        f x y
  where
    n1 = dim x
    n2 = dim y        
{-# INLINE checkOp2 #-}

checkOp3 :: (RVector x, RVector y, RVector z, Storable e, Storable f, Storable g)
         => String
         -> (x e -> y f -> z g -> a)
         -> x e
         -> y f
         -> z g
         -> a
checkOp3 str f x y z
    | n1 /= n2 || n1 /= n3 = error $
        printf ("%s <vector with dim %d> <vector with dim %d>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") str n1 n2 n3
    | otherwise =
        f x y z
  where
    n1 = dim x
    n2 = dim y        
    n3 = dim z
{-# INLINE checkOp3 #-}

newResult :: (RVector v, Storable e, Storable f)
          => (v e -> STVector s f -> ST s a)
          -> v e
          -> ST s (STVector s f)
newResult f v = do
    z <- new_ (dim v)
    _ <- f v z
    return z
{-# INLINE newResult #-}


newResult2 :: (RVector v1, RVector v2, Storable e, Storable f, Storable g)
           => (v1 e -> v2 f -> STVector s g -> ST s a)
           -> v1 e
           -> v2 f
           -> ST s (STVector s g)
newResult2 f v1 v2 = do
    z <- new_ (dim v1)
    _ <- f v1 v2 z
    return z
{-# INLINE newResult2 #-}
