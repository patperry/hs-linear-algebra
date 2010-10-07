{-# LANGUAGE Rank2Types #-}
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
    
    stSlice,
    stDrop,
    stTake,
    stSplitAt,
    
    withSliceView,
    withDropView,
    withTakeView,
    withSplitAtView,
    
    withSTSliceView,
    withSTDropView,
    withSTTakeView,
    withSTSplitAtView,
    
    getSumAbs,
    getNorm2,
    getWhichMaxAbs,
    getDot,
    unsafeGetDot,
    scaleByM_,
    addWithScaleM_,
    unsafeAddWithScaleM_,
    kroneckerTo,
    
    addTo,
    subTo,
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
    unsafeWithSliceView :: (Storable e)
                        => Int -> Int -> v e
                        -> (forall v'. RVector v' => v' e -> ST s a)
                        -> ST s a

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

    unsafeWithSliceView i n' v f =
        f (unsafeSlice i n' v)
    {-# INLINE unsafeWithSliceView #-}

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


stSlice :: (Storable e)
        => Int
        -> Int
        -> STVector s e
        -> STVector s e          
stSlice i n' v
    | i < 0 || n' < 0 || i + n' > n = error $
        printf "slice %d %d <vector with dim %d>: index out of range"
               i n' n
    | otherwise =
        unsafeSlice i n' v
  where
    n = dim v
{-# INLINE stSlice #-}

stDrop :: (Storable e) => Int -> STVector s e -> STVector s e
stDrop i v = stSlice i (dim v - i) v
{-# INLINE stDrop #-}

stTake :: (Storable e) => Int -> STVector s e -> STVector s e
stTake n = stSlice 0 n
{-# INLINE stTake #-}

stSplitAt :: (Storable e) => Int -> STVector s e -> (STVector s e, STVector s e)
stSplitAt k x
    | k < 0 || k > n = error $
        printf "splitAt %d <vector with dim %d>: invalid index" k n
    | otherwise = let
        x1 = unsafeSlice 0 k     x
        x2 = unsafeSlice k (n-k) x
    in (x1,x2)
  where
    n = dim x
{-# INLINE stSplitAt #-}

-- | @withSliceView i n v@ performs an action with a view of the
-- @n@-dimensional subvector of @v@ starting at index @i@.
withSliceView :: (RVector v, Storable e)
              => Int
              -> Int
              -> v e
              -> (forall v'. RVector v' => v' e -> ST s a)
              -> ST s a
withSliceView i n' v f
    | i < 0 || n' < 0 || i + n' > n = error $
        printf "withSliceView %d %d <vector with dim %d>: index out of range"
               i n' n
    | otherwise =
        unsafeWithSliceView i n' v f
  where
    n = dim v

withSTSliceView :: (Storable e)
                => Int
                -> Int
                -> STVector s e
                -> (STVector s e -> ST s a)
                -> ST s a
withSTSliceView i n' v f = f $ stSlice i n' v

withDropView :: (RVector v, Storable e)
             => Int
             -> v e
             -> (forall v'. RVector v' => v' e -> ST s a)
             -> ST s a
withDropView i v f =
    withSliceView i (dim v - i) v f

withSTDropView :: (Storable e)
               => Int
               -> STVector s e
               -> (STVector s e -> ST s a)
               -> ST s a
withSTDropView i v f = f $ stDrop i v
    
withTakeView :: (RVector v, Storable e)
             => Int
             -> v e
             -> (forall v'. RVector v' => v' e -> ST s a)
             -> ST s a
withTakeView n v f =
    withSliceView 0 n v f

withSTTakeView :: (Storable e)
               => Int
               -> STVector s e
               -> (STVector s e -> ST s a)
               -> ST s a
withSTTakeView n v f = f $ stTake n v

withSplitAtView :: (RVector v, Storable e)
                => Int
                -> v e
                -> (forall v1' v2'. (RVector v1', RVector v2') => v1' e -> v2' e -> ST s a)
                -> ST s a
withSplitAtView i v f =
    withSliceView 0 i v $ \v1 ->
    withSliceView i (dim v - i) v $ \v2 ->
        f v1 v2

withSTSplitAtView :: (Storable e)
                  => Int
                  -> STVector s e
                  -> (STVector s e -> STVector s e -> ST s a)
                  -> ST s a
withSTSplitAtView i v f = 
    let (v1,v2) = stSplitAt i v in f v1 v2

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

-- | @copyTo dst src@ replaces the values in @dst@ with those in
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

-- | @mapTo dst f src@ replaces @dst@ elementwise with @f(src)@.
mapTo :: (RVector v, Storable e, Storable f)
      => STVector s f
      -> (e -> f)
      -> v e
      -> ST s ()
mapTo dst f src = (checkOp2 "mapTo _" $ \z x -> unsafeMapTo z f x) dst src
{-# INLINE mapTo #-}
                            
unsafeMapTo :: (RVector v, Storable e, Storable f)
            => STVector s f
            -> (e -> f)
            -> v e
            -> ST s ()
unsafeMapTo dst f src =
    let go end pdst psrc
            | pdst == end =
                return ()
            | otherwise = do
                e <- peek psrc
                poke pdst (f e)
                go end (pdst `advancePtr` 1) (psrc `advancePtr` 1)
    in unsafeIOToST $
           unsafeWith dst $ \pdst ->
           unsafeWith src $ \psrc -> 
               go (pdst `advancePtr` dim dst) pdst psrc
  where

{-# INLINE unsafeMapTo #-}

-- | @zipWithTo dst f x y@ replaces @dst@ elementwise with @f(x,y)@.
zipWithTo :: (RVector v1, RVector v2, Storable e1, Storable e2, Storable f)
          => STVector s f
          -> (e1 -> e2 -> f)
          -> v1 e1
          -> v2 e2
          -> ST s ()
zipWithTo dst f x y = 
    (checkOp3 "zipWithTo _" $ \dst1 x1 y1 -> unsafeZipWithTo dst1 f x1 y1)
        dst x y
{-# INLINE zipWithTo #-}

unsafeZipWithTo :: (RVector v1, RVector v2, Storable e1, Storable e2, Storable f)
                => STVector s f
                -> (e1 -> e2 -> f)
                -> v1 e1
                -> v2 e2
                -> ST s ()
unsafeZipWithTo dst f src1 src2 =
    let go end pdst psrc1 psrc2 
            | pdst == end = 
                return ()
            | otherwise = do
                e1 <- peek psrc1
                e2 <- peek psrc2
                poke pdst (f e1 e2)
                go end (pdst `advancePtr` 1) (psrc1 `advancePtr` 1)
                   (psrc2 `advancePtr` 1)
    in unsafeIOToST $
           unsafeWith dst $ \pdst ->
           unsafeWith src1 $ \psrc1 ->
           unsafeWith src2 $ \psrc2 -> 
               go (pdst `advancePtr` dim dst) pdst psrc1 psrc2
{-# INLINE unsafeZipWithTo #-}


-- | Set every element in the vector to a default value.  For
-- standard numeric types (including 'Double', 'Complex Double', and 'Int'),
-- the default value is '0'.
clear :: (Storable e) => STVector s e -> ST s ()
clear x = unsafeIOToST $ unsafeWith x $ \p -> clearArray p n
  where
    n = dim x

-- | @negateTo dst x@ replaces @dst@ with @negate(x)@.
negateTo :: (RVector v, VNum e) => STVector s e -> v e -> ST s ()
negateTo = checkOp2 "negateTo" $ \dst x -> 
    call2 VMath.vNeg x dst
{-# INLINE negateTo #-}

-- | @absTo dst x@ replaces @dst@ with @abs(x)@.
absTo :: (RVector v, VNum e) => STVector s e -> v e -> ST s ()
absTo = checkOp2 "absTo" $ \dst x -> 
    call2 VMath.vAbs x dst
{-# INLINE absTo #-}

-- | @signumTo dst x@ replaces @dst@ with @signum(x)@.
signumTo :: (RVector v, VNum e) => STVector s e -> v e -> ST s ()
signumTo = checkOp2 "signumTo" $ \dst x ->
    call2 VMath.vSgn x dst
{-# INLINE signumTo #-}


-- | @conjugateTo dst x@ replaces @dst@ with @conjugate(x)@.
conjugateTo :: (RVector v, VNum e) => STVector s e -> v e -> ST s ()
conjugateTo = checkOp2 "conjugateTo" $ \dst x ->
    call2 VMath.vConj x dst
{-# INLINE conjugateTo #-}

-- | @addTo dst x y@ replaces @dst@ with @x+y@.
addTo :: (RVector v1, RVector v2, VNum e)
      =>  STVector s e -> v1 e -> v2 e -> ST s ()
addTo = checkOp3 "addTo" $ \dst x y -> call3 VMath.vAdd x y dst
{-# INLINE addTo #-}

-- | @subTo dst x y@ replaces @dst@ with @x-y@.
subTo :: (RVector v1, RVector v2, VNum e)
      => STVector s e -> v1 e -> v2 e -> ST s ()
subTo = checkOp3 "subTo" $ \dst x y -> call3 VMath.vSub x y dst
{-# INLINE subTo #-}

-- | @mulTo dst x y@ replaces @dst@ with @x*y@.
mulTo :: (RVector v1, RVector v2, VNum e)
      => STVector s e -> v1 e -> v2 e -> ST s ()
mulTo = checkOp3 "mulTo" $ \dst x y -> call3 VMath.vMul x y dst
{-# INLINE mulTo #-}
 
-- | @divTo dst x y@ replaces @dst@ with @x/y@.
divTo :: (RVector v1, RVector v2, VFractional e)
      => STVector s e -> v1 e -> v2 e -> ST s ()
divTo = checkOp3 "divTo" $ \dst x y -> call3 VMath.vDiv x y dst
{-# INLINE divTo #-}

-- | @recipTo dst x@ replaces @dst@ with @1/x@.
recipTo :: (RVector v, VFractional e)
        => STVector s e -> v e -> ST s ()
recipTo = checkOp2 "recipTo" $ \dst x -> call2 VMath.vInv x dst
{-# INLINE recipTo #-}

-- | @sqrtTo dst x@ replaces @dst@ with @sqrt(x)@.
sqrtTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
sqrtTo = checkOp2 "sqrtTo" $ \dst x -> call2 VMath.vSqrt x dst
{-# INLINE sqrtTo #-}

-- | @expTo dst x@ replaces @dst@ with @exp(x)@.
expTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
expTo = checkOp2 "expTo" $ \dst x -> call2 VMath.vExp x dst
{-# INLINE expTo #-}

-- | @logTo dst x@ replaces @dst@ with @log(x)@.
logTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
logTo = checkOp2 "logTo" $ \dst x -> call2 VMath.vLog x dst
{-# INLINE logTo #-}

-- | @powTo dst x y@ replaces @dst@ with @x ** y@.
powTo :: (RVector v1, RVector v2, VFloating e)
      => STVector s e -> v1 e -> v2 e -> ST s ()
powTo = checkOp3 "powTo" $ \dst x y -> call3 VMath.vPow x y dst
{-# INLINE powTo #-}

-- | @sinTo dst x@ replaces @dst@ with @sin(x)@.
sinTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
sinTo = checkOp2 "sinTo" $ \dst x -> call2 VMath.vSin x dst
{-# INLINE sinTo #-}

-- | @cosTo dst x@ replaces @dst@ with @cos(x)@.
cosTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
cosTo = checkOp2 "cosTo" $ \dst x -> call2 VMath.vCos x dst
{-# INLINE cosTo #-}

-- | @tanTo dst x@ replaces @dst@ with @tan(x)@.
tanTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
tanTo = checkOp2 "tanTo" $ \dst x -> call2 VMath.vTan x dst
{-# INLINE tanTo #-}

-- | @asinTo dst x@ replaces @dst@ with @asin(x)@.
asinTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
asinTo = checkOp2 "asinTo" $ \dst x -> call2 VMath.vASin x dst
{-# INLINE asinTo #-}

-- | @acosTo dst x@ replaces @dst@ with @acos(x)@.
acosTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
acosTo = checkOp2 "acosTo" $ \dst x -> call2 VMath.vACos x dst
{-# INLINE acosTo #-}

-- | @atanTo dst x@ replaces @dst@ with @atan(x)@.
atanTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
atanTo = checkOp2 "atanTo" $ \dst x -> call2 VMath.vATan x dst
{-# INLINE atanTo #-}

-- | @sinhTo dst x@ replaces @dst@ with @sinh(x)@.
sinhTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
sinhTo = checkOp2 "sinhTo" $ \dst x -> call2 VMath.vSinh x dst
{-# INLINE sinhTo #-}

-- | @coshTo dst x@ replaces @dst@ with @cosh(x)@.
coshTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
coshTo = checkOp2 "coshTo" $ \dst x -> call2 VMath.vCosh x dst
{-# INLINE coshTo #-}

-- | @tanhTo dst x@ replaces @dst@ with @tanh(x)@.
tanhTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
tanhTo = checkOp2 "tanhTo" $ \dst x -> call2 VMath.vTanh x dst
{-# INLINE tanhTo #-}

-- | @asinhTo dst x@ replaces @dst@ with @asinh(x)@.
asinhTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
asinhTo = checkOp2 "asinhTo" $ \dst x -> call2 VMath.vASinh x dst
{-# INLINE asinhTo #-}

-- | @acoshTo dst x@ replaces @dst@ with @acosh(x)@.
acoshTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
acoshTo = checkOp2 "acoshTo" $ \dst x -> call2 VMath.vACosh x dst
{-# INLINE acoshTo #-}

-- | @atanhTo dst x@ replaces @dst@ with @atanh(x)@.
atanhTo :: (RVector v, VFloating e) => STVector s e -> v e -> ST s ()
atanhTo = checkOp2 "atanhTo" $ \dst x -> call2 VMath.vATanh x dst
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

-- | @scaleByM k x@ sets @x := k * x@.
scaleByM_ :: (Storable e, BLAS1 e) => e -> STVector s e -> ST s ()
scaleByM_ k x =
    unsafeIOToST $
        unsafeWith x $ \px ->
            BLAS.scal (dim x) k px 1
{-# INLINE scaleByM_ #-}

-- | @addWithScaleM_ alpha x y@ sets @y := alpha * x + y@.
addWithScaleM_ :: (RVector v, BLAS1 e) => e -> v e -> STVector s e -> ST s ()
addWithScaleM_ alpha x y =
    (checkOp2 "addWithScaleM_" $ \x1 y1 -> unsafeAddWithScaleM_ alpha x1 y1)
        x y
{-# INLINE addWithScaleM_ #-}

unsafeAddWithScaleM_ :: (RVector v, BLAS1 e)
                     => e -> v e -> STVector s e -> ST s ()
unsafeAddWithScaleM_ alpha x y =
    (strideCall2 $ flip BLAS.axpy alpha) x y
{-# INLINE unsafeAddWithScaleM_ #-}

-- | @kroneckerTo dst x y@ sets @dst := x \otimes y@.
kroneckerTo :: (RVector v1, RVector v2, BLAS2 e)
            => STVector s e -> v1 e -> v2 e -> ST s ()
kroneckerTo dst x y
    | dim dst /= m * n = error $
        printf ("kroneckerTo"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") (dim dst) m n 
    | otherwise = do
        clear dst
        unsafeIOToST $
            unsafeWith dst $ \pdst ->
            unsafeWith x $ \px ->
            unsafeWith y $ \py ->
                BLAS.geru n m 1 py 1 px 1 pdst (max n 1)
  where
    m = dim x
    n = dim y

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
          => (STVector s f -> v e -> ST s a)
          -> v e
          -> ST s (STVector s f)
newResult f v = do
    z <- new_ (dim v)
    _ <- f z v
    return z
{-# INLINE newResult #-}


newResult2 :: (RVector v1, RVector v2, Storable e, Storable f, Storable g)
           => (STVector s g -> v1 e -> v2 f -> ST s a)
           -> v1 e
           -> v2 f
           -> ST s (STVector s g)
newResult2 f v1 v2 = do
    z <- new_ (dim v1)
    _ <- f z v1 v2
    return z
{-# INLINE newResult2 #-}
