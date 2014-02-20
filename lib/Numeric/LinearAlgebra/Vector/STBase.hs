{-# LANGUAGE Rank2Types, DeriveDataTypeable #-}
{-# OPTIONS_HADDOCK hide #-}
{-# OPTIONS_GHC -fno-warn-unused-imports #-}
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
    
    create,
    freeze,
    
    new_,
    new,
    newCopy,
    
    clear,
    copyTo,
    unsafeCopyTo,
    swap,
    unsafeSwap,
    
    getIndices,
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
    
    withSlice,
    withDrop,
    withTake,
    withSplitAt,
    
    withSliceM,
    withDropM,
    withTakeM,
    withSplitAtM,
    
    getSumAbs,
    getNorm2,
    getWhichMaxAbs,
    getDot,
    unsafeGetDot,
    scaleM_,
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

import Control.Monad( when, liftM2 )
import Control.Monad.ST( RealWorld, ST, runST, stToIO )
import Control.Monad.ST.Unsafe( unsafeIOToST )
import Data.Complex( Complex )
import Data.Typeable( Typeable )
import Foreign( Ptr, Storable, advancePtr, peek, poke, peekElemOff,
    pokeElemOff, copyArray, mallocForeignPtrArray )
import System.IO.Unsafe( unsafeInterleaveIO )
import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import Foreign.BLAS( BLAS1, BLAS2 )
import qualified Foreign.BLAS as BLAS
import Foreign.VMath( VNum, VFractional, VFloating )
import qualified Foreign.VMath as VMath

import Numeric.LinearAlgebra.Internal( clearArray )
import Numeric.LinearAlgebra.Vector.Base hiding ( unsafeWith )
import qualified Numeric.LinearAlgebra.Vector.Base as V


-- | Mutable vectors in the 'ST' monad.
newtype STVector s e = STVector { unSTVector :: Vector e }
    deriving (Typeable)

-- | Mutable vectors in the 'IO' monad.  Note that 'IO' operations
-- aren't directly supported; to perform an operation in the 'IO'
-- monad, perform the action in 'ST' 'RealWorld' and then convert
-- it via 'stToIO'.
type IOVector s = STVector RealWorld

-- | A safe way to create and work with a mutable vector before returning 
-- an immutable vector for later perusal. This function avoids copying
-- the vector before returning it - it uses 'unsafeFreeze' internally,
-- but this wrapper is a safe interface to that function.
create :: (Storable e) => (forall s . ST s (STVector s e)) -> Vector e
create stmv = runST $ do
    mv <- stmv
    unsafeFreeze mv
{-# INLINE create #-}

-- | Converts a mutable vector to an immutable one by taking a complete
-- copy of it.
freeze :: (RVector v, Storable e) => v e -> ST s (Vector e)
freeze mv = do
    mv' <- newCopy mv
    unsafeFreeze mv'
{-# INLINE freeze #-}



-- | Read-only vectors
class RVector v where
    -- | Get the dimension of the vector.  This is equal to the number of
    -- elements in the vector.                          
    getDim :: (Storable e) => v e -> ST s Int

    -- | Same as 'withSlice' but does not range-check indices.
    unsafeWithSlice :: (Storable e)
                    => Int -> Int -> v e
                    -> (forall v'. RVector v' => v' e -> ST s a)
                    -> ST s a

    -- | Execute an 'IO' action with a pointer to the first element in the
    -- vector.
    unsafeWith :: (Storable e) => v e -> (Ptr e -> IO a) -> IO a

    -- | Converts a read-only vector into an immutable vector. This simply
    -- casts the vector from one type to the other without copying the vector.
    -- Note that because the vector is possibly not copied, any subsequent
    -- modifications made to the mutable version of the vector may be shared
    -- with the immutable version. It is safe to use, therefore, if the
    -- mutable version is never modified after the freeze operation.
    unsafeFreeze :: (Storable e) => v e -> ST s (Vector e)

    -- | Unsafe cast from a read-only vector to a mutable vector.
    unsafeThaw :: (Storable e)
               => v e -> ST s (STVector s e)


instance RVector Vector where
    getDim = return . dim
    {-# INLINE getDim #-}

    unsafeWith = V.unsafeWith
    {-# INLINE unsafeWith #-}

    unsafeWithSlice i n' v f =
        f (V.unsafeSlice i n' v)
    {-# INLINE unsafeWithSlice #-}

    unsafeFreeze = return . id
    {-# INLINE unsafeFreeze #-}

    unsafeThaw = return . STVector
    {-# INLINE unsafeThaw #-}


instance RVector (STVector s) where
    getDim = return . dim . unSTVector
    {-# INLINE getDim #-}

    unsafeWith v f = V.unsafeWith (unSTVector v) f
    {-# INLINE unsafeWith #-}

    unsafeWithSlice i n' v f =
        f $ unsafeSlice i n' (unSTVector v)
    {-# INLINE unsafeWithSlice #-}

    unsafeFreeze = return . unSTVector
    {-# INLINE unsafeFreeze #-}

    unsafeThaw v = return $ cast v
      where
        cast :: STVector s e -> STVector s' e
        cast = unsafeCoerce
    {-# INLINE unsafeThaw #-}


-- | @withSlice i n v@ performs an action with a view of the
-- @n@-dimensional subvector of @v@ starting at index @i@.
withSlice :: (RVector v, Storable e)
          => Int
          -> Int
          -> v e
          -> (forall v'. RVector v' => v' e -> ST s a)
          -> ST s a
withSlice i n' v f = do
    n <- getDim v
    when (i < 0 || n' < 0 || i + n' > n) $ error $
        printf "withSlice %d %d <vector with dim %d>: index out of range"
               i n' n
    unsafeWithSlice i n' v f

-- | Same as 'withSliceM' but does not range-check indices.
unsafeWithSliceM :: (Storable e)
                 => Int
                 -> Int
                 -> STVector s e
                 -> (STVector s e -> ST s a)
                 -> ST s a
unsafeWithSliceM i n' v f =
    f $ STVector $ unsafeSlice i n' (unSTVector v)

-- | Like 'withSlice', but perform the action with a mutable view 
-- of the vector.
withSliceM :: (Storable e)
           => Int
           -> Int
           -> STVector s e
           -> (STVector s e -> ST s a)
           -> ST s a
withSliceM i n' v f = do
    n <- getDim v
    when (i < 0 || n' < 0 || i + n' > n) $ error $
        printf "withSlice %d %d <vector with dim %d>: index out of range"
               i n' n
    unsafeWithSliceM i n' v f

-- | Like 'withDrop', but perform the action with a mutable view 
-- of the vector.
withDropM :: (Storable e)
          => Int
          -> STVector s e
          -> (STVector s e -> ST s a)
          -> ST s a
withDropM i v f = do
    n <- getDim v
    withSliceM i (n-i) v f

-- | Like 'withTake', but perform the action with a mutable view 
-- of the vector.
withTakeM :: (Storable e)
          => Int
          -> STVector s e
          -> (STVector s e -> ST s a)
          -> ST s a
withTakeM = withSliceM 0

-- | Like 'withSplitAt' but perform the action with mutable views
-- of the vector.
withSplitAtM :: (Storable e)
             => Int
             -> STVector s e
             -> (STVector s e -> STVector s e -> ST s a)
             -> ST s a
withSplitAtM i v f = do
    n <- getDim v
    withSliceM 0 i v $ \v1 ->
        withSliceM i (n-i) v $ \v2 ->
            f v1 v2

-- | Perform an action the a view gotten from dropping the given
-- number of elements from the start of the vector.
withDrop :: (RVector v, Storable e)
         => Int
         -> v e
         -> (forall v'. RVector v' => v' e -> ST s a)
         -> ST s a
withDrop i v f = do
    mv <- unsafeThaw v
    withDropM i mv f

-- | Perform an action with a view gotten from taking the given
-- number of elements from the start of the vector.
withTake :: (RVector v, Storable e)
         => Int
         -> v e
         -> (forall v'. RVector v' => v' e -> ST s a)
         -> ST s a
withTake n v f = do
    mv <- unsafeThaw v
    withTakeM n mv f

-- | Perform an action with views from splitting the vector at the
-- given index.
withSplitAt :: (RVector v, Storable e)
            => Int
            -> v e
            -> (forall v1' v2'. (RVector v1', RVector v2') => v1' e -> v2' e -> ST s a)
            -> ST s a
withSplitAt i v f = do
    mv <- unsafeThaw v
    withSplitAtM i mv f

-- | Creates a new vector of the given length.  The elements will be
-- uninitialized.
new_ :: (Storable e) => Int -> ST s (STVector s e)
new_ n
    | n < 0 =  error $
        printf "new_ %d: invalid dimension" n
    | otherwise = unsafeIOToST $ do
        f <- mallocForeignPtrArray n
        return $ STVector $ V.unsafeFromForeignPtr f 0 n

-- | Create a vector with every element initialized to the same value.
new :: (Storable e) => Int -> e -> ST s (STVector s e)
new n e = do
    x <- new_ n
    setElems x $ replicate n e
    return x

-- | Creates a new vector by copying another one.    
newCopy :: (RVector v, Storable e) => v e -> ST s (STVector s e)
newCopy x = do
    n <- getDim x
    y <- new_ n
    unsafeCopyTo y x
    return y

-- | @copyTo dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyTo :: (RVector v, Storable e) => STVector s e -> v e -> ST s ()
copyTo = checkOp2 "copyTo" unsafeCopyTo
{-# INLINE copyTo #-}

-- | Same as 'copyTo' but does not check the dimensions.
unsafeCopyTo :: (RVector v, Storable e) => STVector s e -> v e -> ST s ()
unsafeCopyTo dst src = do
    n <- getDim dst
    unsafeIOToST $
        unsafeWith dst $ \pdst ->
        unsafeWith src $ \psrc ->
            copyArray pdst psrc n
{-# INLINE unsafeCopyTo #-}

-- | Swap the values stored in two vectors.
swap :: (BLAS1 e) => STVector s e -> STVector s e -> ST s ()
swap = checkOp2 "swap" unsafeSwap
{-# INLINE swap #-}

-- | Same as 'swap' but does not check the dimensions.
unsafeSwap :: (BLAS1 e) => STVector s e -> STVector s e -> ST s ()
unsafeSwap = strideCall2 BLAS.swap
{-# INLINE unsafeSwap #-}

-- | Get the indices of the elements in the vector, @[ 0..n-1 ]@, where
-- @n@ is the dimension of the vector.
getIndices :: (RVector v, Storable e) => v e -> ST s [Int]
getIndices v = do
    n <- getDim v
    return $ [ 0..n-1 ]
{-# INLINE getIndices #-}

-- | Lazily get the elements of the vector.
getElems :: (RVector v, Storable e) => v e -> ST s [e]
getElems v = let
    go end p' | p' == end = do
                  touch v
                  return []
              | otherwise = unsafeInterleaveIO $ do
                  e   <- peek p'
                  es  <- go end (p' `advancePtr` 1)
                  return $ e `seq` (e:es)
    in do
        n <- getDim v
        unsafeIOToST $
            unsafeWith v $ \p -> 
                go (p `advancePtr` n) p
  where
    touch v' = unsafeWith v' $ const (return ())
{-# SPECIALIZE INLINE getElems :: STVector s Double -> ST s [Double] #-}
{-# SPECIALIZE INLINE getElems :: STVector s (Complex Double) -> ST s [Complex Double] #-}

-- | Get the elements of the vector.
getElems' :: (RVector v, Storable e) => v e -> ST s [e]
getElems' v = let
    go end p' es | p' == end =
                     return es
                 | otherwise = do
                     e <- peek p'
                     go end (p' `advancePtr` (-1)) (e:es)
    in do
        n <- getDim v
        unsafeIOToST $
            unsafeWith v $ \p ->
                go (p `advancePtr` (-1)) (p `advancePtr` (n-1)) []
{-# SPECIALIZE INLINE getElems' :: STVector s Double -> ST s [Double] #-}
{-# SPECIALIZE INLINE getElems' :: STVector s (Complex Double) -> ST s [Complex Double] #-}

-- | Lazily get the association list of the vector.
getAssocs :: (RVector v, Storable e) => v e -> ST s [(Int,e)]
getAssocs x = liftM2 zip (getIndices x) (getElems x)
{-# INLINE getAssocs #-}

-- | Get the association list of the vector.
getAssocs' :: (RVector v, Storable e) => v e -> ST s [(Int,e)]
getAssocs' x = liftM2 zip (getIndices x) (getElems' x)
{-# INLINE getAssocs' #-}

-- | Set all of the values of the vector from the elements in the list.
setElems  :: (Storable e) => STVector s e -> [e] -> ST s ()
setElems x es = let
    go n [] i _ | i < n = error $ 
                    printf ("setElems <vector with dim %d>"
                            ++ " <list with length %d>:"
                            ++ " not enough elements") n i
    go n (_:_) i _ | i == n = error $ 
                    printf ("setElems <vector with dim %d>"
                            ++ " <list with length at least %d>:"
                            ++ " too many elements") n (i+1)
    go _ []     _ _ = return ()
    go n (f:fs) i p = do
        poke p f
        go n fs (i+1) (p `advancePtr` 1)
    
    in do
        n <- getDim x
        unsafeIOToST $ unsafeWith x $ go n es 0

-- | Set the given values in the vector.  If an index is repeated twice,
-- the value is implementation-defined.
setAssocs :: (Storable e) => STVector s e -> [(Int,e)] -> ST s ()
setAssocs x ies =
    let go n p ((i,e):ies') = do
            when (i < 0 || i >= n) $ error $
                printf ("setAssocs <vector with dim %d>"
                        ++ " [ ..., (%d,_), ... ]: invalid index") n i
            pokeElemOff p i e
            go n p ies'
        go _ _ [] = return ()
    in do
        n <- getDim x
        unsafeIOToST $ unsafeWith x $ \p -> go n p ies

unsafeSetAssocs :: (Storable e) => STVector s e -> [(Int,e)] -> ST s ()
unsafeSetAssocs x ies =
    let go p ((i,e):ies') = do
            pokeElemOff p i e
            go p ies'
        go _ [] = return ()
    in unsafeIOToST $ unsafeWith x $ \p -> go p ies

-- | Get the element stored at the given index.
read :: (RVector v, Storable e) => v e -> Int -> ST s e
read x i = do
    n <- getDim x
    when (i < 0 || i >= n) $ error $
        printf ("read <vector with dim %d> %d:"
                ++ " invalid index") n i
    unsafeRead x i
{-# SPECIALIZE INLINE read :: STVector s Double -> Int -> ST s (Double) #-}
{-# SPECIALIZE INLINE read :: STVector s (Complex Double) -> Int -> ST s (Complex Double) #-}

-- | Same as 'read' but does not range check the index.
unsafeRead :: (RVector v, Storable e) => v e -> Int -> ST s e
unsafeRead x i =
    unsafeIOToST $ unsafeWith x $ \p -> peekElemOff p i
{-# SPECIALIZE INLINE unsafeRead :: STVector s Double -> Int -> ST s (Double) #-}
{-# SPECIALIZE INLINE unsafeRead :: STVector s (Complex Double) -> Int -> ST s (Complex Double) #-}

-- | Set the element stored at the given index.
write :: (Storable e) => STVector s e -> Int -> e -> ST s ()
write x i e = do
    n <- getDim x
    when (i < 0 || i >= n) $ error $
        printf ("write <vector with dim %d> %d:"
                ++ " invalid index") n i
    unsafeWrite x i e
{-# SPECIALIZE INLINE write :: STVector s Double -> Int -> Double -> ST s () #-}
{-# SPECIALIZE INLINE write :: STVector s (Complex Double) -> Int -> Complex Double -> ST s () #-}

-- | Same as 'write' but does not range check the index.
unsafeWrite :: (Storable e) => STVector s e -> Int -> e -> ST s ()
unsafeWrite x i e =
    unsafeIOToST $ unsafeWith x $ \p -> pokeElemOff p i e
{-# SPECIALIZE INLINE unsafeWrite :: STVector s Double -> Int -> Double -> ST s () #-}
{-# SPECIALIZE INLINE unsafeWrite :: STVector s (Complex Double) -> Int -> Complex Double -> ST s () #-}

-- | Modify the element stored at the given index.
modify :: (Storable e) => STVector s e -> Int -> (e -> e) -> ST s ()
modify x i f = do
    n <- getDim x
    when (i < 0 || i >= n) $ error $
        printf ("modify <vector with dim %d> %d:"
                ++ " invalid index") n i
    unsafeModify x i f
{-# SPECIALIZE INLINE modify :: STVector s Double -> Int -> (Double -> Double) -> ST s () #-}
{-# SPECIALIZE INLINE modify :: STVector s (Complex Double) -> Int -> (Complex Double -> Complex Double) -> ST s () #-}

-- | Same as 'modify' but does not range check the index.
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

-- | Same as 'mapTo' but does not check dimensions.
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
    in do
        ndst <- getDim dst
        unsafeIOToST $
           unsafeWith dst $ \pdst ->
           unsafeWith src $ \psrc -> 
               go (pdst `advancePtr` ndst) pdst psrc
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

-- | Same as 'zipWithTo' but does not range-check dimensions.
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
    in do
        ndst <- getDim dst
        unsafeIOToST $
           unsafeWith dst $ \pdst ->
           unsafeWith src1 $ \psrc1 ->
           unsafeWith src2 $ \psrc2 -> 
               go (pdst `advancePtr` ndst) pdst psrc1 psrc2
{-# INLINE unsafeZipWithTo #-}


-- | Set every element in the vector to a default value.  For
-- standard numeric types (including 'Double', 'Complex Double', and 'Int'),
-- the default value is '0'.
clear :: (Storable e) => STVector s e -> ST s ()
clear x = do
    n <- getDim x
    unsafeIOToST $ unsafeWith x $ \p -> clearArray p n

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
getWhichMaxAbs x = do
    n <- getDim x
    when (n == 0) $ error $
        "getWhichMaxAbs <vector with dim 0>: empty vector"

    i <- strideCall BLAS.iamax x
    e <- unsafeRead x i
    return (i,e)
{-# INLINE getWhichMaxAbs #-}

-- | Computes the dot product of two vectors.
getDot :: (RVector v, RVector v', BLAS1 e)
       => v e -> v' e -> ST s e
getDot = checkOp2 "getDot" unsafeGetDot
{-# INLINE getDot #-}

-- | Same as 'getDot' but does not check dimensions.
unsafeGetDot :: (RVector x, RVector y, BLAS1 e)
             => x e -> y e -> ST s e
unsafeGetDot x y = (strideCall2 BLAS.dotc) y x
{-# INLINE unsafeGetDot #-}

-- | @scaleM k x@ sets @x := k * x@.
scaleM_ :: (Storable e, BLAS1 e) => e -> STVector s e -> ST s ()
scaleM_ k x = do
    n <- getDim x
    unsafeIOToST $
        unsafeWith x $ \px ->
            BLAS.scal n k px 1
{-# INLINE scaleM_ #-}

-- | @addWithScaleM_ alpha x y@ sets @y := alpha * x + y@.
addWithScaleM_ :: (RVector v, BLAS1 e) => e -> v e -> STVector s e -> ST s ()
addWithScaleM_ alpha x y =
    (checkOp2 "addWithScaleM_" $ \x1 y1 -> unsafeAddWithScaleM_ alpha x1 y1)
        x y
{-# INLINE addWithScaleM_ #-}

-- | Same as 'addWithScaleM_' but does not check dimensions.
unsafeAddWithScaleM_ :: (RVector v, BLAS1 e)
                     => e -> v e -> STVector s e -> ST s ()
unsafeAddWithScaleM_ alpha x y =
    (strideCall2 $ flip BLAS.axpy alpha) x y
{-# INLINE unsafeAddWithScaleM_ #-}

-- | @kroneckerTo dst x y@ sets @dst := x \otimes y@.
kroneckerTo :: (RVector v1, RVector v2, BLAS2 e)
            => STVector s e -> v1 e -> v2 e -> ST s ()
kroneckerTo dst x y = do
    m <- getDim x
    n <- getDim y
    dimdst <- getDim dst
    
    when (dimdst /= m * n) $ error $
        printf ("kroneckerTo"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") dimdst m n 

    clear dst
    unsafeIOToST $
        unsafeWith dst $ \pdst ->
        unsafeWith x $ \px ->
        unsafeWith y $ \py ->
            BLAS.geru n m 1 py 1 px 1 pdst (max n 1)

call2 :: (RVector x, RVector y, Storable e, Storable f)
      => (Int -> Ptr e -> Ptr f -> IO a) 
      -> x e -> y f -> ST s a
call2 f x y = do
    n <- getDim x
    unsafeIOToST $
        unsafeWith x $ \pX ->
        unsafeWith y $ \pY ->
            f n pX pY
{-# INLINE call2 #-}    

call3 :: (RVector x, RVector y, RVector z, Storable e, Storable f, Storable g)
      => (Int -> Ptr e -> Ptr f -> Ptr g -> IO a) 
      -> x e -> y f -> z g -> ST s a
call3 f x y z = do
    n <- getDim x
    unsafeIOToST $
        unsafeWith x $ \pX ->
        unsafeWith y $ \pY ->
        unsafeWith z $ \pZ ->           
            f n pX pY pZ
{-# INLINE call3 #-}   

strideCall :: (RVector x, Storable e)
           => (Int -> Ptr e -> Int -> IO a) 
           ->  x e -> ST s a
strideCall f x = do
    n <- getDim x
    unsafeIOToST $
        unsafeWith x $ \pX ->
            f n pX incX
  where
    incX = 1
{-# INLINE strideCall #-}

strideCall2 :: (RVector x, RVector y, Storable e, Storable f)
            => (Int -> Ptr e -> Int -> Ptr f -> Int -> IO a) 
            -> x e -> y f -> ST s a
strideCall2 f x y = do
    n <- getDim x
    unsafeIOToST $
       unsafeWith x $ \pX ->
       unsafeWith y $ \pY ->
           f n pX incX pY incY
  where
    incX = 1
    incY = 1
{-# INLINE strideCall2 #-}    

checkOp2 :: (RVector x, RVector y, Storable e, Storable f)
         => String
         -> (x e -> y f -> ST s a)
         -> x e
         -> y f
         -> ST s a
checkOp2 str f x y = do
    n1 <- getDim x
    n2 <- getDim y
    
    when (n1 /= n2) $ error $
        printf ("%s <vector with dim %d> <vector with dim %d>:"
                ++ " dimension mismatch") str n1 n2

    f x y
{-# INLINE checkOp2 #-}

checkOp3 :: (RVector x, RVector y, RVector z, Storable e, Storable f, Storable g)
         => String
         -> (x e -> y f -> z g -> ST s a)
         -> x e
         -> y f
         -> z g
         -> ST s a
checkOp3 str f x y z = do
    n1 <- getDim x
    n2 <- getDim y        
    n3 <- getDim z
    
    when (n1 /= n2 || n1 /= n3) $ error $
        printf ("%s <vector with dim %d> <vector with dim %d>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") str n1 n2 n3

    f x y z
{-# INLINE checkOp3 #-}

