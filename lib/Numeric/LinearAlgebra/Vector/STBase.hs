{-# LANGUAGE TypeSynonymInstances #-}
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
    RVector(..),
    
    sliceVector,
    dropVector,
    takeVector,
    splitVectorAt,
    
    vectorViewArray,
    
    newVector_,
    newVector,
    newCopyVector,
    newResultVector,
    newResultVector2,
    
    clearVector,
    copyToVector,
    unsafeCopyToVector,
    swapVector,
    unsafeSwapVector,
    
    indicesVector,
    getElemsVector,
    getElemsVector',
    getAssocsVector,
    getAssocsVector',
    setElemsVector,
    setAssocsVector,
    unsafeSetAssocsVector,
    
    readVector,
    unsafeReadVector,
    writeVector,
    unsafeWriteVector,
    updateVector,
    unsafeUpdateVector,
    unsafeSwapElemsVector,

    mapToVector,
    unsafeMapToVector,
    zipWithToVector,
    unsafeZipWithToVector,
    
    getSumVector,    
    getSumAbsVector,
    getNorm2Vector,
    getWhichMaxAbsVector,
    getDotVector,
    unsafeGetDotVector,
    kroneckerToVector,
    
    shiftToVector,
    addToVector,
    addToVectorWithScale,
    subToVector,
    scaleToVector,
    mulToVector,
    negateToVector,
    conjToVector,
    absToVector,
    signumToVector,
    divToVector,
    recipToVector,        
    sqrtToVector,
    expToVector,
    logToVector,
    powToVector,
    sinToVector,
    cosToVector,
    tanToVector,
    asinToVector,
    acosToVector,
    atanToVector,
    sinhToVector,
    coshToVector,
    tanhToVector,
    asinhToVector,
    acoshToVector,
    atanhToVector,
    ) where

import Control.Monad
import Control.Monad.ST
import Foreign
import System.IO.Unsafe( unsafeInterleaveIO )
import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import Data.Vector.Storable.Mutable( STVector )
import qualified Data.Vector.Storable.Mutable as STVector

import Numeric.LinearAlgebra.Internal( clearArray )
import Numeric.LinearAlgebra.Elem( Complex, VNum, VFractional, VFloating, BLAS1 )
import qualified Numeric.LinearAlgebra.Elem.BLAS as BLAS
import qualified Numeric.LinearAlgebra.Elem.VNum as VNum
import qualified Numeric.LinearAlgebra.Elem.VFractional as VFractional
import qualified Numeric.LinearAlgebra.Elem.VFloating as VFloating


-- | Read-only vectors
class RVector v where
    -- | Get the dimension of the vector.  This is equal to the number of
    -- elements in the vector.                          
    dimVector :: (Storable e) => v e -> Int

    unsafeSliceVector :: (Storable e) => Int -> Int -> v e -> v e

    -- | Execute an 'IO' action with a pointer to the first element in the
    -- vector.
    unsafeWithVector :: (Storable e) => v e -> (Ptr e -> IO a) -> IO a



instance RVector (STVector s) where
    dimVector = STVector.length
    {-# INLINE dimVector #-}

    unsafeSliceVector = STVector.unsafeSlice
    {-# INLINE unsafeSliceVector #-}

    unsafeWithVector v f =
        STVector.unsafeWith (unsafeCoerce v) f
    {-# INLINE unsafeWithVector #-}



-- | @sliceVector i n v@ creates a subvector view of @v@ starting at
-- index @i@ and having dimension @n@.
sliceVector :: (RVector v, Storable e)
            => Int
            -> Int
            -> v e
            -> v e            
sliceVector i n' v
    | i < 0 || n' < 0 || i + n' > n = error $
        printf "sliceVector %d %d <vector with dim %d>: index out of range"
               i n' n
    | otherwise =
        unsafeSliceVector i n' v
  where
    n = dimVector v
{-# INLINE sliceVector #-}

-- | @dropVector i v@ is equal to @sliceVector i (n-i) v@, where @n@ is
-- the dimension of the vector.
dropVector :: (RVector v, Storable e) => Int -> v e -> v e
dropVector i v = sliceVector i (dimVector v - i) v
{-# INLINE dropVector #-}

-- | @takeVector n v@ is equal to @sliceVector 0 n v@.
takeVector :: (RVector v, Storable e) => Int -> v e -> v e
takeVector n = sliceVector 0 n
{-# INLINE takeVector #-}

-- | View an array in memory as a vector.
vectorViewArray :: (Storable e)
                => ForeignPtr e -- ^ pointer
                -> Int          -- ^ offset
                -> Int          -- ^ length
                -> STVector s e
vectorViewArray = STVector.unsafeFromForeignPtr
{-# INLINE vectorViewArray #-}

-- | Creates a new vector of the given length.  The elements will be
-- uninitialized.
newVector_ :: (Storable e) => Int -> ST s (STVector s e)
newVector_ n
    | n < 0 =  error $
        printf "newVector_ %d: invalid dimension" n
    | otherwise = unsafeIOToST $ do
        f <- mallocForeignPtrArray n
        return $ vectorViewArray f 0 n

-- | Create a vector with every element initialized to the same value.
newVector :: (Storable e) => Int -> e -> ST s (STVector s e)
newVector n e = do
    x <- newVector_ n
    setElemsVector x $ replicate n e
    return x

-- | Creates a new vector by copying another one.    
newCopyVector :: (RVector v, Storable e) => v e -> ST s (STVector s e)
newCopyVector x = do
    y <- newVector_ (dimVector x)
    unsafeCopyToVector x y
    return y

-- | @copyToVector src dst@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyToVector :: (RVector v, Storable e) => v e -> STVector s e -> ST s ()
copyToVector = checkVectorOp2 "copyToVector" unsafeCopyToVector
{-# INLINE copyToVector #-}

unsafeCopyToVector :: (RVector v, Storable e) => v e -> STVector s e -> ST s ()
unsafeCopyToVector x y = unsafeIOToST $
    unsafeWithVector x $ \px ->
    unsafeWithVector y $ \py ->
        when (px /= py) $ copyArray py px n
  where
    n = dimVector y
{-# INLINE unsafeCopyToVector #-}

-- | Swap the values stored in two vectors.
swapVector :: (BLAS1 e) => STVector s e -> STVector s e -> ST s ()
swapVector = checkVectorOp2 "swapVector" unsafeSwapVector
{-# INLINE swapVector #-}

unsafeSwapVector :: (BLAS1 e) => STVector s e -> STVector s e -> ST s ()
unsafeSwapVector = vectorStrideCall2 BLAS.swap
{-# INLINE unsafeSwapVector #-}

-- | Split a vector into two blocks and returns views into the blocks.  In
-- @(x1, x2) = splitVectorAt k x@, we have
-- @x1 = sliceVector 0 k x@ and
-- @x2 = sliceVector k (dimVector x - k) x@.
splitVectorAt :: (RVector v, Storable e) => Int -> v e -> (v e, v e)
splitVectorAt k x
    | k < 0 || k > n = error $
        printf "splitVectorAt %d <vector with dim %d>: invalid index" k n
    | otherwise = let
        x1 = unsafeSliceVector 0 k     x
        x2 = unsafeSliceVector k (n-k) x
    in (x1,x2)
  where
    n = dimVector x
{-# INLINE splitVectorAt #-}

-- | Get the indices of the elements in the vector, @[ 0..n-1 ]@, where
-- @n@ is the dimension of the vector.
indicesVector :: (RVector v, Storable e) => v e -> [Int]
indicesVector x = [ 0..n-1 ] where n = dimVector x
{-# INLINE indicesVector #-}

-- | Lazily get the elements of the vector.
getElemsVector :: (RVector v, Storable e) => v e -> ST s [e]
getElemsVector v = unsafeIOToST $ unsafeWithVector v $ \p -> let
    n   = dimVector v
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
    touchVector v' = unsafeWithVector v' $ const (return ())
{-# SPECIALIZE INLINE getElemsVector :: STVector s Double -> ST s [Double] #-}
{-# SPECIALIZE INLINE getElemsVector :: STVector s (Complex Double) -> ST s [Complex Double] #-}

-- | Get the elements of the vector.
getElemsVector' :: (RVector v, Storable e) => v e -> ST s [e]
getElemsVector' v = unsafeIOToST $ unsafeWithVector v $ \p -> let
    n   = dimVector v
    end = p `advancePtr` (-1)
    go p' es | p' == end =
                 return es
             | otherwise = do
                 e <- peek p'
                 go (p' `advancePtr` (-1)) (e:es) in
    go (p `advancePtr` (n-1)) []
{-# SPECIALIZE INLINE getElemsVector' :: STVector s Double -> ST s [Double] #-}
{-# SPECIALIZE INLINE getElemsVector' :: STVector s (Complex Double) -> ST s [Complex Double] #-}

-- | Lazily get the association list of the vector.
getAssocsVector :: (RVector v, Storable e) => v e -> ST s [(Int,e)]
getAssocsVector x = liftM (zip (indicesVector x)) (getElemsVector x)
{-# INLINE getAssocsVector #-}

-- | Get the association list of the vector.
getAssocsVector' :: (RVector v, Storable e) => v e -> ST s [(Int,e)]
getAssocsVector' x = liftM (zip (indicesVector x)) (getElemsVector' x)
{-# INLINE getAssocsVector' #-}

-- | Set all of the values of the vector from the elements in the list.
setElemsVector  :: (Storable e) => STVector s e -> [e] -> ST s ()
setElemsVector x es = let
    n = dimVector x
    go [] i _ | i < n = error $ 
                    printf ("setElemsVector <vector with dim %d>"
                            ++ " <list with length %d>:"
                            ++ " not enough elements") n i
    go (_:_) i _ | i == n = error $ 
                    printf ("setElemsVector <vector with dim %d>"
                            ++ " <list with length at least %d>:"
                            ++ " too many elements") n (i+1)
    go []     _ _ = return ()
    go (f:fs) i p = do
        poke p f
        go fs (i+1) (p `advancePtr` 1)
    
    in unsafeIOToST $ unsafeWithVector x $ go es 0

-- | Set the given values in the vector.  If an index is repeated twice,
-- the value is implementation-defined.
setAssocsVector :: (Storable e) => STVector s e -> [(Int,e)] -> ST s ()
setAssocsVector x ies =
    let n = dimVector x
        go p ((i,e):ies') = do
            when (i < 0 || i >= n) $ error $
                printf ("setAssocsVector <vector with dim %d>"
                        ++ " [ ..., (%d,_), ... ]: invalid index") n i
            pokeElemOff p i e
            go p ies'
        go _ [] = return ()
    in unsafeIOToST $ unsafeWithVector x $ \p -> go p ies

unsafeSetAssocsVector :: (Storable e) => STVector s e -> [(Int,e)] -> ST s ()
unsafeSetAssocsVector x ies =
    let go p ((i,e):ies') = do
            pokeElemOff p i e
            go p ies'
        go _ [] = return ()
    in unsafeIOToST $ unsafeWithVector x $ \p -> go p ies

-- | Get the element stored at the given index.
readVector :: (RVector v, Storable e) => v e -> Int -> ST s e
readVector x i
    | i < 0 || i >= n = error $
        printf ("readVector <vector with dim %d> %d:"
                ++ " invalid index") n i
    | otherwise =
        unsafeReadVector x i
  where
    n = dimVector x
{-# SPECIALIZE INLINE readVector :: STVector s Double -> Int -> ST s (Double) #-}
{-# SPECIALIZE INLINE readVector :: STVector s (Complex Double) -> Int -> ST s (Complex Double) #-}

unsafeReadVector :: (RVector v, Storable e) => v e -> Int -> ST s e
unsafeReadVector x i =
    unsafeIOToST $ unsafeWithVector x $ \p -> peekElemOff p i
{-# SPECIALIZE INLINE unsafeReadVector :: STVector s Double -> Int -> ST s (Double) #-}
{-# SPECIALIZE INLINE unsafeReadVector :: STVector s (Complex Double) -> Int -> ST s (Complex Double) #-}

-- | Set the element stored at the given index.
writeVector :: (Storable e) => STVector s e -> Int -> e -> ST s ()
writeVector x i e
    | i < 0 || i >= n = error $
        printf ("writeVector <vector with dim %d> %d:"
                ++ " invalid index") n i
    | otherwise =
        unsafeWriteVector x i e
  where
    n = dimVector x
{-# SPECIALIZE INLINE writeVector :: STVector s Double -> Int -> Double -> ST s () #-}
{-# SPECIALIZE INLINE writeVector :: STVector s (Complex Double) -> Int -> Complex Double -> ST s () #-}

unsafeWriteVector :: (Storable e) => STVector s e -> Int -> e -> ST s ()
unsafeWriteVector x i e =
    unsafeIOToST $ unsafeWithVector x $ \p -> pokeElemOff p i e
{-# SPECIALIZE INLINE unsafeWriteVector :: STVector s Double -> Int -> Double -> ST s () #-}
{-# SPECIALIZE INLINE unsafeWriteVector :: STVector s (Complex Double) -> Int -> Complex Double -> ST s () #-}

-- | Modify the element stored at the given index.
updateVector :: (Storable e) => STVector s e -> Int -> (e -> e) -> ST s ()
updateVector x i f
    | i < 0 || i >= n = error $
        printf ("updateVector <vector with dim %d> %d:"
                ++ " invalid index") n i
    | otherwise =
        unsafeUpdateVector x i f
  where
    n = dimVector x
{-# SPECIALIZE INLINE updateVector :: STVector s Double -> Int -> (Double -> Double) -> ST s () #-}
{-# SPECIALIZE INLINE updateVector :: STVector s (Complex Double) -> Int -> (Complex Double -> Complex Double) -> ST s () #-}

unsafeUpdateVector :: (Storable e) => STVector s e -> Int -> (e -> e) -> ST s ()
unsafeUpdateVector x i f =
    unsafeIOToST $ unsafeWithVector x $ \p -> do
        e <- peekElemOff p i
        pokeElemOff p i $ f e
{-# SPECIALIZE INLINE unsafeUpdateVector :: STVector s Double -> Int -> (Double -> Double) -> ST s () #-}
{-# SPECIALIZE INLINE unsafeUpdateVector :: STVector s (Complex Double) -> Int -> (Complex Double -> Complex Double) -> ST s () #-}

unsafeSwapElemsVector :: (Storable e) => STVector s e -> Int -> Int -> ST s ()
unsafeSwapElemsVector x i1 i2 = unsafeIOToST $ unsafeWithVector x $ \p ->
    let p1 = p `advancePtr` i1
        p2 = p `advancePtr` i2
    in do
        e1  <- peek p1
        e2  <- peek p2
        poke p2 e1
        poke p1 e2
{-# SPECIALIZE INLINE unsafeSwapElemsVector :: STVector s Double -> Int -> Int -> ST s () #-}
{-# SPECIALIZE INLINE unsafeSwapElemsVector :: STVector s (Complex Double) -> Int -> Int -> ST s () #-}

-- | @mapToVector f x z@ replaces @z@ elementwise with @f(x)@.
mapToVector :: (RVector v, Storable e, Storable f)
            => (e -> f)
            -> v e
            -> STVector s f
            -> ST s ()
mapToVector f = checkVectorOp2 "mapToVector _" $ unsafeMapToVector f
{-# INLINE mapToVector #-}
                            
unsafeMapToVector :: (RVector v, Storable e, Storable f)
                  => (e -> f)
                  -> v e
                  -> STVector s f
                  -> ST s ()
unsafeMapToVector f v mv =
    let go psrc pdst end
            | pdst == end =
                return ()
            | otherwise = do
                e <- peek psrc
                poke pdst (f e)
                go (psrc `advancePtr` 1) (pdst `advancePtr` 1) end
    in unsafeIOToST $
           unsafeWithVector v $ \psrc -> 
           unsafeWithVector mv $ \pdst ->
               go psrc pdst (pdst `advancePtr` dimVector mv)
  where

{-# INLINE unsafeMapToVector #-}

-- | @zipWithToVector f x y z@ replaces @z@ elementwise with @f(x, y)@.
zipWithToVector :: (RVector v1, RVector v2, Storable e1, Storable e2, Storable f)
                => (e1 -> e2 -> f)
                -> v1 e1
                -> v2 e2
                -> STVector s f
                -> ST s ()
zipWithToVector f = checkVectorOp3 "zipWithToVector _" $
    unsafeZipWithToVector f
{-# INLINE zipWithToVector #-}

unsafeZipWithToVector :: (RVector v1, RVector v2, Storable e1, Storable e2, Storable f)
                      => (e1 -> e2 -> f)
                      -> v1 e1
                      -> v2 e2
                      -> STVector s f
                      -> ST s ()
unsafeZipWithToVector f v1 v2 mv =
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
           unsafeWithVector v1 $ \psrc1 ->
           unsafeWithVector v2 $ \psrc2 -> 
           unsafeWithVector mv $ \pdst ->
               go psrc1 psrc2 pdst (pdst `advancePtr` dimVector mv)
{-# INLINE unsafeZipWithToVector #-}


-- | Set every element in the vector to a default value.  For
-- standard numeric types (including 'Double', 'Complex Double', and 'Int'),
-- the default value is '0'.
clearVector :: (Storable e) => STVector s e -> ST s ()
clearVector x = unsafeIOToST $ unsafeWithVector x $ \p -> clearArray p n
  where
    n = dimVector x

-- | @negateToVector x z@ replaces @z@ with @negate(x)@.
negateToVector :: (RVector v, VNum e) => v e -> STVector s e -> ST s ()
negateToVector = checkVectorOp2 "negateToVector" $
    vectorCall2 VNum.vNeg
{-# INLINE negateToVector #-}

-- | @absToVector x z@ replaces @z@ with @abs(x)@.
absToVector :: (RVector v, VNum e) => v e -> STVector s e -> ST s ()
absToVector = checkVectorOp2 "absToVector" $
    vectorCall2 VNum.vAbs
{-# INLINE absToVector #-}

-- | @signumToVector x z@ replaces @z@ with @signum(x)@.
signumToVector :: (RVector v, VNum e) => v e -> STVector s e -> ST s ()
signumToVector = checkVectorOp2 "signumToVector" $
    vectorCall2 VNum.vSgn
{-# INLINE signumToVector #-}



-- | @conjToVector x z@ replaces @z@ with @conj(x)@.
conjToVector :: (RVector v, VNum e) => v e -> STVector s e -> ST s ()
conjToVector = checkVectorOp2 "conjToVector" $
    vectorCall2 VNum.vConj
{-# INLINE conjToVector #-}

-- | @shiftToVector alpha x z@ replaces @z@ with @alpha + x@.
shiftToVector :: (RVector v, VNum e) => e -> v e -> STVector s e -> ST s ()
shiftToVector alpha = checkVectorOp2 "shiftToVector _" $
    vectorCall2 (flip VNum.vShift alpha)
{-# INLINE shiftToVector #-}

-- | @scaleToVector alpha x z@ replaces @z@ with @alpha * x@.
scaleToVector :: (RVector v, VNum e) => e -> v e -> STVector s e -> ST s ()
scaleToVector alpha = checkVectorOp2 "scaleToVector _" $
    vectorCall2 (flip VNum.vScale alpha)
{-# INLINE scaleToVector #-}

-- | @addToVectorWithScale alpha x beta y z@ replaces @z@ with
-- @alpha * x + beta * y@.
addToVectorWithScale :: (RVector v1, RVector v2, VNum e)
                     => e -> v1 e -> e -> v2 e -> STVector s e -> ST s ()
addToVectorWithScale alpha x beta y z =
    checkVectorOp3 "addToVectorWithScale"
        (vectorCall3 (\n v1 v2 v3 -> VNum.vAxpby n alpha v1 beta v2 v3))
        x y z
{-# INLINE addToVectorWithScale #-}

-- | @addToVector x y z@ replaces @z@ with @x+y@.
addToVector :: (RVector v1, RVector v2, VNum e)
            => v1 e -> v2 e -> STVector s e -> ST s ()
addToVector = checkVectorOp3 "addToVector" $ vectorCall3 VNum.vAdd
{-# INLINE addToVector #-}

-- | @subToVector x y z@ replaces @z@ with @x-y@.
subToVector :: (RVector v1, RVector v2, VNum e)
            => v1 e -> v2 e -> STVector s e -> ST s ()
subToVector = checkVectorOp3 "subToVector" $ vectorCall3 VNum.vSub
{-# INLINE subToVector #-}

-- | @mulToVector x y z@ replaces @z@ with @x*y@.
mulToVector :: (RVector v1, RVector v2, VNum e)
            => v1 e -> v2 e -> STVector s e -> ST s ()
mulToVector = checkVectorOp3 "mulToVector" $ vectorCall3 VNum.vMul
{-# INLINE mulToVector #-}
 
-- | @divToVector x y z@ replaces @z@ with @x/y@.
divToVector :: (RVector v1, RVector v2, VFractional e)
            => v1 e -> v2 e -> STVector s e -> ST s ()
divToVector = checkVectorOp3 "divToVector" $ vectorCall3 VFractional.vDiv
{-# INLINE divToVector #-}

-- | @recipToVector x z@ replaces @z@ with @1/x@.
recipToVector :: (RVector v, VFractional e)
              => v e -> STVector s e -> ST s ()
recipToVector = checkVectorOp2 "recipToVector" $ vectorCall2 VFractional.vInv
{-# INLINE recipToVector #-}

-- | @sqrtToVector x z@ replaces @z@ with @sqrt(x)@.
sqrtToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
sqrtToVector = checkVectorOp2 "sqrtToVector" $
    vectorCall2 VFloating.vSqrt
{-# INLINE sqrtToVector #-}

-- | @expToVector x z@ replaces @z@ with @exp(x)@.
expToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
expToVector = checkVectorOp2 "expToVector" $
    vectorCall2 VFloating.vExp
{-# INLINE expToVector #-}

-- | @logToVector x z@ replaces @z@ with @log(x)@.
logToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
logToVector = checkVectorOp2 "logToVector" $
    vectorCall2 VFloating.vLog
{-# INLINE logToVector #-}

-- | @powToVector x y z@ replaces @z@ with @x ** y@.
powToVector :: (RVector v1, RVector v2, VFloating e)
            => v1 e -> v2 e -> STVector s e -> ST s ()
powToVector = checkVectorOp3 "powToVector" $
    vectorCall3 VFloating.vPow
{-# INLINE powToVector #-}

-- | @sinToVector x z@ replaces @z@ with @sin(x)@.
sinToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
sinToVector = checkVectorOp2 "sinToVector" $
    vectorCall2 VFloating.vSin
{-# INLINE sinToVector #-}

-- | @cosToVector x z@ replaces @z@ with @cos(x)@.
cosToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
cosToVector = checkVectorOp2 "cosToVector" $
    vectorCall2 VFloating.vCos
{-# INLINE cosToVector #-}

-- | @tanToVector x z@ replaces @z@ with @tan(x)@.
tanToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
tanToVector = checkVectorOp2 "tanToVector" $
    vectorCall2 VFloating.vTan
{-# INLINE tanToVector #-}

-- | @asinToVector x z@ replaces @z@ with @asin(x)@.
asinToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
asinToVector = checkVectorOp2 "asinToVector" $
    vectorCall2 VFloating.vASin
{-# INLINE asinToVector #-}

-- | @acosToVector x z@ replaces @z@ with @acos(x)@.
acosToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
acosToVector = checkVectorOp2 "acosToVector" $
    vectorCall2 VFloating.vACos
{-# INLINE acosToVector #-}

-- | @atanToVector x z@ replaces @z@ with @atan(x)@.
atanToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
atanToVector = checkVectorOp2 "atanToVector" $
    vectorCall2 VFloating.vATan
{-# INLINE atanToVector #-}

-- | @sinhToVector x z@ replaces @z@ with @sinh(x)@.
sinhToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
sinhToVector = checkVectorOp2 "sinhToVector" $
    vectorCall2 VFloating.vSinh
{-# INLINE sinhToVector #-}

-- | @coshToVector x z@ replaces @z@ with @cosh(x)@.
coshToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
coshToVector = checkVectorOp2 "coshToVector" $
    vectorCall2 VFloating.vCosh
{-# INLINE coshToVector #-}

-- | @tanhToVector x z@ replaces @z@ with @tanh(x)@.
tanhToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
tanhToVector = checkVectorOp2 "tanhToVector" $
    vectorCall2 VFloating.vTanh
{-# INLINE tanhToVector #-}

-- | @asinhToVector x z@ replaces @z@ with @asinh(x)@.
asinhToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
asinhToVector = checkVectorOp2 "asinhToVector" $
    vectorCall2 VFloating.vASinh
{-# INLINE asinhToVector #-}

-- | @acoshToVector x z@ replaces @z@ with @acosh(x)@.
acoshToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
acoshToVector = checkVectorOp2 "acoshToVector" $
    vectorCall2 VFloating.vACosh
{-# INLINE acoshToVector #-}

-- | @atanhToVector x z@ replaces @z@ with @atanh(x)@.
atanhToVector :: (RVector v, VFloating e) => v e -> STVector s e -> ST s ()
atanhToVector = checkVectorOp2 "atanhToVector" $
    vectorCall2 VFloating.vATanh
{-# INLINE atanhToVector #-}

-- | Gets the sum of the vector entries.
getSumVector :: (RVector v, VNum e) => v e -> ST s e
getSumVector = vectorCall VNum.vSum
{-# INLINE getSumVector #-}

-- | Gets the sum of the absolute values of the vector entries.
getSumAbsVector :: (RVector v, BLAS1 e) => v e -> ST s Double
getSumAbsVector = vectorStrideCall BLAS.asum
{-# INLINE getSumAbsVector #-}

-- | Gets the 2-norm of a vector.
getNorm2Vector :: (RVector v, BLAS1 e) => v e -> ST s Double
getNorm2Vector = vectorStrideCall BLAS.nrm2
{-# INLINE getNorm2Vector #-}

-- | Gets the index and norm of the element with maximum magnitude.  This is 
-- undefined if any of the elements are @NaN@.  It will throw an exception if 
-- the dimension of the vector is 0.
getWhichMaxAbsVector :: (RVector v, BLAS1 e) => v e -> ST s (Int, e)
getWhichMaxAbsVector x =
    case (dimVector x) of
        0 -> error $ "getWhichMaxAbs <vector with dim 0>: empty vector"
        _ -> do
            i <- vectorStrideCall BLAS.iamax x
            e <- unsafeReadVector x i
            return (i,e)
{-# INLINE getWhichMaxAbsVector #-}

-- | Computes the dot product of two vectors.
getDotVector :: (RVector v, RVector v', BLAS1 e)
             => v e -> v' e -> ST s e
getDotVector = checkVectorOp2 "getDotVector" unsafeGetDotVector
{-# INLINE getDotVector #-}

unsafeGetDotVector :: (RVector x, RVector y, BLAS1 e)
                   => x e -> y e -> ST s e
unsafeGetDotVector x y = (vectorStrideCall2 BLAS.dotc) y x
{-# INLINE unsafeGetDotVector #-}

-- | @kroneckerToVector x y z@ sets @z := x \otimes y@.
kroneckerToVector :: (RVector v1, RVector v2, VNum e)
                  => v1 e -> v2 e -> STVector s e -> ST s ()
kroneckerToVector x y z
    | dimVector z /= m * n = error $
        printf ("kroneckerToVector"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (dimVector z)
    | otherwise = do
        ies <- getAssocsVector x
        forM_ ies $ \(i,e) ->
            scaleToVector e y (unsafeSliceVector (i*n) n z)
  where
    m = dimVector x
    n = dimVector y

vectorCall :: (RVector x, Storable e)
           => (Int -> Ptr e -> IO a) 
           ->  x e -> ST s a
vectorCall f x = 
    let n    = dimVector x
    in unsafeIOToST $
           unsafeWithVector x $ \pX ->
               f n pX
{-# INLINE vectorCall #-}

vectorCall2 :: (RVector x, RVector y, Storable e, Storable f)
            => (Int -> Ptr e -> Ptr f -> IO a) 
            -> x e -> y f -> ST s a
vectorCall2 f x y =
    let n    = dimVector x
    in unsafeIOToST $
           unsafeWithVector x $ \pX ->
           unsafeWithVector y $ \pY ->
               f n pX pY
{-# INLINE vectorCall2 #-}    

vectorCall3 :: (RVector x, RVector y, RVector z, Storable e, Storable f, Storable g)
            => (Int -> Ptr e -> Ptr f -> Ptr g -> IO a) 
            -> x e -> y f -> z g -> ST s a
vectorCall3 f x y z =
    let n    = dimVector x
    in unsafeIOToST $
           unsafeWithVector x $ \pX ->
           unsafeWithVector y $ \pY ->
           unsafeWithVector z $ \pZ ->           
               f n pX pY pZ
{-# INLINE vectorCall3 #-}   

vectorStrideCall :: (RVector x, Storable e)
                 => (Int -> Ptr e -> Int -> IO a) 
                 ->  x e -> ST s a
vectorStrideCall f x = 
    let n    = dimVector x
        incX = 1
    in unsafeIOToST $
           unsafeWithVector x $ \pX ->
               f n pX incX
{-# INLINE vectorStrideCall #-}

vectorStrideCall2 :: (RVector x, RVector y, Storable e, Storable f)
                  => (Int -> Ptr e -> Int -> Ptr f -> Int -> IO a) 
                  -> x e -> y f -> ST s a
vectorStrideCall2 f x y =
    let n    = dimVector x
        incX = 1
        incY = 1
    in unsafeIOToST $
           unsafeWithVector x $ \pX ->
           unsafeWithVector y $ \pY ->
               f n pX incX pY incY
{-# INLINE vectorStrideCall2 #-}    

checkVectorOp2 :: (RVector x, RVector y, Storable e, Storable f)
               => String
               -> (x e -> y f -> a)
               -> x e
               -> y f
               -> a
checkVectorOp2 str f x y
    | n1 /= n2 = error $
        printf ("%s <vector with dim %d> <vector with dim %d>:"
                ++ " dimension mismatch") str n1 n2
    | otherwise =
        f x y
  where
    n1 = dimVector x
    n2 = dimVector y        
{-# INLINE checkVectorOp2 #-}

checkVectorOp3 :: (RVector x, RVector y, RVector z, Storable e, Storable f, Storable g)
               => String
               -> (x e -> y f -> z g -> a)
               -> x e
               -> y f
               -> z g
               -> a
checkVectorOp3 str f x y z
    | n1 /= n2 || n1 /= n3 = error $
        printf ("%s <vector with dim %d> <vector with dim %d>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") str n1 n2 n3
    | otherwise =
        f x y z
  where
    n1 = dimVector x
    n2 = dimVector y        
    n3 = dimVector z
{-# INLINE checkVectorOp3 #-}

newResultVector :: (RVector v, Storable e, Storable f)
                => (v e -> STVector s f -> ST s a)
                -> v e
                -> ST s (STVector s f)
newResultVector f v = do
    z <- newVector_ (dimVector v)
    _ <- f v z
    return z
{-# INLINE newResultVector #-}


newResultVector2 :: (RVector v1, RVector v2, Storable e, Storable f, Storable g)
                 => (v1 e -> v2 f -> STVector s g -> ST s a)
                 -> v1 e
                 -> v2 f
                 -> ST s (STVector s g)
newResultVector2 f v1 v2 = do
    z <- newVector_ (dimVector v1)
    _ <- f v1 v2 z
    return z
{-# INLINE newResultVector2 #-}
