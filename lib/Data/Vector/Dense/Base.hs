{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances,
        Rank2Types #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Base
    where

import Control.Monad
import Control.Monad.Interleave
import Control.Monad.ST
import Data.AEq
import Foreign
import Text.Printf

import BLAS.Internal( checkBinaryOp, inlinePerformIO,
    checkedSubvector, checkedSubvectorWithStride, checkVecVecOp )
import BLAS.Types( ConjEnum(..) )

import Data.Elem.BLAS
import qualified Data.Elem.BLAS.Base   as BLAS
import qualified Data.Elem.BLAS.Level1 as BLAS

import Data.Tensor.Class
import Data.Tensor.Class.ITensor
import Data.Tensor.Class.MTensor

import Data.Vector.Dense.IOBase

infixl 7 <.>

-- | Immutable dense vectors. The type arguments are as follows:
--
--     * @e@: the element type of the vector.  Only certain element types
--       are supported.
--
newtype Vector e = Vector (IOVector e)

freezeIOVector :: IOVector e -> IO (Vector e)
freezeIOVector = liftM Vector . newCopyIOVector
{-# INLINE freezeIOVector #-}

thawIOVector :: Vector e -> IO (IOVector e)
thawIOVector (Vector x) = newCopyIOVector x
{-# INLINE thawIOVector #-}

unsafeFreezeIOVector :: IOVector e -> IO (Vector e)
unsafeFreezeIOVector = return . Vector
{-# INLINE unsafeFreezeIOVector #-}

unsafeThawIOVector :: Vector e -> IO (IOVector e)
unsafeThawIOVector (Vector x) = return x
{-# INLINE unsafeThawIOVector #-}

-- | Common functionality for all vector types.
class (Shaped x Int) => BaseVector x where
    -- | Get the dimension (length) of the vector.
    dim :: x e -> Int
    
    -- | Get the memory stride (in elements) between consecutive elements.
    stride :: x e -> Int

    -- | Indicate whether or not internally the vector stores the complex
    -- conjugates of its elements.
    isConj :: x e -> Bool
    isConj x = conjEnum x == Conj
    {-# INLINE isConj #-}

    -- | Get the storage type.
    conjEnum :: x e -> ConjEnum
    
    -- | Get a view into the complex conjugate of a vector.
    conj :: x e -> x e

    unsafeSubvectorViewWithStride :: Int -> x e -> Int -> Int -> x e

    -- | Unsafe cast from a vector to an 'IOVector'.
    unsafeVectorToIOVector :: x e -> IOVector e
    unsafeIOVectorToVector :: IOVector e -> x e


-- | Vectors that can be read in a monad.
class (BaseVector x, MonadInterleave m, ReadTensor x Int m) => ReadVector x m where
    -- | Convert a mutable vector to an immutable one by taking a complete
    -- copy of it.
    freezeVector :: x e -> m (Vector e)
    unsafeFreezeVector :: x e -> m (Vector e)

    -- | Cast the vector to an 'IOVector', perform an @IO@ action, and
    -- convert the @IO@ action to an action in the monad @m@.
    unsafePerformIOWithVector :: x e -> (IOVector e -> IO a) -> m a


-- | Vectors that can be created or modified in a monad.
class (ReadVector x m, WriteTensor x Int m) => WriteVector x m where
    -- | Creates a new vector of the given length.  The elements will be
    -- uninitialized.
    newVector_ :: (Elem e) => Int -> m (x e)

    -- | Convert an immutable vector to a mutable one by taking a complete
    -- copy of it.
    thawVector :: Vector e -> m (x e)
    unsafeThawVector :: Vector e -> m (x e)

    -- | Unsafely convert an 'IO' action that creates an 'IOVector' into
    -- an action in @m@ that creates a vector.
    unsafeConvertIOVector :: IO (IOVector e) -> m (x e)




-- | Creates a new vector with the given association list.  Unspecified
-- indices will get initialized to zero.
newVector :: (WriteVector x m, Elem e) => Int -> [(Int,e)] -> m (x e)
newVector n ies = 
    unsafeConvertIOVector $ newIOVector n ies
{-# INLINE newVector #-}

unsafeNewVector :: (WriteVector x m, Elem e) => Int -> [(Int,e)] -> m (x e)
unsafeNewVector n ies = 
    unsafeConvertIOVector $ unsafeNewIOVector n ies
{-# INLINE unsafeNewVector #-}

-- | Creates a new vector of the given dimension with the given elements.
-- If the list has length less than the passed-in dimenson, the tail of
-- the vector will be uninitialized.
newListVector :: (WriteVector x m, Elem e) => Int -> [e] -> m (x e)
newListVector n es = 
    unsafeConvertIOVector $ newListIOVector n es
{-# INLINE newListVector #-}

-- | Create a zero vector of the specified length.
newZeroVector :: (WriteVector x m, Elem e) => Int -> m (x e)
newZeroVector n =
    unsafeConvertIOVector $ newZeroIOVector n
{-# INLINE newZeroVector #-}

-- | Set every element in the vector to zero.
setZeroVector :: (WriteVector x m) => x e -> m ()
setZeroVector x =
    unsafePerformIOWithVector x $ setZeroIOVector
{-# INLINE setZeroVector #-}

-- | Create a vector with every element initialized to the same value.
newConstantVector :: (WriteVector x m, Elem e) => Int -> e -> m (x e)
newConstantVector n e = 
    unsafeConvertIOVector $ newConstantIOVector n e
{-# INLINE newConstantVector #-}

-- | Set every element in the vector to a constant.
setConstantVector :: (WriteVector x m) => e -> x e -> m ()
setConstantVector e x =
    unsafePerformIOWithVector x $ setConstantIOVector e
{-# INLINE setConstantVector #-}

-- | @newBasisVector n i@ creates a vector of length @n@ that is all zero 
-- except for at position @i@, where it equal to one.
newBasisVector :: (WriteVector x m, Elem e) => Int -> Int -> m (x e)
newBasisVector n i = unsafeConvertIOVector $
    newBasisIOVector n i
{-# INLINE newBasisVector #-}

-- | @setBasis x i@ sets the @i@th coordinate of @x@ to @1@, and all other
-- coordinates to @0@.  If the vector has been scaled, it is possible that
-- @readVector x i@ will not return exactly @1@.  See 'setElem'.
setBasisVector :: (WriteVector x m) => Int -> x e -> m ()
setBasisVector i x =
    unsafePerformIOWithVector x $ setBasisIOVector i
{-# INLINE setBasisVector #-}

-- | Creates a new vector by copying another one.
newCopyVector :: (ReadVector x m, WriteVector y m) =>
    x e -> m (y e)
newCopyVector x = unsafeConvertIOVector $
    newCopyIOVector (unsafeVectorToIOVector x)
{-# INLINE newCopyVector #-}

-- | Creates a new vector by copying another one.  The returned vector
-- is gauranteed not to be a view into another vector.  That is, the
-- returned vector will have @isConj@ to be @False@.
newCopyVector' :: (ReadVector x m, WriteVector y m) => x e -> m (y e)
newCopyVector' x = unsafeConvertIOVector $
    newCopyIOVector' (unsafeVectorToIOVector x)
{-# INLINE newCopyVector' #-}

-- | @copyVector dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyVector :: (WriteVector y m, ReadVector x m) =>
    y e -> x e -> m ()
copyVector y x = checkBinaryOp (shape x) (shape y) $ unsafeCopyVector y x
{-# INLINE copyVector #-}

unsafeCopyVector :: (WriteVector y m, ReadVector x m) =>
    y e -> x e -> m ()
unsafeCopyVector y x =
    unsafePerformIOWithVector y $
        (`unsafeCopyIOVector` (unsafeVectorToIOVector x))
{-# INLINE unsafeCopyVector #-}

-- | Swap the values stored in two vectors.
swapVector :: (WriteVector x m, WriteVector y m) => 
    x e -> y e -> m ()
swapVector x y = checkBinaryOp (shape x) (shape y) $ unsafeSwapVector x y
{-# INLINE swapVector #-}

unsafeSwapVector :: (WriteVector x m, WriteVector y m) =>
    x e -> y e -> m ()
unsafeSwapVector x y =
    unsafePerformIOWithVector x $ 
        (`unsafeSwapIOVector` (unsafeVectorToIOVector y))
{-# INLINE unsafeSwapVector #-}

-- | @subvectorView x o n@ creates a subvector view of @x@ starting at index @o@ 
-- and having length @n@.
subvectorView :: (BaseVector x) => 
    x e -> Int -> Int -> x e
subvectorView x = checkedSubvector (dim x) (unsafeSubvectorView x)
{-# INLINE subvectorView #-}

unsafeSubvectorView :: (BaseVector x) => 
    x e -> Int -> Int -> x e
unsafeSubvectorView = unsafeSubvectorViewWithStride 1
{-# INLINE unsafeSubvectorView #-}

-- | @subvectorViewWithStride s x o n@ creates a subvector view of @x@ starting 
-- at index @o@, having length @n@ and stride @s@.
subvectorViewWithStride :: (BaseVector x) => 
    Int -> x e -> Int -> Int -> x e
subvectorViewWithStride s x = 
    checkedSubvectorWithStride s (dim x) (unsafeSubvectorViewWithStride s x)
{-# INLINE subvectorViewWithStride #-}

-- | Split a vector into two blocks and returns views into the blocks.  In
-- @(x1, x2) = splitElemsAt k x@, we have that @x1 = subvectorView x 0 k'@
-- and @x2 = subvectorView x k (dim x - k')@, where @k'@ is defined
-- as @k' = max 0 (min (dim x) k)@.
splitElemsAt :: (BaseVector x) => Int -> x e -> (x e, x e)
splitElemsAt k x =
    let n  = dim x
        k' = max 0 $ min n k
        x1 = unsafeSubvectorView x 0  k'
        x2 = unsafeSubvectorView x k' (n-k')
    in (x1,x2)

-- | Get a new vector with elements with the conjugates of the elements
-- of the given vector
getConjVector :: (ReadVector x m, WriteVector y m) =>
    x e -> m (y e)
getConjVector = getUnaryVectorOp doConjVector
{-# INLINE getConjVector #-}

-- | Conjugate every element of the vector.
doConjVector :: (WriteVector y m) => y e -> m ()
doConjVector x =
    unsafePerformIOWithVector x $ doConjIOVector
{-# INLINE doConjVector #-}

-- | Get a new vector by scaling the elements of another vector
-- by a given value.
getScaledVector :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> x e -> m (y e)
getScaledVector e = getUnaryVectorOp (scaleByVector e)
{-# INLINE getScaledVector #-}

-- | Scale every element by the given value.
scaleByVector :: (WriteVector y m, BLAS1 e) => e -> y e -> m ()
scaleByVector k x =
    unsafePerformIOWithVector x $ scaleByIOVector k
{-# INLINE scaleByVector #-}

-- | Get a new vector by shifting the elements of another vector
-- by a given value.
getShiftedVector :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> x e -> m (y e)
getShiftedVector e = getUnaryVectorOp (shiftByVector e)
{-# INLINE getShiftedVector #-}

-- | Add the given value to every element.
shiftByVector :: (WriteVector y m, BLAS1 e) => e -> y e -> m ()
shiftByVector k x =
    unsafePerformIOWithVector x $ shiftByIOVector k
{-# INLINE shiftByVector #-}

-- | @getAddVector x y@ creates a new vector equal to the sum @x+y@.  The 
-- operands must have the same dimension.
getAddVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x e -> y e -> m (z e)
getAddVector = checkVectorOp2 unsafeGetAddVector
{-# INLINE getAddVector #-}

unsafeGetAddVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x e -> y e -> m (z e)
unsafeGetAddVector = unsafeGetBinaryVectorOp unsafeAddVector
{-# INLINE unsafeGetAddVector #-}

-- | @addVector y x@ replaces @y@ with @y+x@.
addVector :: (WriteVector y m, ReadVector x m, BLAS1 e) => 
    y e -> x e -> m ()
addVector y x = checkBinaryOp (dim y) (dim x) $ unsafeAddVector y x
{-# INLINE addVector #-}

unsafeAddVector :: (WriteVector y m, ReadVector x m, BLAS1 e) =>
    y e -> x e -> m ()
unsafeAddVector y x = unsafeAxpyVector 1 x y
{-# INLINE unsafeAddVector #-}

-- | @getSubVector x y@ creates a new tensor equal to the difference @x-y@.  
-- The operands must have the same dimension.
getSubVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x e -> y e -> m (z e)
getSubVector = checkVectorOp2 unsafeGetSubVector
{-# INLINE getSubVector #-}

unsafeGetSubVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x e -> y e -> m (z e)
unsafeGetSubVector = unsafeGetBinaryVectorOp unsafeSubVector
{-# INLINE unsafeGetSubVector #-}

-- | @subVector y x@ replaces @y@ with @y-x@.
subVector :: (WriteVector y m, ReadVector x m, BLAS1 e) => 
    y e -> x e -> m ()
subVector y x = checkBinaryOp (dim y) (dim x) $ unsafeSubVector y x
{-# INLINE subVector #-}

unsafeSubVector :: (WriteVector y m, ReadVector x m, BLAS1 e) =>
    y e -> x e -> m ()
unsafeSubVector y x = unsafeAxpyVector (-1) x y
{-# INLINE unsafeSubVector #-}

-- | @axpyVector alpha x y@ replaces @y@ with @alpha * x + y@.
axpyVector :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> x e -> y e -> m ()
axpyVector alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyVector alpha x y
{-# INLINE axpyVector #-}

unsafeAxpyVector :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> x e -> y e -> m ()
unsafeAxpyVector alpha x y
    | isConj y =
        unsafeAxpyVector (conjugate alpha) (conj x) (conj y)
    | isConj x =
        vectorCall2 (flip BLAS.acxpy alpha) x y
    | otherwise =
        vectorCall2 (flip BLAS.axpy alpha) x y
{-# INLINE unsafeAxpyVector #-}

-- | @getMulVector x y@ creates a new vector equal to the elementwise product 
-- @x*y@.  The operands must have the same dimensino
getMulVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x e -> y e -> m (z e)
getMulVector = checkVectorOp2 unsafeGetMulVector
{-# INLINE getMulVector #-}

unsafeGetMulVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x e -> y e -> m (z e)
unsafeGetMulVector = unsafeGetBinaryVectorOp unsafeMulVector
{-# INLINE unsafeGetMulVector #-}

-- | @mulVector y x@ replaces @y@ with @y*x@.
mulVector :: (WriteVector y m, ReadVector x m, BLAS1 e) => 
    y e -> x e -> m ()
mulVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeMulVector y x
{-# INLINE mulVector #-}
 
unsafeMulVector :: (WriteVector y m, ReadVector x m, BLAS1 e) =>
    y e -> x e -> m ()
unsafeMulVector y x
    | isConj y =
        unsafeMulVector (conj y) (conj x)
    | isConj x =
        vectorCall2 BLAS.vcmul x y
    | otherwise =
        vectorCall2 BLAS.vmul x y
{-# INLINE unsafeMulVector #-}

-- | @getDivVector x y@ creates a new vector equal to the elementwise 
-- ratio @x/y@.  The operands must have the same shape.
getDivVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x e -> y e -> m (z e)
getDivVector = checkVectorOp2 unsafeGetDivVector
{-# INLINE getDivVector #-}

unsafeGetDivVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x e -> y e -> m (z e)
unsafeGetDivVector = unsafeGetBinaryVectorOp unsafeDivVector
{-# INLINE unsafeGetDivVector #-}

-- | @divVector y x@ replaces @y@ with @y/x@.
divVector :: (WriteVector y m, ReadVector x m, BLAS1 e) => 
    y e -> x e -> m ()
divVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeDivVector y x
{-# INLINE divVector #-}

unsafeDivVector :: (WriteVector y m, ReadVector x m, BLAS1 e) =>
    y e -> x e -> m ()
unsafeDivVector y x
    | isConj y =
        unsafeDivVector (conj y) (conj x)
    | isConj x =
        vectorCall2 BLAS.vcdiv x y
    | otherwise =
        vectorCall2 BLAS.vdiv x y
{-# INLINE unsafeDivVector #-}

-- | Gets the sum of the absolute values of the vector entries.
getSumAbs :: (ReadVector x m, BLAS1 e) => x e -> m Double
getSumAbs = vectorCall BLAS.asum
{-# INLINE getSumAbs #-}
    
-- | Gets the 2-norm of a vector.
getNorm2 :: (ReadVector x m, BLAS1 e) => x e -> m Double
getNorm2 = vectorCall BLAS.nrm2
{-# INLINE getNorm2 #-}

-- | Gets the index and norm of the element with maximum magnitude.  This is 
-- undefined if any of the elements are @NaN@.  It will throw an exception if 
-- the dimension of the vector is 0.
getWhichMaxAbs :: (ReadVector x m, BLAS1 e) => x e -> m (Int, e)
getWhichMaxAbs x =
    case (dim x) of
        0 -> fail $ "getWhichMaxAbs of an empty vector"
        _ -> do
            i <- vectorCall BLAS.iamax x
            e <- unsafeReadElem x i
            return (i,e)
{-# INLINE getWhichMaxAbs #-}

-- | Computes the dot product of two vectors.
getDot :: (ReadVector x m, ReadVector y m, BLAS1 e) => 
    x e -> y e -> m e
getDot x y = checkVecVecOp "getDot" (dim x) (dim y) $ unsafeGetDot x y
{-# INLINE getDot #-}

unsafeGetDot :: (ReadVector x m, ReadVector y m, BLAS1 e) => 
    x e -> y e -> m e
unsafeGetDot x y =
    vectorCall2 (BLAS.dot (conjEnum x) (conjEnum y)) x y
{-# INLINE unsafeGetDot #-}

instance BaseVector IOVector where
    dim = dimIOVector
    {-# INLINE dim #-}
    stride = strideIOVector
    {-# INLINE stride #-}
    conjEnum = conjEnumIOVector
    {-# INLINE conjEnum #-}
    conj = conjIOVector
    {-# INLINE conj #-}
    unsafeSubvectorViewWithStride = unsafeSubvectorViewWithStrideIOVector
    {-# INLINE unsafeSubvectorViewWithStride #-}
    unsafeVectorToIOVector = id
    {-# INLINE unsafeVectorToIOVector #-}
    unsafeIOVectorToVector = id
    {-# INLINE unsafeIOVectorToVector #-}
    
instance ReadVector IOVector IO where
    freezeVector = freezeIOVector
    {-# INLINE freezeVector #-}
    unsafeFreezeVector = unsafeFreezeIOVector
    {-# INLINE unsafeFreezeVector #-}
    unsafePerformIOWithVector x f = f x
    {-# INLINE unsafePerformIOWithVector #-}

instance WriteVector IOVector IO where
    newVector_ = newIOVector_
    {-# INLINE newVector_ #-}
    thawVector = thawIOVector
    {-# INLINE thawVector #-}
    unsafeThawVector = unsafeThawIOVector
    {-# INLINE unsafeThawVector #-}
    unsafeConvertIOVector = id
    {-# INLINE unsafeConvertIOVector #-}

-- | Create a vector with the given dimension and elements.  The elements
-- given in the association list must all have unique indices, otherwise
-- the result is undefined.
vector :: (Elem e) => Int -> [(Int, e)] -> Vector e
vector n ies = unsafePerformIO $
    unsafeFreezeIOVector =<< newIOVector n ies
{-# INLINE vector #-}

-- Same as 'vector', but does not range-check the indices.
unsafeVector :: (Elem e) => Int -> [(Int, e)] -> Vector e
unsafeVector n ies = unsafePerformIO $
    unsafeFreezeIOVector =<< unsafeNewIOVector n ies
{-# INLINE unsafeVector #-}

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.  @listVector n es@ is equivalent to 
-- @vector n (zip [0..(n-1)] es)@, except that the result is undefined 
-- if @length es@ is less than @n@.
listVector :: (Elem e) => Int -> [e] -> Vector e
listVector n es = unsafePerformIO $
    unsafeFreezeIOVector =<< newListIOVector n es
{-# INLINE listVector #-}

replaceVector :: Vector e -> [(Int,e)] -> Vector e
replaceVector (Vector x@(IOVector _ n _ _ _)) ies =
    let go y ((i,e):ies') = do
              when (i < 0 || i >= n) $ error $ printf
                  "(//) <vector of dim %d> [ ..., (%d,_), ... ]: invalid index"
                  n i
              unsafeWriteElemIOVector y i e
              go y ies'
        go _ [] = return ()
    in unsafePerformIO $ do
        y  <- newCopyIOVector x
        go y ies
        unsafeFreezeIOVector y
{-# INLINE replaceVector #-}

unsafeReplaceVector :: Vector e -> [(Int,e)] -> Vector e
unsafeReplaceVector (Vector x@(IOVector _ _ _ _ _)) ies =
    let go y ((i,e):ies') = do
              unsafeWriteElemIOVector y i e
              go y ies'
        go _ [] = return ()
    in unsafePerformIO $ do
        y  <- newCopyIOVector x
        go y ies
        unsafeFreezeIOVector y
{-# INLINE unsafeReplaceVector #-}

-- | @zeroVector n@ creates a vector of dimension @n@ with all values
-- set to zero.
zeroVector :: (Elem e) => Int -> Vector e
zeroVector n = unsafePerformIO $ 
    unsafeFreezeIOVector =<< newZeroVector n
{-# INLINE zeroVector #-}

-- | @constantVector n e@ creates a vector of dimension @n@ with all values
-- set to @e@.
constantVector :: (Elem e) => Int -> e -> Vector e
constantVector n e = unsafePerformIO $
    unsafeFreezeIOVector =<< newConstantVector n e
{-# INLINE constantVector #-}

-- | @basisVector n i@ creates a vector of dimension @n@ with zeros 
-- everywhere but position @i@, where there is a one.
basisVector :: (Elem e) => Int -> Int -> Vector e
basisVector n i = unsafePerformIO $
    unsafeFreezeIOVector =<< newBasisVector n i
{-# INLINE basisVector #-}

-- | @subvector x o n@ creates a subvector of @x@ starting at index @o@ 
-- and having length @n@.
subvector :: Vector e -> Int -> Int -> Vector e
subvector = subvectorView
{-# INLINE subvector #-}

unsafeSubvector :: Vector e -> Int -> Int -> Vector e
unsafeSubvector = unsafeSubvectorView
{-# INLINE unsafeSubvector #-}

unsafeSubvectorWithStride :: 
    Int -> Vector e -> Int -> Int -> Vector e
unsafeSubvectorWithStride = unsafeSubvectorViewWithStride
{-# INLINE unsafeSubvectorWithStride #-}

-- | @subvectorWithStride s x o n@ creates a subvector of @x@ starting 
-- at index @o@, having length @n@ and stride @s@.
subvectorWithStride ::
    Int -> Vector e -> Int -> Int -> Vector e
subvectorWithStride = subvectorViewWithStride
{-# INLINE subvectorWithStride #-}
    
sizeVector :: Vector e -> Int
sizeVector (Vector x) = sizeIOVector x
{-# INLINE sizeVector #-}

indicesVector :: Vector e -> [Int]
indicesVector (Vector x) = indicesIOVector x
{-# INLINE indicesVector #-}

elemsVector :: Vector e -> [e]
elemsVector x@(Vector (IOVector _ _ _ _ _))
    | isConj x  = (map conjugate . elemsVector . conj) x
    | otherwise = case x of { (Vector (IOVector _ n f p incX)) ->
    let end = p `advancePtr` (n*incX)
        go p' | p' == end = inlinePerformIO $ do
                                touchForeignPtr f
                                return []
              | otherwise = let e  = inlinePerformIO (peek p')
                                es = go (p' `advancePtr` incX)
                            in e `seq` (e:es)
    in go p }
{-# SPECIALIZE INLINE elemsVector :: Vector Double -> [Double] #-}
{-# SPECIALIZE INLINE elemsVector :: Vector (Complex Double) -> [Complex Double] #-}

assocsVector :: Vector e -> [(Int,e)]
assocsVector x = zip (indicesVector x) (elemsVector x)
{-# INLINE assocsVector #-}

unsafeAtVector :: Vector e -> Int -> e
unsafeAtVector (Vector (IOVector c _ f p inc)) i
    | c == Conj = inlinePerformIO $ do
        e <- peekElemOff p (i*inc)
        touchForeignPtr f
        return $! conjugate e
    | otherwise = inlinePerformIO $ do
        e  <- peekElemOff p (i*inc)
        touchForeignPtr f
        return $! e
{-# INLINE unsafeAtVector #-}

tmapVector :: (e -> e) -> Vector e -> Vector e
tmapVector f x@(Vector (IOVector _ _ _ _ _)) = 
    listVector (dim x) (map f $ elemsVector x)
{-# INLINE tmapVector #-}

tzipWithVector :: 
    (e -> e -> e) -> Vector e -> Vector e -> Vector e
tzipWithVector f x@(Vector (IOVector _ _ _ _ _)) y
    | dim y /= n =
        error ("tzipWith: vector lengths differ; first has length `" ++
                show n ++ "' and second has length `" ++
                show (dim y) ++ "'")
    | otherwise =
        listVector n (zipWith f (elems x) (elems y))
  where
    n = dim x
{-# INLINE tzipWithVector #-}

-- | Compute the sum of absolute values of entries in the vector.
sumAbs :: (BLAS1 e) => Vector e -> Double
sumAbs (Vector x) = unsafePerformIO $ getSumAbs x
{-# INLINE sumAbs #-}

-- | Compute the 2-norm of a vector.
norm2 :: (BLAS1 e) => Vector e -> Double
norm2 (Vector x) = unsafePerformIO $ getNorm2 x
{-# INLINE norm2 #-}

-- | Get the index and norm of the element with absulte value.  Not valid 
-- if any of the vector entries are @NaN@.  Raises an exception if the 
-- vector has length @0@.
whichMaxAbs :: (BLAS1 e) => Vector e -> (Int, e)
whichMaxAbs (Vector x) = unsafePerformIO $ getWhichMaxAbs x
{-# INLINE whichMaxAbs #-}

-- | Compute the dot product of two vectors.
(<.>) :: (BLAS1 e) => Vector e -> Vector e -> e
(<.>) x y = unsafePerformIO $ getDot x y
{-# INLINE (<.>) #-}

unsafeDot :: (BLAS1 e) => Vector e -> Vector e -> e
unsafeDot x y = unsafePerformIO $ unsafeGetDot x y
{-# INLINE unsafeDot #-}

instance Shaped Vector Int where
    shape (Vector x) = shapeIOVector x
    {-# INLINE shape #-}
    bounds (Vector x) = boundsIOVector x
    {-# INLINE bounds #-}

instance ITensor Vector Int where
    (//) = replaceVector
    {-# INLINE (//) #-}
    unsafeReplace = unsafeReplaceVector
    {-# INLINE unsafeReplace #-}
    unsafeAt = unsafeAtVector
    {-# INLINE unsafeAt #-}
    size = sizeVector
    {-# INLINE size #-}
    elems = elemsVector
    {-# INLINE elems #-}
    indices = indicesVector
    {-# INLINE indices #-}
    assocs = assocsVector
    {-# INLINE assocs #-}
    tmap = tmapVector
    {-# INLINE tmap #-}
    (*>) k x = unsafePerformIO $ unsafeFreezeIOVector =<< getScaledVector k x
    {-# INLINE (*>) #-}
    shift k x = unsafePerformIO $ unsafeFreezeIOVector =<< getShiftedVector k x
    {-# INLINE shift #-}

instance (Monad m) => ReadTensor Vector Int m where
    getSize = return . size
    {-# INLINE getSize #-}
    getAssocs = return . assocs
    {-# INLINE getAssocs #-}
    getIndices = return . indices
    {-# INLINE getIndices #-}
    getElems = return . elems
    {-# INLINE getElems #-}
    getAssocs' = return . assocs
    {-# INLINE getAssocs' #-}
    getIndices' = return . indices
    {-# INLINE getIndices' #-}
    getElems' = return . elems
    {-# INLINE getElems' #-}
    unsafeReadElem x i = return $ unsafeAt x i
    {-# INLINE unsafeReadElem #-}
    
instance BaseVector Vector where    
    dim (Vector x) = dimIOVector x
    {-# INLINE dim #-}
    stride (Vector x) = strideIOVector x
    {-# INLINE stride #-}
    conjEnum (Vector x) = conjEnumIOVector x
    {-# INLINE conjEnum #-}
    conj (Vector x) = (Vector (conjIOVector x))
    {-# INLINE conj #-}
    unsafeSubvectorViewWithStride s (Vector x) o n = 
        Vector (unsafeSubvectorViewWithStrideIOVector s x o n)
    {-# INLINE unsafeSubvectorViewWithStride #-}
    unsafeVectorToIOVector (Vector x) = x
    {-# INLINE unsafeVectorToIOVector #-}
    unsafeIOVectorToVector = Vector
    {-# INLINE unsafeIOVectorToVector #-}

-- The NOINLINE pragmas and the strictness annotations here are *really*
-- important.  Otherwise, the compiler might think that certain actions
-- don't need to be run.
instance (MonadInterleave m) => ReadVector Vector m where
    freezeVector (Vector x) = return $! unsafePerformIO $ freezeIOVector x
    {-# NOINLINE freezeVector #-}
    unsafeFreezeVector = return
    {-# INLINE unsafeFreezeVector #-}
    unsafePerformIOWithVector (Vector x) f = return $! unsafePerformIO $ f x
    {-# NOINLINE unsafePerformIOWithVector #-}

instance (Show e) => Show (Vector e) where
    show x | isConj x  = "conj (" ++ show (conj x) ++ ")"
           | otherwise = "listVector " ++ show (dim x) ++ " " ++ show (elemsVector x)

instance (Eq e) => Eq (Vector e) where
    (==) = compareVectorWith (==)
    {-# INLINE (==) #-}

instance (AEq e) => AEq (Vector e) where
    (===) = compareVectorWith (===)
    {-# INLINE (===) #-}
    (~==) = compareVectorWith (~==)
    {-# INLINE (~==) #-}

compareVectorWith :: (e -> e -> Bool) -> Vector e -> Vector e -> Bool
compareVectorWith cmp x y
    | isConj x && isConj y =
        compareVectorWith cmp (conj x) (conj y)
    | otherwise =
        (dim x == dim y) && (and $ zipWith cmp (elemsVector x) (elemsVector y))
{-# INLINE compareVectorWith #-}

instance (BLAS1 e) => Num (Vector e) where
    (+) x y = unsafePerformIO $ unsafeFreezeIOVector =<< getAddVector x y
    {-# INLINE (+) #-}
    (-) x y = unsafePerformIO $ unsafeFreezeIOVector =<< getSubVector x y
    {-# INLINE (-) #-}
    (*) x y = unsafePerformIO $ unsafeFreezeIOVector =<< getMulVector x y
    {-# INLINE (*) #-}
    negate = ((-1) *>)
    {-# INLINE negate #-}
    abs = tmap abs
    {-# INLINE abs #-}
    signum = tmap signum
    {-# INLINE signum #-}
    fromInteger n = listVector 1 [fromInteger n]
    {-# INLINE fromInteger #-}
    
instance (BLAS1 e) => Fractional (Vector e) where
    (/) x y = unsafePerformIO $ unsafeFreezeIOVector =<< getDivVector x y
    {-# INLINE (/) #-}
    recip          = tmap recip
    {-# INLINE recip #-}
    fromRational q = listVector 1 [fromRational q]
    {-# INLINE fromRational #-}
    
instance (BLAS1 e, Floating e) => Floating (Vector e) where
    pi    = listVector 1 [pi]
    exp   = tmap exp
    {-# INLINE exp #-}
    sqrt  = tmap sqrt
    {-# INLINE sqrt #-}
    log   = tmap log
    {-# INLINE log #-}
    (**)  = tzipWithVector (**)
    {-# INLINE (**) #-}
    sin   = tmap sin
    {-# INLINE sin #-}
    cos   = tmap cos
    {-# INLINE cos #-}
    tan   = tmap tan
    {-# INLINE tan #-}
    asin  = tmap asin
    {-# INLINE asin #-}
    acos  = tmap acos
    {-# INLINE acos #-}
    atan  = tmap atan
    {-# INLINE atan #-}
    sinh  = tmap sinh
    {-# INLINE sinh #-}
    cosh  = tmap cosh
    {-# INLINE cosh #-}
    tanh  = tmap tanh
    {-# INLINE tanh #-}
    asinh = tmap asinh
    {-# INLINE asinh #-}
    acosh = tmap acosh
    {-# INLINE acosh #-}
    atanh = tmap atanh
    {-# INLINE atanh #-}

vectorCall :: (ReadVector x m)
           => (Int -> Ptr e -> Int -> IO a) 
           ->  x e -> m a
vectorCall f x = 
    unsafePerformIOWithVector x $ \x' ->
        let n    = dimIOVector x'
            incX = strideIOVector x'
        in withIOVector x' $ \pX ->
               f n pX incX
{-# INLINE vectorCall #-}

vectorCall2 :: (ReadVector x m, ReadVector y m)
            => (Int -> Ptr e -> Int -> Ptr f -> Int -> IO a) 
            -> x e -> y f -> m a
vectorCall2 f x y =
    unsafePerformIOWithVector y $ \y' ->
        let x'   = unsafeVectorToIOVector x
            n    = dimIOVector x'
            incX = strideIOVector x'
            incY = strideIOVector y'
        in withIOVector x' $ \pX ->
           withIOVector y' $ \pY ->
               f n pX incX pY incY
{-# INLINE vectorCall2 #-}    

checkVectorOp2 :: (BaseVector x, BaseVector y) => 
    (x e -> y f -> a) ->
        x e -> y f -> a
checkVectorOp2 f x y = 
    checkBinaryOp (dim x) (dim y) $ f x y
{-# INLINE checkVectorOp2 #-}

getUnaryVectorOp :: (ReadVector x m, WriteVector y m) =>
    (y e -> m ()) -> x e -> m (y e)
getUnaryVectorOp f x = do
    y  <- newCopyVector x
    f y
    return y
{-# INLINE getUnaryVectorOp #-}

unsafeGetBinaryVectorOp :: 
    (WriteVector z m, ReadVector x m, ReadVector y m) =>
    (z e -> y e -> m ()) ->
        x e -> y e -> m (z e)
unsafeGetBinaryVectorOp f x y = do
    z  <- newCopyVector x
    f z y
    return z
{-# INLINE unsafeGetBinaryVectorOp #-}

-- | Dense vectors in the 'ST' monad.  The type arguments are as follows:
--
--     * @s@: the state variable argument for the 'ST' type
--
--     * @e@: the element type of the vector.  Only certain element types
--       are supported.
--
newtype STVector s e = STVector (IOVector e)

-- | A safe way to create and work with a mutable vector before returning 
-- an immutable vector for later perusal. This function avoids copying
-- the vector before returning it - it uses unsafeFreezeVector internally,
-- but this wrapper is a safe interface to that function. 
runSTVector :: (forall s . ST s (STVector s e)) -> Vector e
runSTVector mx = 
    runST $ mx >>= \(STVector x) -> return (Vector x)
{-# INLINE runSTVector #-}

instance Shaped (STVector s) Int where
    shape (STVector x) = shapeIOVector x
    {-# INLINE shape #-}
    bounds (STVector x) = boundsIOVector x
    {-# INLINE bounds #-}

instance ReadTensor (STVector s) Int (ST s) where
    getSize (STVector x) = unsafeIOToST $ getSizeIOVector x
    {-# INLINE getSize #-}
    unsafeReadElem (STVector x) i = unsafeIOToST $ unsafeReadElemIOVector x i
    {-# INLINE unsafeReadElem #-}
    getIndices (STVector x) = unsafeIOToST $ getIndicesIOVector x
    {-# INLINE getIndices #-}
    getIndices' (STVector x) = unsafeIOToST $ getIndicesIOVector' x
    {-# INLINE getIndices' #-}
    getElems (STVector x) = unsafeIOToST $ getElemsIOVector x
    {-# INLINE getElems #-}
    getElems' (STVector x) = unsafeIOToST $ getElemsIOVector' x
    {-# INLINE getElems' #-}
    getAssocs (STVector x) = unsafeIOToST $ getAssocsIOVector x
    {-# INLINE getAssocs #-}
    getAssocs' (STVector x) = unsafeIOToST $ getAssocsIOVector' x
    {-# INLINE getAssocs' #-}

instance WriteTensor (STVector s) Int (ST s) where
    getMaxSize (STVector x) = unsafeIOToST $ getMaxSizeIOVector x
    {-# INLINE getMaxSize #-}
    setZero (STVector x) = unsafeIOToST $ setZeroIOVector x
    {-# INLINE setZero #-}
    setConstant e (STVector x) = unsafeIOToST $ setConstantIOVector e x
    {-# INLINE setConstant #-}
    canModifyElem (STVector x) i = unsafeIOToST $ canModifyElemIOVector x i
    {-# INLINE canModifyElem #-}
    unsafeWriteElem (STVector x) i e= unsafeIOToST $ unsafeWriteElemIOVector x i e
    {-# INLINE unsafeWriteElem #-}
    unsafeModifyElem (STVector x) i f = unsafeIOToST $ unsafeModifyElemIOVector x i f
    {-# INLINE unsafeModifyElem #-}
    modifyWith f (STVector x) = unsafeIOToST $ modifyWithIOVector f x
    {-# INLINE modifyWith #-}
    doConj (STVector x) = unsafeIOToST $ doConjIOVector x
    {-# INLINE doConj #-}
    scaleBy k (STVector x) = unsafeIOToST $ scaleByIOVector k x
    {-# INLINE scaleBy #-}
    shiftBy k (STVector x) = unsafeIOToST $ shiftByIOVector k x
    {-# INLINE shiftBy #-}

instance BaseVector (STVector s) where
    dim (STVector x) = dimIOVector x
    {-# INLINE dim #-}
    stride (STVector x) = strideIOVector x
    {-# INLINE stride #-}
    conjEnum (STVector x) = conjEnumIOVector x
    {-# INLINE conjEnum #-}
    conj (STVector x) = STVector (conjIOVector x)
    {-# INLINE conj #-}
    unsafeSubvectorViewWithStride s (STVector x) o n = 
        STVector (unsafeSubvectorViewWithStrideIOVector s x o n)
    {-# INLINE unsafeSubvectorViewWithStride #-}    
    unsafeVectorToIOVector (STVector x) = x
    {-# INLINE unsafeVectorToIOVector #-}
    unsafeIOVectorToVector = STVector
    {-# INLINE unsafeIOVectorToVector #-}

instance ReadVector (STVector s) (ST s) where
    freezeVector (STVector x) = unsafeIOToST $ freezeIOVector x
    {-# INLINE freezeVector #-}
    unsafeFreezeVector (STVector x) = unsafeIOToST $ unsafeFreezeIOVector x
    {-# INLINE unsafeFreezeVector #-}
    unsafePerformIOWithVector (STVector x) f = unsafeIOToST $ f x
    {-# INLINE unsafePerformIOWithVector #-}

instance WriteVector (STVector s) (ST s) where
    newVector_ = liftM STVector . unsafeIOToST . newIOVector_
    {-# INLINE newVector_ #-}
    thawVector = liftM STVector . unsafeIOToST . thawIOVector
    {-# INLINE thawVector #-}
    unsafeThawVector = liftM STVector . unsafeIOToST . unsafeThawIOVector
    {-# INLINE unsafeThawVector #-}
    unsafeConvertIOVector = unsafeIOToST . liftM STVector
    {-# INLINE unsafeConvertIOVector #-}
