{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances #-}
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
import Data.AEq
import Foreign
import Text.Printf
import Unsafe.Coerce

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
--     * @n@: a phantom type for the dimension of the vector
--
--     * @e@: the element type of the vector.  Only certain element types
--       are supported.
--
newtype Vector n e = Vector (IOVector n e)

freezeIOVector :: IOVector n e -> IO (Vector n e)
freezeIOVector = liftM Vector . newCopyIOVector
{-# INLINE freezeIOVector #-}

thawIOVector :: Vector n e -> IO (IOVector n e)
thawIOVector (Vector x) = newCopyIOVector x
{-# INLINE thawIOVector #-}

unsafeFreezeIOVector :: IOVector n e -> IO (Vector n e)
unsafeFreezeIOVector = return . Vector
{-# INLINE unsafeFreezeIOVector #-}

unsafeThawIOVector :: Vector n e -> IO (IOVector n e)
unsafeThawIOVector (Vector x) = return x
{-# INLINE unsafeThawIOVector #-}

-- | Common functionality for all vector types.
class (Shaped x Int) => BaseVector x where
    -- | Get the dimension (length) of the vector.
    dim :: x n e -> Int
    
    -- | Get the memory stride (in elements) between consecutive elements.
    stride :: x n e -> Int

    -- | Indicate whether or not internally the vector stores the complex
    -- conjugates of its elements.
    isConj :: x n e -> Bool
    isConj x = conjEnum x == Conj
    {-# INLINE isConj #-}

    -- | Get the storage type.
    conjEnum :: x n e -> ConjEnum
    
    -- | Get a view into the complex conjugate of a vector.
    conj :: x n e -> x n e

    -- | Cast the shape type of the vector.
    coerceVector :: x n e -> x n' e
    coerceVector = unsafeCoerce
    {-# INLINE coerceVector #-}

    unsafeSubvectorViewWithStride :: Int -> x n e -> Int -> Int -> x n' e

    -- | Unsafe cast from a vector to an 'IOVector'.
    unsafeVectorToIOVector :: x n e -> IOVector n e
    unsafeIOVectorToVector :: IOVector n e -> x n e


-- | Vectors that can be read in a monad.
class (BaseVector x, Monad m, ReadTensor x Int m) => ReadVector x m where
    -- | Cast the vector to an 'IOVector', perform an @IO@ action, and
    -- convert the @IO@ action to an action in the monad @m@.  This
    -- operation is /very/ unsafe.
    unsafePerformIOWithVector :: x n e -> (IOVector n e -> IO a) -> m a

    -- | Convert a mutable vector to an immutable one by taking a complete
    -- copy of it.
    freezeVector :: x n e -> m (Vector n e)
    unsafeFreezeVector :: x n e -> m (Vector n e)


-- | Vectors that can be created or modified in a monad.
class (ReadVector x m, WriteTensor x Int m) => WriteVector x m where
    -- | Creates a new vector of the given length.  The elements will be
    -- uninitialized.
    newVector_ :: (Elem e) => Int -> m (x n e)

    -- | Convert an immutable vector to a mutable one by taking a complete
    -- copy of it.
    thawVector :: Vector n e -> m (x n e)
    unsafeThawVector :: Vector n e -> m (x n e)

    -- | Unsafely convert an 'IO' action that creates an 'IOVector' into
    -- an action in @m@ that creates a vector.
    unsafeConvertIOVector :: IO (IOVector n e) -> m (x n e)


-- | Creates a new vector with the given association list.  Unspecified
-- indices will get initialized to zero.
newVector :: (WriteVector x m, Elem e) => Int -> [(Int,e)] -> m (x n e)
newVector n ies = 
    unsafeConvertIOVector $ newIOVector n ies
{-# INLINE newVector #-}

unsafeNewVector :: (WriteVector x m, Elem e) => Int -> [(Int,e)] -> m (x n e)
unsafeNewVector n ies = 
    unsafeConvertIOVector $ unsafeNewIOVector n ies
{-# INLINE unsafeNewVector #-}

-- | Creates a new vector of the given dimension with the given elements.
-- If the list has length less than the passed-in dimenson, the tail of
-- the vector will be uninitialized.
newListVector :: (WriteVector x m, Elem e) => Int -> [e] -> m (x n e)
newListVector n es = 
    unsafeConvertIOVector $ newListIOVector n es
{-# INLINE newListVector #-}

-- | Create a zero vector of the specified length.
newZeroVector :: (WriteVector x m, Elem e) => Int -> m (x n e)
newZeroVector n =
    unsafeConvertIOVector $ newZeroIOVector n
{-# INLINE newZeroVector #-}

-- | Set every element in the vector to zero.
setZeroVector :: (WriteVector x m) => x n e -> m ()
setZeroVector x =
    unsafePerformIOWithVector x $ setZeroIOVector
{-# INLINE setZeroVector #-}

-- | Create a vector with every element initialized to the same value.
newConstantVector :: (WriteVector x m, Elem e) => Int -> e -> m (x n e)
newConstantVector n e = 
    unsafeConvertIOVector $ newConstantIOVector n e
{-# INLINE newConstantVector #-}

-- | Set every element in the vector to a constant.
setConstantVector :: (WriteVector x m) => e -> x n e -> m ()
setConstantVector e x =
    unsafePerformIOWithVector x $ setConstantIOVector e
{-# INLINE setConstantVector #-}

-- | @newBasisVector n i@ creates a vector of length @n@ that is all zero 
-- except for at position @i@, where it equal to one.
newBasisVector :: (WriteVector x m, Elem e) => Int -> Int -> m (x n e)
newBasisVector n i = unsafeConvertIOVector $
    newBasisIOVector n i
{-# INLINE newBasisVector #-}

-- | @setBasis x i@ sets the @i@th coordinate of @x@ to @1@, and all other
-- coordinates to @0@.  If the vector has been scaled, it is possible that
-- @readVector x i@ will not return exactly @1@.  See 'setElem'.
setBasisVector :: (WriteVector x m) => Int -> x n e -> m ()
setBasisVector i x =
    unsafePerformIOWithVector x $ setBasisIOVector i
{-# INLINE setBasisVector #-}

-- | Creates a new vector by copying another one.
newCopyVector :: (ReadVector x m, WriteVector y m) =>
    x n e -> m (y n e)
newCopyVector x = unsafeConvertIOVector $
    newCopyIOVector (unsafeVectorToIOVector x)
{-# INLINE newCopyVector #-}

-- | Creates a new vector by copying another one.  The returned vector
-- is gauranteed not to be a view into another vector.  That is, the
-- returned vector will have @isConj@ to be @False@.
newCopyVector' :: (ReadVector x m, WriteVector y m) => x n e -> m (y n e)
newCopyVector' x = unsafeConvertIOVector $
    newCopyIOVector' (unsafeVectorToIOVector x)
{-# INLINE newCopyVector' #-}

-- | @copyVector dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyVector :: (WriteVector y m, ReadVector x m) =>
    y n e -> x n e -> m ()
copyVector y x = checkBinaryOp (shape x) (shape y) $ unsafeCopyVector y x
{-# INLINE copyVector #-}

unsafeCopyVector :: (WriteVector y m, ReadVector x m) =>
    y n e -> x n e -> m ()
unsafeCopyVector y x =
    unsafePerformIOWithVector y $
        (`unsafeCopyIOVector` (unsafeVectorToIOVector x))
{-# INLINE unsafeCopyVector #-}

-- | Swap the values stored in two vectors.
swapVector :: (WriteVector x m, WriteVector y m) => 
    x n e -> y n e -> m ()
swapVector x y = checkBinaryOp (shape x) (shape y) $ unsafeSwapVector x y
{-# INLINE swapVector #-}

unsafeSwapVector :: (WriteVector x m, WriteVector y m) =>
    x n e -> y n e -> m ()
unsafeSwapVector x y =
    unsafePerformIOWithVector x $ 
        (`unsafeSwapIOVector` (unsafeVectorToIOVector y))
{-# INLINE unsafeSwapVector #-}

-- | @subvectorView x o n@ creates a subvector view of @x@ starting at index @o@ 
-- and having length @n@.
subvectorView :: (BaseVector x) => 
    x n e -> Int -> Int -> x n' e
subvectorView x = checkedSubvector (dim x) (unsafeSubvectorView x)
{-# INLINE subvectorView #-}

unsafeSubvectorView :: (BaseVector x) => 
    x n e -> Int -> Int -> x n' e
unsafeSubvectorView = unsafeSubvectorViewWithStride 1
{-# INLINE unsafeSubvectorView #-}

-- | @subvectorViewWithStride s x o n@ creates a subvector view of @x@ starting 
-- at index @o@, having length @n@ and stride @s@.
subvectorViewWithStride :: (BaseVector x) => 
    Int -> x n e -> Int -> Int -> x n' e
subvectorViewWithStride s x = 
    checkedSubvectorWithStride s (dim x) (unsafeSubvectorViewWithStride s x)
{-# INLINE subvectorViewWithStride #-}

-- | Get a new vector with elements with the conjugates of the elements
-- of the given vector
getConjVector :: (ReadVector x m, WriteVector y m) =>
    x n e -> m (y n e)
getConjVector = getUnaryVectorOp doConjVector
{-# INLINE getConjVector #-}

-- | Conjugate every element of the vector.
doConjVector :: (WriteVector y m) => y n e -> m ()
doConjVector x =
    unsafePerformIOWithVector x $ doConjIOVector
{-# INLINE doConjVector #-}

-- | Get a new vector by scaling the elements of another vector
-- by a given value.
getScaledVector :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> x n e -> m (y n e)
getScaledVector e = getUnaryVectorOp (scaleByVector e)
{-# INLINE getScaledVector #-}

-- | Scale every element by the given value.
scaleByVector :: (WriteVector y m, BLAS1 e) => e -> y n e -> m ()
scaleByVector k x =
    unsafePerformIOWithVector x $ scaleByIOVector k
{-# INLINE scaleByVector #-}

-- | Get a new vector by shifting the elements of another vector
-- by a given value.
getShiftedVector :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> x n e -> m (y n e)
getShiftedVector e = getUnaryVectorOp (shiftByVector e)
{-# INLINE getShiftedVector #-}

-- | Add the given value to every element.
shiftByVector :: (WriteVector y m, BLAS1 e) => e -> y n e -> m ()
shiftByVector k x =
    unsafePerformIOWithVector x $ shiftByIOVector k
{-# INLINE shiftByVector #-}

-- | @getAddVector x y@ creates a new vector equal to the sum @x+y@.  The 
-- operands must have the same dimension.
getAddVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x n e -> y n e -> m (z n e)
getAddVector = checkVectorOp2 unsafeGetAddVector
{-# INLINE getAddVector #-}

unsafeGetAddVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x n e -> y n e -> m (z n e)
unsafeGetAddVector = unsafeGetBinaryVectorOp unsafeAddVector
{-# INLINE unsafeGetAddVector #-}

-- | @addVector y x@ replaces @y@ with @y+x@.
addVector :: (WriteVector y m, ReadVector x m, BLAS1 e) => 
    y n e -> x n e -> m ()
addVector y x = checkBinaryOp (dim y) (dim x) $ unsafeAddVector y x
{-# INLINE addVector #-}

unsafeAddVector :: (WriteVector y m, ReadVector x m, BLAS1 e) =>
    y n e -> x n e -> m ()
unsafeAddVector y x = unsafeAxpyVector 1 x y
{-# INLINE unsafeAddVector #-}

-- | @getSubVector x y@ creates a new tensor equal to the difference @x-y@.  
-- The operands must have the same dimension.
getSubVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x n e -> y n e -> m (z n e)
getSubVector = checkVectorOp2 unsafeGetSubVector
{-# INLINE getSubVector #-}

unsafeGetSubVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x n e -> y n e -> m (z n e)
unsafeGetSubVector = unsafeGetBinaryVectorOp unsafeSubVector
{-# INLINE unsafeGetSubVector #-}

-- | @subVector y x@ replaces @y@ with @y-x@.
subVector :: (WriteVector y m, ReadVector x m, BLAS1 e) => 
    y n e -> x n e -> m ()
subVector y x = checkBinaryOp (dim y) (dim x) $ unsafeSubVector y x
{-# INLINE subVector #-}

unsafeSubVector :: (WriteVector y m, ReadVector x m, BLAS1 e) =>
    y n e -> x n e -> m ()
unsafeSubVector y x = unsafeAxpyVector (-1) x y
{-# INLINE unsafeSubVector #-}

-- | @axpyVector alpha x y@ replaces @y@ with @alpha * x + y@.
axpyVector :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> x n e -> y n e -> m ()
axpyVector alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyVector alpha x y
{-# INLINE axpyVector #-}

unsafeAxpyVector :: (ReadVector x m, ReadVector y m, BLAS1 e) =>
    e -> x n e -> y n e -> m ()
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
    x n e -> y n e -> m (z n e)
getMulVector = checkVectorOp2 unsafeGetMulVector
{-# INLINE getMulVector #-}

unsafeGetMulVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x n e -> y n e -> m (z n e)
unsafeGetMulVector = unsafeGetBinaryVectorOp unsafeMulVector
{-# INLINE unsafeGetMulVector #-}

-- | @mulVector y x@ replaces @y@ with @y*x@.
mulVector :: (WriteVector y m, ReadVector x m, BLAS1 e) => 
    y n e -> x n e -> m ()
mulVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeMulVector y x
{-# INLINE mulVector #-}
 
unsafeMulVector :: (WriteVector y m, ReadVector x m, BLAS1 e) =>
    y n e -> x n e -> m ()
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
    x n e -> y n e -> m (z n e)
getDivVector = checkVectorOp2 unsafeGetDivVector
{-# INLINE getDivVector #-}

unsafeGetDivVector :: 
    (ReadVector x m, ReadVector y m, WriteVector z m, BLAS1 e) =>
    x n e -> y n e -> m (z n e)
unsafeGetDivVector = unsafeGetBinaryVectorOp unsafeDivVector
{-# INLINE unsafeGetDivVector #-}

-- | @divVector y x@ replaces @y@ with @y/x@.
divVector :: (WriteVector y m, ReadVector x m, BLAS1 e) => 
    y n e -> x n e -> m ()
divVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeDivVector y x
{-# INLINE divVector #-}

unsafeDivVector :: (WriteVector y m, ReadVector x m, BLAS1 e) =>
    y n e -> x n e -> m ()
unsafeDivVector y x
    | isConj y =
        unsafeDivVector (conj y) (conj x)
    | isConj x =
        vectorCall2 BLAS.vcdiv x y
    | otherwise =
        vectorCall2 BLAS.vdiv x y
{-# INLINE unsafeDivVector #-}

-- | Gets the sum of the absolute values of the vector entries.
getSumAbs :: (ReadVector x m, BLAS1 e) => x n e -> m Double
getSumAbs = vectorCall BLAS.asum
{-# INLINE getSumAbs #-}
    
-- | Gets the 2-norm of a vector.
getNorm2 :: (ReadVector x m, BLAS1 e) => x n e -> m Double
getNorm2 = vectorCall BLAS.nrm2
{-# INLINE getNorm2 #-}

-- | Gets the index and norm of the element with maximum magnitude.  This is 
-- undefined if any of the elements are @NaN@.  It will throw an exception if 
-- the dimension of the vector is 0.
getWhichMaxAbs :: (ReadVector x m, BLAS1 e) => x n e -> m (Int, e)
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
    x n e -> y n e -> m e
getDot x y = checkVecVecOp "getDot" (dim x) (dim y) $ unsafeGetDot x y
{-# INLINE getDot #-}

unsafeGetDot :: (ReadVector x m, ReadVector y m, BLAS1 e) => 
    x n e -> y n e -> m e
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
    unsafePerformIOWithVector x f = f x
    {-# INLINE unsafePerformIOWithVector #-}
    freezeVector = freezeIOVector
    {-# INLINE freezeVector #-}
    unsafeFreezeVector = unsafeFreezeIOVector
    {-# INLINE unsafeFreezeVector #-}

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
vector :: (Elem e) => Int -> [(Int, e)] -> Vector n e
vector n ies = unsafePerformIO $
    unsafeFreezeIOVector =<< newIOVector n ies
{-# INLINE vector #-}

-- Same as 'vector', but does not range-check the indices.
unsafeVector :: (Elem e) => Int -> [(Int, e)] -> Vector n e
unsafeVector n ies = unsafePerformIO $
    unsafeFreezeIOVector =<< unsafeNewIOVector n ies
{-# INLINE unsafeVector #-}

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.  @listVector n es@ is equivalent to 
-- @vector n (zip [0..(n-1)] es)@, except that the result is undefined 
-- if @length es@ is less than @n@.
listVector :: (Elem e) => Int -> [e] -> Vector n e
listVector n es = unsafePerformIO $
    unsafeFreezeIOVector =<< newListIOVector n es
{-# INLINE listVector #-}

replaceVector :: Vector n e -> [(Int,e)] -> Vector n e
replaceVector (Vector x@(IOVector _ n _ _ _)) ies =
    let go y ((i,e):ies') = do
              when (i < 0 || i >= n) $ error $ printf
                  "(//) <vector of dim %d> [ ..., (%d,_), ... ]: invalid index"
                  n i
              io <- unsafeWriteElemIOVector y i e
              io `seq` go y ies'
        go _ [] = return ()
    in unsafePerformIO $ do
        y  <- newCopyIOVector x
        io <- go y ies
        io `seq` unsafeFreezeIOVector y
{-# INLINE replaceVector #-}

unsafeReplaceVector :: Vector n e -> [(Int,e)] -> Vector n e
unsafeReplaceVector (Vector x@(IOVector _ _ _ _ _)) ies =
    let go y ((i,e):ies') = do
              io <- unsafeWriteElemIOVector y i e
              io `seq` go y ies'
        go _ [] = return ()
    in unsafePerformIO $ do
        y  <- newCopyIOVector x
        io <- go y ies
        io `seq` unsafeFreezeIOVector y
{-# INLINE unsafeReplaceVector #-}

-- | @zeroVector n@ creates a vector of dimension @n@ with all values
-- set to zero.
zeroVector :: (Elem e) => Int -> Vector n e
zeroVector n = unsafePerformIO $ 
    unsafeFreezeIOVector =<< newZeroVector n
{-# INLINE zeroVector #-}

-- | @constantVector n e@ creates a vector of dimension @n@ with all values
-- set to @e@.
constantVector :: (Elem e) => Int -> e -> Vector n e
constantVector n e = unsafePerformIO $
    unsafeFreezeIOVector =<< newConstantVector n e
{-# INLINE constantVector #-}

-- | @basisVector n i@ creates a vector of dimension @n@ with zeros 
-- everywhere but position @i@, where there is a one.
basisVector :: (Elem e) => Int -> Int -> Vector n e
basisVector n i = unsafePerformIO $
    unsafeFreezeIOVector =<< newBasisVector n i
{-# INLINE basisVector #-}

-- | @subvector x o n@ creates a subvector of @x@ starting at index @o@ 
-- and having length @n@.
subvector :: Vector n e -> Int -> Int -> Vector n' e
subvector = subvectorView
{-# INLINE subvector #-}

unsafeSubvector :: Vector n e -> Int -> Int -> Vector n' e
unsafeSubvector = unsafeSubvectorView
{-# INLINE unsafeSubvector #-}

unsafeSubvectorWithStride :: 
    Int -> Vector n e -> Int -> Int -> Vector n' e
unsafeSubvectorWithStride = unsafeSubvectorViewWithStride
{-# INLINE unsafeSubvectorWithStride #-}

-- | @subvectorWithStride s x o n@ creates a subvector of @x@ starting 
-- at index @o@, having length @n@ and stride @s@.
subvectorWithStride ::
    Int -> Vector n e -> Int -> Int -> Vector n' e
subvectorWithStride = subvectorViewWithStride
{-# INLINE subvectorWithStride #-}
    
sizeVector :: Vector n e -> Int
sizeVector (Vector x) = sizeIOVector x
{-# INLINE sizeVector #-}

indicesVector :: Vector n e -> [Int]
indicesVector (Vector x) = indicesIOVector x
{-# INLINE indicesVector #-}

elemsVector :: Vector n e -> [e]
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
{-# SPECIALIZE INLINE elemsVector :: Vector n Double -> [Double] #-}
{-# SPECIALIZE INLINE elemsVector :: Vector n (Complex Double) -> [Complex Double] #-}

assocsVector :: Vector n e -> [(Int,e)]
assocsVector x = zip (indicesVector x) (elemsVector x)
{-# INLINE assocsVector #-}

unsafeAtVector :: Vector n e -> Int -> e
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

tmapVector :: (e -> e) -> Vector n e -> Vector n e
tmapVector f x@(Vector (IOVector _ _ _ _ _)) = 
    listVector (dim x) (map f $ elemsVector x)
{-# INLINE tmapVector #-}

tzipWithVector :: 
    (e -> e -> e) -> Vector n e -> Vector n e -> Vector n e
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
sumAbs :: (BLAS1 e) => Vector n e -> Double
sumAbs (Vector x) = unsafePerformIO $ getSumAbs x
{-# INLINE sumAbs #-}

-- | Compute the 2-norm of a vector.
norm2 :: (BLAS1 e) => Vector n e -> Double
norm2 (Vector x) = unsafePerformIO $ getNorm2 x
{-# INLINE norm2 #-}

-- | Get the index and norm of the element with absulte value.  Not valid 
-- if any of the vector entries are @NaN@.  Raises an exception if the 
-- vector has length @0@.
whichMaxAbs :: (BLAS1 e) => Vector n e -> (Int, e)
whichMaxAbs (Vector x) = unsafePerformIO $ getWhichMaxAbs x
{-# INLINE whichMaxAbs #-}

-- | Compute the dot product of two vectors.
(<.>) :: (BLAS1 e) => Vector n e -> Vector n e -> e
(<.>) x y = unsafePerformIO $ getDot x y
{-# INLINE (<.>) #-}

unsafeDot :: (BLAS1 e) => Vector n e -> Vector n e -> e
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

instance (Monad m) => ReadVector Vector m where
    unsafePerformIOWithVector (Vector x) f = (return . unsafePerformIO) $ do
        r <- f x
        return $! r
    {-# INLINE unsafePerformIOWithVector #-}    
    freezeVector (Vector x) = (return . unsafePerformIO) $ freezeIOVector x
    {-# INLINE freezeVector #-}
    unsafeFreezeVector = return
    {-# INLINE unsafeFreezeVector #-}

instance (Show e) => Show (Vector n e) where
    show x | isConj x  = "conj (" ++ show (conj x) ++ ")"
           | otherwise = "listVector " ++ show (dim x) ++ " " ++ show (elemsVector x)

instance (Eq e) => Eq (Vector n e) where
    (==) = compareVectorWith (==)
    {-# INLINE (==) #-}

instance (AEq e) => AEq (Vector n e) where
    (===) = compareVectorWith (===)
    {-# INLINE (===) #-}
    (~==) = compareVectorWith (~==)
    {-# INLINE (~==) #-}

compareVectorWith :: (e -> e -> Bool) -> Vector n e -> Vector n e -> Bool
compareVectorWith cmp x y
    | isConj x && isConj y =
        compareVectorWith cmp (conj x) (conj y)
    | otherwise =
        (dim x == dim y) && (and $ zipWith cmp (elemsVector x) (elemsVector y))
{-# INLINE compareVectorWith #-}

instance (BLAS1 e) => Num (Vector n e) where
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
    
instance (BLAS1 e) => Fractional (Vector n e) where
    (/) x y = unsafePerformIO $ unsafeFreezeIOVector =<< getDivVector x y
    {-# INLINE (/) #-}
    recip          = tmap recip
    {-# INLINE recip #-}
    fromRational q = listVector 1 [fromRational q]
    {-# INLINE fromRational #-}
    
instance (BLAS1 e, Floating e) => Floating (Vector n e) where
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
           ->  x n e -> m a
vectorCall f x = 
    unsafePerformIOWithVector x $ \x' ->
        let n    = dimIOVector x'
            incX = strideIOVector x'
        in withIOVector x' $ \pX ->
               f n pX incX
{-# INLINE vectorCall #-}

vectorCall2 :: (ReadVector x m, ReadVector y m)
            => (Int -> Ptr e -> Int -> Ptr f -> Int -> IO a) 
            -> x n e -> y n' f -> m a
vectorCall2 f x y =
    unsafePerformIOWithVector x $ \x' ->
        let y'   = unsafeVectorToIOVector y
            n    = dimIOVector x'
            incX = strideIOVector x'
            incY = strideIOVector y'
        in withIOVector x' $ \pX ->
           withIOVector y' $ \pY ->
               f n pX incX pY incY
{-# INLINE vectorCall2 #-}    

checkVectorOp2 :: (BaseVector x, BaseVector y) => 
    (x n e -> y n f -> a) ->
        x n e -> y n f -> a
checkVectorOp2 f x y = 
    checkBinaryOp (dim x) (dim y) $ f x y
{-# INLINE checkVectorOp2 #-}

getUnaryVectorOp :: (ReadVector x m, WriteVector y m) =>
    (y n e -> m ()) -> x n e -> m (y n e)
getUnaryVectorOp f x = do
    y  <- newCopyVector x
    io <- f y
    io `seq` return y
{-# INLINE getUnaryVectorOp #-}

unsafeGetBinaryVectorOp :: 
    (WriteVector z m, ReadVector x m, ReadVector y m) =>
    (z n e -> y n e -> m ()) ->
        x n e -> y n e -> m (z n e)
unsafeGetBinaryVectorOp f x y = do
    z  <- newCopyVector x
    io <- f z y
    io `seq` return z
{-# INLINE unsafeGetBinaryVectorOp #-}
