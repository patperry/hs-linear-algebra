{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances,
        FunctionalDependencies #-}
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
import Control.Monad.ST
import Data.AEq
import Data.Complex( Complex )
import Foreign
import Unsafe.Coerce

import BLAS.Internal( checkBinaryOp, clearArray, inlinePerformIO )

import Data.Elem.BLAS ( Elem, BLAS1, conjugate )
import qualified Data.Elem.BLAS as BLAS

import Data.Tensor.Class
import Data.Tensor.Class.ITensor
import Data.Tensor.Class.MTensor

import Data.Vector.Dense.IOBase

-- | The immutable dense vector data type.
newtype Vector n e = Vector (IOVector n e)

unsafeFreezeIOVector :: IOVector n e -> IO (Vector n e)
unsafeFreezeIOVector = return . Vector

unsafeThawIOVector :: Vector n e -> IO (IOVector n e)
unsafeThawIOVector (Vector x) = return x

newtype STVector s n e = STVector (IOVector n e)

-- | Common functionality between all vector types.
class (Shaped x Int e, Elem e) => BaseVector x e where
    -- | Get the dimension (length) of the vector.
    dim :: x n e -> Int
    
    -- | Get the memory stride (in elements) between consecutive elements.
    stride :: x n e -> Int

    -- | Indicate whether or not the vector stores the complex conjugates
    -- of its elements.
    isConj :: x n e -> Bool

    -- | Get a view into the complex conjugate of a vector.
    conj :: x n e -> x n e

    -- | Execute an @IO@ action with a pointer to the first element in the
    -- vector.
    withVectorPtrIO :: x n e -> (Ptr e -> IO a) -> IO a

-- | An overloaded interface for vectors that can be read or created in
-- a monad.
class (BaseVector x e, BLAS1 e, Monad m, ReadTensor x Int e m) => ReadVector x e m where
    -- | Execture an @IO@ action with a pointer to the first element and then
    -- convert the @IO@ action to an action in the monad @m@.
    withVectorPtr :: x n e -> (Ptr e -> IO a) -> m a

    unsafeFreezeVector :: x n e -> m (Vector n e)

-- | An overloaded interface for vectors that can be modified in a monad.
class (ReadVector x e m, WriteTensor x Int e m) => WriteVector x e m | x -> m, m -> x where
    -- | Creates a new vector of the given length.  The elements will be 
    -- uninitialized.
    newVector_ :: Int -> m (x n e)

    unsafeThawVector :: Vector n e -> m (x n e)

-- | Cast the shape type of the vector.
coerceVector :: (BaseVector x e) => x n e -> x n' e
coerceVector = unsafeCoerce
{-# INLINE coerceVector #-}

-- | Creates a new vector with the given association list.  Unspecified
-- indices will get initialized to zero.
newVector :: (WriteVector x e m) => Int -> [(Int,e)] -> m (x n e)
newVector n ies = do
    x <- newZeroVector n
    mapM_ (uncurry $ writeElem x) ies
    return x
{-# INLINE newVector #-}

unsafeNewVector :: (WriteVector x e m) => Int -> [(Int,e)] -> m (x n e)
unsafeNewVector n ies = do
    x <- newZeroVector n
    mapM_ (uncurry $ unsafeWriteElem x) ies
    return x
{-# INLINE unsafeNewVector #-}

-- | Create a zero vector of the specified length.
newZeroVector :: (WriteVector x e m) => Int -> m (x n e)
newZeroVector n = do
    x <- newVector_ n
    withVectorPtr x $ \p ->
        clearArray p n
    return x
{-# INLINE newZeroVector #-}

-- | Create a vector with every element initialized to the same value.
newConstantVector :: (WriteVector x e m) => Int -> e -> m (x n e)
newConstantVector n e = do
    x <- newVector_ n
    withVectorPtr x $ \p ->
        pokeArray p (replicate n e)
    return x
{-# INLINE newConstantVector #-}

-- | Creates a new vector of the given dimension with the given elements.
-- If the list has length less than the passed-in dimenson, the tail of
-- the vector will be uninitialized.
newListVector :: (WriteVector x e m) => Int -> [e] -> m (x n e)
newListVector n es = do
    x <- newVector_ n
    withVectorPtr x $ \p ->
        pokeArray p $ take n $ es ++ (repeat 0)
    return x
{-# INLINE newListVector #-}

-- | Set every element in the vector to zero.
setZeroVector :: (WriteVector x e m) => x n e -> m ()
setZeroVector = setZero
{-# INLINE setZeroVector #-}

-- | Set every element in the vector to a constant.
setConstantVector :: (WriteVector x e m) => e -> x n e -> m ()
setConstantVector = setConstant
{-# INLINE setConstantVector #-}

-- | Creats a new vector by copying another one.
newCopyVector :: (ReadVector x e m, WriteVector y e m) =>
    x n e -> m (y n e)
newCopyVector x
    | isConj x =
        newCopyVector (conj x) >>= return . conj
    | otherwise = do
        y <- newVector_ (dim x)
        unsafeCopyVector y x
        return y
{-# INLINE newCopyVector #-}

-- | Same as 'copyVector' but the sizes of the arguments are not checked.
unsafeCopyVector :: (WriteVector y e m, ReadVector x e m) =>
    y n e -> x n e -> m ()
unsafeCopyVector y x
    | isConj x && isConj y =
        unsafeCopyVector (conj y) (conj x)
    | isConj x || isConj y =
        forM_ [0..(dim x - 1)] $ \i -> do
            unsafeReadElem x i >>= unsafeWriteElem y i
    | otherwise =
        vectorCall2 BLAS.copy x y
{-# INLINE unsafeCopyVector #-}

-- | Same as 'swapVector' but the sizes of the arguments are not checked.
unsafeSwapVector :: (WriteVector y e m, BLAS1 e) =>
    y n e -> y n e -> m ()
unsafeSwapVector x y
    | isConj x && isConj y =
        unsafeSwapVector (conj x) (conj y)
    | isConj x || isConj y =
        forM_ [0..(dim x - 1)] $ \i -> do
            tmp <- unsafeReadElem x i
            unsafeReadElem y i >>= unsafeWriteElem x i
            unsafeWriteElem y i tmp
    | otherwise =
        vectorCall2 BLAS.swap x y
{-# INLINE unsafeSwapVector #-}

-- | Get a new vector with elements with the conjugates of the elements
-- of the given vector
getConjVector :: (ReadVector x e m, WriteVector y e m) =>
    x n e -> m (y n e)
getConjVector = getUnaryOp doConjVector
{-# INLINE getConjVector #-}

-- | Conjugate every element of the vector.
doConjVector :: (WriteVector y e m) => y n e -> m ()
doConjVector = doConj
{-# INLINE doConjVector #-}

-- | Get a new vector by scaling the elements of another vector
-- by a given value.
getScaledVector :: (ReadVector x e m, WriteVector y e m) =>
    e -> x n e -> m (y n e)
getScaledVector e = getUnaryOp (scaleByVector e)
{-# INLINE getScaledVector #-}

-- | Scale every element by the given value.
scaleByVector :: (WriteVector y e m) => e -> y n e -> m ()
scaleByVector = scaleBy
{-# INLINE scaleByVector #-}

-- | Get a new vector by shifting the elements of another vector
-- by a given value.
getShiftedVector :: (ReadVector x e m, WriteVector y e m) =>
    e -> x n e -> m (y n e)
getShiftedVector e = getUnaryOp (shiftByVector e)
{-# INLINE getShiftedVector #-}

-- | Add the given value to every element.
shiftByVector :: (WriteVector y e m) => e -> y n e -> m ()
shiftByVector = shiftBy
{-# INLINE shiftByVector #-}

-- | @getAddVector x y@ creates a new vector equal to the sum @x+y@.  The 
-- operands must have the same dimension.
getAddVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m) => 
    x n e -> y n e -> m (z n e)
getAddVector = checkTensorOp2 unsafeGetAddVector
{-# INLINE getAddVector #-}

unsafeGetAddVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m) => 
    x n e -> y n e -> m (z n e)
unsafeGetAddVector = unsafeGetBinaryOp unsafeAddVector
{-# INLINE unsafeGetAddVector #-}

-- | @addVector y x@ replaces @y@ with @y+x@.
addVector :: (WriteVector y e m, ReadVector x e m) => 
    y n e -> x n e -> m ()
addVector y x = checkBinaryOp (dim y) (dim x) $ unsafeAddVector y x
{-# INLINE addVector #-}

unsafeAddVector :: (WriteVector y e m, ReadVector x e m) => 
    y n e -> x n e -> m ()
unsafeAddVector y x = unsafeAxpyVector 1 x y
{-# INLINE unsafeAddVector #-}

-- | @getSubVector x y@ creates a new tensor equal to the difference @x-y@.  
-- The operands must have the same dimension.
getSubVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m) => 
    x n e -> y n e -> m (z n e)
getSubVector = checkTensorOp2 unsafeGetSubVector
{-# INLINE getSubVector #-}

unsafeGetSubVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m) => 
    x n e -> y n e -> m (z n e)
unsafeGetSubVector = unsafeGetBinaryOp unsafeSubVector
{-# INLINE unsafeGetSubVector #-}

-- | @subVector y x@ replaces @y@ with @y-x@.
subVector :: (WriteVector y e m, ReadVector x e m) => 
    y n e -> x n e -> m ()
subVector y x = checkBinaryOp (dim y) (dim x) $ unsafeSubVector y x
{-# INLINE subVector #-}

unsafeSubVector :: (WriteVector y e m, ReadVector x e m) => 
    y n e -> x n e -> m ()
unsafeSubVector y x = unsafeAxpyVector (-1) x y
{-# INLINE unsafeSubVector #-}

-- | @axpyVector alpha x y@ replaces @y@ with @alpha * x + y@.
axpyVector :: (ReadVector x e m, WriteVector y e m) =>
    e -> x n e -> y n e -> m ()
axpyVector alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyVector alpha x y
{-# INLINE axpyVector #-}

unsafeAxpyVector :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
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
    (ReadVector x e m, ReadVector y e m, WriteVector z e m) => 
    x n e -> y n e -> m (z n e)
getMulVector = checkTensorOp2 unsafeGetMulVector
{-# INLINE getMulVector #-}

unsafeGetMulVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m) => 
    x n e -> y n e -> m (z n e)
unsafeGetMulVector = unsafeGetBinaryOp unsafeMulVector
{-# INLINE unsafeGetMulVector #-}

-- | @mulVector y x@ replaces @y@ with @y*x@.
mulVector :: (WriteVector y e m, ReadVector x e m) => 
    y n e -> x n e -> m ()
mulVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeMulVector y x
{-# INLINE mulVector #-}
 
unsafeMulVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) =>
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
    (ReadVector x e m, ReadVector y e m, WriteVector z e m) => 
    x n e -> y n e -> m (z n e)
getDivVector = checkTensorOp2 unsafeGetDivVector
{-# INLINE getDivVector #-}

unsafeGetDivVector :: 
    (ReadVector x e m, ReadVector y e m, WriteVector z e m) => 
    x n e -> y n e -> m (z n e)
unsafeGetDivVector = unsafeGetBinaryOp unsafeDivVector
{-# INLINE unsafeGetDivVector #-}

-- | @divVector y x@ replaces @y@ with @y/x@.
divVector :: (WriteVector y e m, ReadVector x e m) => 
    y n e -> x n e -> m ()
divVector y x =
    checkBinaryOp (shape y) (shape x) $ unsafeDivVector y x
{-# INLINE divVector #-}

unsafeDivVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
unsafeDivVector y x
    | isConj y =
        unsafeDivVector (conj y) (conj x)
    | isConj x =
        vectorCall2 BLAS.vcdiv x y
    | otherwise =
        vectorCall2 BLAS.vdiv x y
{-# INLINE unsafeDivVector #-}


instance (Elem e) => BaseVector IOVector e where
    dim = dimIOVector
    {-# INLINE dim #-}
    stride = strideIOVector
    {-# INLINE stride #-}
    isConj = isConjIOVector
    {-# INLINE isConj #-}
    conj = conjIOVector
    {-# INLINE conj #-}
    withVectorPtrIO = withIOVectorPtr
    {-# INLINE withVectorPtrIO #-}
    
instance (BLAS1 e) => ReadVector IOVector e IO where
    withVectorPtr = withIOVectorPtr
    {-# INLINE withVectorPtr #-}
    unsafeFreezeVector = unsafeFreezeIOVector
    {-# INLINE unsafeFreezeVector #-}

instance (BLAS1 e) => WriteVector IOVector e IO where
    newVector_ = newIOVector_
    {-# INLINE newVector_ #-}
    unsafeThawVector = unsafeThawIOVector
    {-# INLINE unsafeThawVector #-}

-- | Create a vector of the given dimension with elements initialized
-- to the values from the list.  @listVector n es@ is equivalent to 
-- @vector n (zip [0..(n-1)] es)@, except that the result is undefined 
-- if @length es@ is less than @n@.
listVector :: (BLAS1 e) => Int -> [e] -> Vector n e
listVector n es = Vector $ unsafePerformIO $ newListVector n es
{-# NOINLINE listVector #-}
    
sizeVector :: Vector n e -> Int
sizeVector (Vector x) = sizeIOVector x
{-# INLINE sizeVector #-}

replaceVector :: (BLAS1 e) => Vector n e -> [(Int,e)] -> Vector n e
replaceVector (Vector x) ies =
    unsafePerformIO $ do
        y <- newCopyVector x
        mapM_ (uncurry $ writeElem y) ies
        return (Vector y)
{-# NOINLINE replaceVector #-}

unsafeReplaceVector :: (BLAS1 e) => Vector n e -> [(Int,e)] -> Vector n e
unsafeReplaceVector (Vector x) ies =
    unsafePerformIO $ do
        y <- newCopyVector x
        mapM_ (uncurry $ unsafeWriteElem y) ies
        return (Vector y)
{-# NOINLINE unsafeReplaceVector #-}

indicesVector :: Vector n e -> [Int]
indicesVector (Vector x) = indicesIOVector x
{-# INLINE indicesVector #-}

elemsVector :: (Elem e) => Vector n e -> [e]
elemsVector x | isConj x  = (map conjugate . elemsVector . conj) x
              | otherwise = case x of { (Vector (IOVector f p n incX _)) ->
    let end   = p `advancePtr` (n*incX)
        go p' | p' == end = inlinePerformIO $ touchForeignPtr f >> return []
              | otherwise = let e  = inlinePerformIO (peek p')
                                es = go (p' `advancePtr` incX)
                            in e:es
    in go p }
{-# SPECIALIZE INLINE elemsVector :: Vector n Double -> [Double] #-}
{-# SPECIALIZE INLINE elemsVector :: Vector n (Complex Double) -> [Complex Double] #-}

assocsVector :: (Elem e) => Vector n e -> [(Int,e)]
assocsVector x = zip (indicesVector x) (elemsVector x)
{-# INLINE assocsVector #-}

unsafeAtVector :: (Elem e) => Vector n e -> Int -> e
unsafeAtVector (Vector x) i = inlinePerformIO (unsafeReadElemIOVector x i)
{-# INLINE unsafeAtVector #-}

tmapVector :: (BLAS1 e) => (e -> e) -> Vector n e -> Vector n e
tmapVector f x = listVector (dim x) (map f $ elemsVector x)
{-# INLINE tmapVector #-}

tzipWithVector :: (BLAS1 e) =>
    (e -> e -> e) -> Vector n e -> Vector n e -> Vector n e
tzipWithVector f x y
    | dim y /= n =
        error ("tzipWith: vector lengths differ; first has length `" ++
                show n ++ "' and second has length `" ++
                show (dim y) ++ "'")
    | otherwise =
        listVector n (zipWith f (elems x) (elems y))
  where
    n = dim x
{-# INLINE tzipWithVector #-}

scaleVector :: (BLAS1 e) => e -> Vector n e -> Vector n e
scaleVector e (Vector x) = 
    unsafePerformIO $ unsafeFreezeIOVector =<< getScaledVector e x
{-# NOINLINE scaleVector #-}

shiftVector :: (BLAS1 e) => e -> Vector n e -> Vector n e
shiftVector e (Vector x) = 
    unsafePerformIO $ unsafeFreezeIOVector =<< getShiftedVector e x
{-# NOINLINE shiftVector #-}


instance Shaped Vector Int e where
    shape (Vector x) = shapeIOVector x
    {-# INLINE shape #-}
    bounds (Vector x) = boundsIOVector x
    {-# INLINE bounds #-}

instance (BLAS1 e) => ITensor Vector Int e where
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
    (*>) = scaleVector
    {-# INLINE (*>) #-}
    shift = shiftVector
    {-# INLINE shift #-}

instance (BLAS1 e, Monad m) => ReadTensor Vector Int e m where
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
    
instance (Elem e) => BaseVector Vector e where    
    dim (Vector x) = dimIOVector x
    {-# INLINE dim #-}
    stride (Vector x) = strideIOVector x
    {-# INLINE stride #-}
    isConj (Vector x) = isConjIOVector x
    {-# INLINE isConj #-}
    conj (Vector x) = (Vector (conjIOVector x))
    {-# INLINE conj #-}
    withVectorPtrIO (Vector x) = withIOVectorPtr x
    {-# INLINE withVectorPtrIO #-}

instance (BLAS1 e, Monad m) => ReadVector Vector e m where
    withVectorPtr (Vector x) f = (return . inlinePerformIO) $ withVectorPtr x f
    {-# INLINE withVectorPtr #-}
    unsafeFreezeVector = return
    {-# INLINE unsafeFreezeVector #-}

instance (Elem e, Show e) => Show (Vector n e) where
    show x
        | isConj x  = "conj (" ++ show (conj x) ++ ")"
        | otherwise = "listVector " ++ show (dim x) ++ " " ++ show (elemsVector x)

instance (BLAS1 e, Eq e) => Eq (Vector n e) where
    (==) = compareVectorWith (==)

instance (BLAS1 e, AEq e) => AEq (Vector n e) where
    (===) = compareVectorWith (===)
    (~==) = compareVectorWith (~==)

compareVectorWith :: (Elem e) => 
    (e -> e -> Bool) -> 
        Vector n e -> Vector n e -> Bool
compareVectorWith cmp x y
    | isConj x && isConj y =
        compareVectorWith cmp (conj x) (conj y)
    | otherwise =
        (dim x == dim y) && (and $ zipWith cmp (elemsVector x) (elemsVector y))

instance (BLAS1 e) => Num (Vector n e) where
    (+) x y = unsafePerformIO $ unsafeFreezeIOVector =<< getAddVector x y
    {-# NOINLINE (+) #-}
    (-) x y = unsafePerformIO $ unsafeFreezeIOVector =<< getSubVector x y
    {-# NOINLINE (-) #-}
    (*) x y = unsafePerformIO $ unsafeFreezeIOVector =<< getMulVector x y
    {-# NOINLINE (*) #-}
    negate = ((-1) *>)
    {-# INLINE negate #-}
    abs           = tmap abs
    signum        = tmap signum
    fromInteger n = listVector 1 [fromInteger n]
    
instance (BLAS1 e) => Fractional (Vector n e) where
    (/) x y      = unsafePerformIO $ unsafeFreezeIOVector =<< getDivVector x y
    {-# NOINLINE (/) #-}
    recip          = tmap recip
    fromRational q = listVector 1 [fromRational q]
    
instance (BLAS1 e, Floating e) => Floating (Vector n e) where
    pi    = listVector 1 [pi]
    exp   = tmap exp
    sqrt  = tmap sqrt
    log   = tmap log
    (**)  = tzipWithVector (**)
    sin   = tmap sin
    cos   = tmap cos
    tan   = tmap tan
    asin  = tmap asin
    acos  = tmap acos
    atan  = tmap atan
    sinh  = tmap sinh
    cosh  = tmap cosh
    tanh  = tmap tanh
    asinh = tmap asinh
    acosh = tmap acosh
    atanh = tmap atanh

vectorCall :: (ReadVector x e m) => 
    (Int -> Ptr e -> Int -> IO a) 
        ->  x n e -> m a
vectorCall f x =
    let n    = dim x
        incX = stride x
    in withVectorPtr x $ \pX -> 
           f n pX incX
{-# INLINE vectorCall #-}

vectorCall2 :: (ReadVector x e m, ReadVector y f m) => 
       (Int -> Ptr e -> Int -> Ptr f -> Int -> IO a) 
    -> x n e -> y n' f -> m a
vectorCall2 f x y =
    let n    = dim x
        incX = stride x
        incY = stride y
    in withVectorPtr x $ \pX ->
       withVectorPtrIO y $ \pY ->
           f n pX incX pY incY
{-# INLINE vectorCall2 #-}    

checkTensorOp2 :: (Shaped x i e, Shaped y i e) => 
    (x n e -> y n e -> a) ->
        x n e -> y n e -> a
checkTensorOp2 f x y = 
    checkBinaryOp (shape x) (shape y) $ f x y
{-# INLINE checkTensorOp2 #-}

getUnaryOp :: (ReadVector x e m, WriteVector y e m) =>
    (y n e -> m ()) -> x n e -> m (y n e)
getUnaryOp f x = do
    y <- newCopyVector x
    f y
    return y
{-# INLINE getUnaryOp #-}

unsafeGetBinaryOp :: 
    (WriteVector z e m, ReadVector x e m, ReadVector y e m) => 
    (z n e -> y n e -> m ()) ->
        x n e -> y n e -> m (z n e)
unsafeGetBinaryOp f x y = do
    z <- newCopyVector x
    f z y
    return z
{-# INLINE unsafeGetBinaryOp #-}
