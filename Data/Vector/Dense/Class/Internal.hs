{-# LANGUAGE BangPatterns, MultiParamTypeClasses, FunctionalDependencies,
        FlexibleContexts, FlexibleInstances #-}
{-# OPTIONS_HADDOCK hide, prune #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Internal (
    -- * Vector types
    IOVector,
    STVector,
    unsafeIOVectorToSTVector,
    unsafeSTVectorToIOVector,

    -- * Vector type classes
    BaseVector(..),
    ReadVector,
    WriteVector,

    -- * Basic vector propreties
    dim,
    stride,
    isConj,

    -- * Coercing the vector shape
    coerceVector,

    -- * BaseTensor functions
    shapeVector,
    boundsVector,

    -- * ReadTensor functions
    getSizeVector,
    getAssocsVector,
    getIndicesVector,
    getElemsVector,
    getAssocsVector',
    getIndicesVector',
    getElemsVector',
    unsafeReadElemVector,

    -- * WriteTensor functions
    newVector_,
    newZeroVector,
    setZeroVector,
    newConstantVector,
    setConstantVector,
    canModifyElemVector,
    unsafeWriteElemVector,
    modifyWithVector,

    -- * Numeric functions
    doConjVector,
    scaleByVector,
    shiftByVector,
    
    -- * CopyTensor functions
    newCopyVector,
    unsafeCopyVector,
    
    -- * SwapTensor functions
    unsafeSwapVector,
    
    -- * Numeric2 functions
    unsafeAxpyVector,
    unsafeMulVector,
    unsafeDivVector,
     
    -- * Utility functions
    withVectorPtr,
    indexOfVector,
    indicesVector,
    vectorCall,
    vectorCall2,

    ) where

import Control.Monad( forM_ )
import Control.Monad.ST
import Foreign
import Unsafe.Coerce

import BLAS.Internal ( clearArray, inlinePerformIO )
import BLAS.Elem
import qualified BLAS.C as BLAS
import BLAS.Tensor
import BLAS.UnsafeIOToM

import Data.Vector.Dense.Class.Internal.Base

class (BaseVector x e, BLAS1 e, UnsafeIOToM m, ReadTensor x Int e m) => 
    ReadVector x e m
class (ReadVector x e m, WriteTensor x Int e m) => 
    WriteVector x e m | x -> m, m -> x


-- | Cast the shape type of the vector.
coerceVector :: (BaseVector x e) => x n e -> x n' e
coerceVector = unsafeCoerce
{-# INLINE coerceVector #-}

-------------------------- BaseTensor functions -----------------------------

shapeVector :: (BaseVector x e) => x n e -> Int
shapeVector = dim
{-# INLINE shapeVector #-}

boundsVector :: (BaseVector x e) => x n e -> (Int,Int)
boundsVector x = (0, dim x - 1)
{-# INLINE boundsVector #-}

-------------------------- ReadTensor functions -----------------------------
    
getSizeVector :: (ReadVector x e m) => x n e -> m Int
getSizeVector = return . dim
{-# INLINE getSizeVector #-}

getIndicesVector :: (ReadVector x e m) => x n e -> m [Int]
getIndicesVector = return . indicesVector
{-# INLINE getIndicesVector #-}

getIndicesVector' :: (ReadVector x e m) => x n e -> m [Int]
getIndicesVector' = getIndicesVector
{-# INLINE getIndicesVector' #-}

getElemsVector :: (ReadVector x e m) => x n e -> m [e]
getElemsVector x = do
    ies <- getAssocsVector x
    return $ (snd . unzip) ies

getElemsVector' :: (ReadVector x e m) => x n e -> m [e]
getElemsVector' x = do
    ies <- getAssocsVector' x
    return $ (snd . unzip) ies

getAssocsVector :: (ReadVector x e m) => x n e -> m [(Int,e)]
getAssocsVector x
    | isConj x =
        getAssocsVector (conj x) 
            >>= return . map (\(i,e) -> (i,conj e))
    | otherwise =
        let (f,p,n,incX,_) = arrayFromVector x
        in return $ go n f incX p 0
  where
        go !n !f !incX !ptr !i 
            | i >= n = 
                 -- This is very important since we are doing unsafe IO.
                 -- Otherwise, the DVector might get discared and the
                 -- memory freed before all of the elements are read
                 inlinePerformIO $ do
                     touchForeignPtr f
                     return []
            | otherwise =
                let e    = inlinePerformIO $ peek ptr
                    ptr' = ptr `advancePtr` incX
                    i'   = i + 1
                    ies  = go n f incX ptr' i'
                in e `seq` ((i,e):ies)
{-# INLINE getAssocsVector #-}

getAssocsVector' :: (ReadVector x e m) => x n e -> m [(Int,e)]
getAssocsVector' x
    | isConj x =
        getAssocsVector' (conj x) 
            >>= return . map (\(i,e) -> (i,conj e))
    | otherwise =
        unsafeIOToM $
            withVectorPtr x $ \ptr ->
                go (dim x) (stride x) ptr 0
  where
        go !n !incX !ptr !i 
            | i >= n = 
                return []
            | otherwise = do
                e <- peek ptr
                let ptr' = ptr `advancePtr` incX
                    i'   = i + 1
                ies <- go n incX ptr' i'
                return $ (i,e):ies

unsafeReadElemVector :: (ReadVector x e m) => x n e -> Int -> m e
unsafeReadElemVector x i
    | isConj x = 
        unsafeReadElemVector (conj x) i >>= return . conj
    | otherwise =
        unsafeIOToM $
            withVectorPtr x $ \ptr ->
                peekElemOff ptr (indexOfVector x i) 


------------------------- WriteTensor functions -----------------------------

-- | Creates a new vector of the given length.  The elements will be 
-- uninitialized.
newVector_ :: (WriteVector x e m) => Int -> m (x n e)
newVector_ n
    | n < 0 = 
        fail $ "Tried to create a vector with `" ++ show n ++ "' elements."
    | otherwise = unsafeIOToM $ do
        arr <- mallocForeignPtrArray n
        return $ vectorViewArray arr (unsafeForeignPtrToPtr arr) n 1 False

-- | Create a zero vector of the specified length.
newZeroVector :: (WriteVector y e m) => Int -> m (y n e)
newZeroVector n = do
    x <- newVector_ n
    setZeroVector x
    return x

-- | Set every element in the vector to zero.
setZeroVector :: (WriteVector y e m) => y n e -> m ()
setZeroVector x 
    | stride x == 1 = unsafeIOToM $
                          withVectorPtr x $ 
                              flip clearArray (dim x)
    | otherwise     = setConstantVector 0 x

-- | Create a constant vector of the specified length.
newConstantVector :: (WriteVector y e m) => Int -> e -> m (y n e)
newConstantVector n e = do
    x <- newVector_ n
    setConstantVector e x
    return x
        
-- | Set every element in the vector to a constant.
setConstantVector :: (WriteVector y e m) => e -> y n e -> m ()
setConstantVector e x 
    | isConj x  = setConstantVector (conj e) (conj x)
    | otherwise = unsafeIOToM $ withVectorPtr x $ go (dim x)
  where
    go !n !ptr | n <= 0 = return ()
               | otherwise = let ptr' = ptr `advancePtr` (stride x)
                                 n'   = n - 1
                             in poke ptr e >> 
                                go n' ptr'

canModifyElemVector :: (WriteVector y e m) => y n e -> Int -> m Bool
canModifyElemVector _ _ = return True
{-# INLINE canModifyElemVector #-}

unsafeWriteElemVector :: (WriteVector y e m) => y n e -> Int -> e -> m ()
unsafeWriteElemVector x i e =
    let e' = if isConj x then conj e else e
    in unsafeIOToM $ withVectorPtr x $ \ptr -> 
           pokeElemOff ptr (indexOfVector x i) e'
                    
modifyWithVector :: (WriteVector y e m) => (e -> e) -> y n e -> m ()
modifyWithVector f x
    | isConj x = modifyWithVector (conj . f . conj) (conj x)
    | otherwise = unsafeIOToM $
                      withVectorPtr x $ \ptr ->
                          go (dim x) ptr
  where
    go !n !ptr | n <= 0 = return ()
               | otherwise = do
                   peek ptr >>= poke ptr . f
                   go (n-1) (ptr `advancePtr` incX)

    incX = stride x


------------------------- CopyTensor functions ------------------------------

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


------------------------- SwapTensor functions ------------------------------

-- | Same as 'swapVector' but the sizes of the arguments are not checked.
unsafeSwapVector :: (WriteVector y e m) => 
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


--------------------------- Numeric functions -------------------------------

doConjVector :: (WriteVector y e m) => y n e -> m ()
doConjVector = vectorCall BLAS.conj

scaleByVector :: (WriteVector y e m) => e -> y n e -> m ()
scaleByVector 1 _ = return ()
scaleByVector k x | isConj x    = scaleByVector (conj k) (conj x)
                    | otherwise = vectorCall (flip BLAS.scal k) x
                    
shiftByVector :: (WriteVector y e m) => e -> y n e -> m ()
shiftByVector k x | isConj x  = shiftByVector (conj k) (conj x)
                  | otherwise = modifyWithVector (k+) x


-------------------------- Numeric2 functions -------------------------------

unsafeAxpyVector :: (ReadVector x e m, WriteVector y e m) => 
    e -> x n e -> y n e -> m ()
unsafeAxpyVector alpha x y
    | isConj y =
        unsafeAxpyVector (conj alpha) (conj x) (conj y)
    | isConj x =
        vectorCall2 (flip BLAS.acxpy alpha) x y
    | otherwise =
        vectorCall2 (flip BLAS.axpy alpha) x y

unsafeMulVector :: (WriteVector y e m, ReadVector x e m) => 
    y n e -> x n e -> m ()
unsafeMulVector y x
    | isConj y =
        unsafeMulVector (conj y) (conj x)
    | isConj x =
        vectorCall2 BLAS.cmul x y
    | otherwise =
        vectorCall2 BLAS.mul x y

unsafeDivVector :: (WriteVector y e m, ReadVector x e m) => 
    y n e -> x n e -> m ()
unsafeDivVector y x
    | isConj y =
        unsafeDivVector (conj y) (conj x)
    | isConj x =
        vectorCall2 BLAS.cdiv x y
    | otherwise =
        vectorCall2 BLAS.div x y


--------------------------- Utility functions -------------------------------

indexOfVector :: (BaseVector x e) => x n e -> Int -> Int
indexOfVector x i = i * stride x
{-# INLINE indexOfVector #-}

indicesVector :: (BaseVector x e) => x n e -> [Int]
indicesVector x = [0..(n-1)] where n = dim x
{-# INLINE indicesVector #-}

vectorCall :: (ReadVector x e m) => 
    (Int -> Ptr e -> Int -> IO a) 
        ->  x n e -> m a
vectorCall f x =
    let n    = dim x
        incX = stride x
    in unsafeIOToM $
           withVectorPtr x $ \pX -> 
               f n pX incX
{-# INLINE vectorCall #-}

vectorCall2 :: (ReadVector x e m, ReadVector y f m) => 
       (Int -> Ptr e -> Int -> Ptr f -> Int -> IO a) 
    -> x n e -> y n' f -> m a
vectorCall2 f x y =
    let n    = dim x
        incX = stride x
        incY = stride y
    in unsafeIOToM $
           withVectorPtr x $ \pX ->
               withVectorPtr y $ \pY ->
                   f n pX incX pY incY
{-# INLINE vectorCall2 #-}    


------------------------------------ Instances ------------------------------

data IOVector n e = 
      DV {-# UNPACK #-} !(ForeignPtr e) -- memory owner
         {-# UNPACK #-} !(Ptr e)        -- pointer to the first element
         {-# UNPACK #-} !Int            -- the length of the vector
         {-# UNPACK #-} !Int            -- the stride (in elements, not bytes) between elements.
         {-# UNPACK #-} !Bool           -- indicates whether or not the vector is conjugated

newtype STVector s n e = ST (IOVector n e)

unsafeIOVectorToSTVector :: IOVector n e -> STVector s n e
unsafeIOVectorToSTVector = ST

unsafeSTVectorToIOVector :: STVector s n e -> IOVector n e
unsafeSTVectorToIOVector (ST x) = x


instance (Storable e) => BaseVector IOVector e where
    vectorViewArray                = DV
    {-# INLINE vectorViewArray #-}
    
    arrayFromVector (DV f p n s c) = (f,p,n,s,c)
    {-# INLINE arrayFromVector #-}

instance (Storable e) => BaseVector (STVector s) e where
    vectorViewArray f p n s c = ST $ DV f p n s c
    {-# INLINE vectorViewArray #-}    
    
    arrayFromVector (ST x)    = arrayFromVector x
    {-# INLINE arrayFromVector #-}

instance (Storable e) => BaseTensor IOVector Int e where
    bounds = boundsVector
    shape  = shapeVector

instance (Storable e) => BaseTensor (STVector s) Int e where
    bounds = boundsVector
    shape  = shapeVector
        
instance (BLAS1 e) => ReadTensor IOVector Int e IO where
    getSize        = getSizeVector
    getAssocs      = getAssocsVector
    getIndices     = getIndicesVector
    getElems       = getElemsVector
    getAssocs'     = getAssocsVector'
    getIndices'    = getIndicesVector'
    getElems'      = getElemsVector'
    unsafeReadElem = unsafeReadElemVector

instance (BLAS1 e) => ReadTensor (STVector s) Int e (ST s) where
    getSize        = getSizeVector
    getAssocs      = getAssocsVector
    getIndices     = getIndicesVector
    getElems       = getElemsVector
    getAssocs'     = getAssocsVector'
    getIndices'    = getIndicesVector'
    getElems'      = getElemsVector'
    unsafeReadElem = unsafeReadElemVector

instance (BLAS1 e) => ReadVector IOVector     e IO where
instance (BLAS1 e) => ReadVector (STVector s) e (ST s) where    

instance (BLAS1 e) => WriteTensor IOVector Int e IO where
    setConstant     = setConstantVector
    setZero         = setZeroVector
    canModifyElem   = canModifyElemVector
    unsafeWriteElem = unsafeWriteElemVector
    modifyWith      = modifyWithVector
    doConj          = doConjVector
    scaleBy         = scaleByVector
    shiftBy         = shiftByVector
    
instance (BLAS1 e) => WriteTensor (STVector s) Int e (ST s) where
    setConstant     = setConstantVector
    setZero         = setZeroVector
    canModifyElem   = canModifyElemVector
    unsafeWriteElem = unsafeWriteElemVector
    modifyWith      = modifyWithVector
    doConj          = doConjVector
    scaleBy         = scaleByVector
    shiftBy         = shiftByVector

instance (BLAS1 e) => WriteVector IOVector     e IO where
instance (BLAS1 e) => WriteVector (STVector s) e (ST s) where
