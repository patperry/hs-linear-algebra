{-# LANGUAGE BangPatterns, FlexibleInstances, MultiParamTypeClasses #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.ST
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Internal (
    -- * The IOVector data type
    IOVector,

    -- * Vector Properties
    getSumAbs,
    getNorm2,
    getWhichMaxAbs,
    
    -- * Vector Operations
    getDot,
    unsafeGetDot,
    
    -- * Internal Operations
    unsafeAxpyVector,
    unsafeMulVector,
    unsafeDivVector,
    
    module BLAS.Tensor,
    module BLAS.Numeric,
    module Data.Vector.Dense.Class,

    ) where
  
import BLAS.Elem.Base ( Elem )

import BLAS.Internal ( clearArray, inlinePerformIO, checkVecVecOp )

import BLAS.C ( BLAS1 )
import qualified BLAS.C as BLAS

import BLAS.Tensor hiding ( ITensor(..) )
import BLAS.Numeric hiding ( INumeric(..) )

import Control.Monad
import Data.Complex

import Data.Vector.Dense.Class hiding ( conjVector )
import Data.Vector.Dense.Class.Base

import Foreign 


-- | A dense mutable vector.  @s@ is the state thread,
-- @n@ is a phantom type for the dimension of the vector, and @e@ is the 
-- element type.  An @IOVector@ @x@ stores @dimIOVector x@ elements.  Indices into
-- the vector are @0@-based.
data IOVector n e = 
      DV { storageIOVector :: {-# UNPACK #-} !(ForeignPtr e) -- ^ a pointer to the storage region
         , offsetIOVector  :: {-# UNPACK #-} !Int            -- ^ an offset (in elements, not bytes) to the first element in the vector. 
         , dimIOVector     :: {-# UNPACK #-} !Int            -- ^ the length of the vector
         , strideIOVector  :: {-# UNPACK #-} !Int            -- ^ the stride (in elements, not bytes) between elements.
         , isConjIOVector  :: {-# UNPACK #-} !Bool           -- ^ indicates whether or not the vector is conjugated
         }

unsafeSubvectorWithStrideIOVector :: Int -> IOVector n e -> Int -> Int -> IOVector m e
unsafeSubvectorWithStrideIOVector s x o n =
    let f  = storageIOVector x
        o' = indexOfIOVector x o
        n' = n
        s' = s * (strideIOVector x)
        c  = isConjIOVector x
    in DV f o' n' s' c

newIOVector_ :: (Elem e) => Int -> IO (IOVector n e)
newIOVector_ n
    | n < 0 = 
        fail $ "Tried to create a vector with `" ++ show n ++ "' elements."
    | otherwise = do
        arr <- mallocForeignPtrArray n
        return $ DV arr 0 n 1 False

indexOfIOVector :: IOVector n e -> Int -> Int
indexOfIOVector x i = offsetIOVector x + i * strideIOVector x
{-# INLINE indexOfIOVector #-}

-- | Evaluate a function with a pointer to the value stored at the given
-- index.  Note that the value may need to conjugated before using it.  See
-- 'isConj'.
unsafeWithElemPtrIOVector :: (Storable e) => IOVector n e -> Int -> (Ptr e -> IO a) -> IO a
unsafeWithElemPtrIOVector x i f =
    withForeignPtr (storageIOVector x) $ \ptr ->
        let elemPtr = ptr `advancePtr` (indexOfIOVector x i)
        in f elemPtr
{-# INLINE unsafeWithElemPtrIOVector #-}

withIOVectorPtr :: (Storable e) => IOVector n e -> (Ptr e -> IO a) -> IO a
withIOVectorPtr = flip unsafeWithElemPtrIOVector 0

conjIOVector :: IOVector n e -> IOVector n e
conjIOVector x = let c' = (not . isConjIOVector) x 
                 in x { isConjIOVector=c' }
{-# INLINE conjIOVector #-}

getSizeIOVector :: IOVector n e -> IO Int
getSizeIOVector = return . dimIOVector
{-# INLINE getSizeIOVector #-}

shapeIOVector :: IOVector n e -> Int
shapeIOVector = dimIOVector

boundsIOVector :: IOVector n e -> (Int,Int)
boundsIOVector x = (0, dimIOVector x - 1)

indicesIOVector :: IOVector n e -> [Int]
indicesIOVector x = [0..(n-1)] where n = dimIOVector x
{-# INLINE indicesIOVector #-}

getIndicesIOVector :: IOVector n e -> IO [Int]
getIndicesIOVector = return . indicesIOVector
{-# INLINE getIndicesIOVector #-}

getIndicesIOVector' :: IOVector n e -> IO [Int]
getIndicesIOVector' = getIndicesIOVector
{-# INLINE getIndicesIOVector' #-}

unsafeReadElemIOVector :: (Elem e) => IOVector n e -> Int -> IO e
unsafeReadElemIOVector x i
    | isConjIOVector x = 
        unsafeReadElemIOVector (conjIOVector x) i >>= return . conj
    | otherwise =
        withForeignPtr (storageIOVector x) $ \ptr ->
            peekElemOff ptr (indexOfIOVector x i) 

getElemsIOVector :: (Elem e) => IOVector n e -> IO [e]
getElemsIOVector x = do
    ies <- getAssocsIOVector x
    return $ (snd . unzip) ies

getElemsIOVector' :: (Elem e) => IOVector n e -> IO [e]
getElemsIOVector' x = do
    ies <- getAssocsIOVector' x
    return $ (snd . unzip) ies

getAssocsIOVector :: (Elem e) => IOVector n e -> IO [(Int,e)]
getAssocsIOVector x
    | isConjIOVector x =
        getAssocsIOVector (conjIOVector x) 
            >>= return . map (\(i,e) -> (i,conj e))
    | otherwise =
        let (DV f o n incX _) = x
            ptr = (unsafeForeignPtrToPtr f) `advancePtr` o
        in return $ go n f incX ptr 0
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
{-# INLINE getAssocsIOVector #-}

getAssocsIOVector' :: (Elem e) => IOVector n e -> IO [(Int,e)]
getAssocsIOVector' x
    | isConjIOVector x =
        getAssocsIOVector' (conjIOVector x) 
            >>= return . map (\(i,e) -> (i,conj e))
    | otherwise =
        withIOVectorPtr x $ \ptr ->
            go (dimIOVector x) (strideIOVector x) ptr 0
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


newZeroIOVector :: (Elem e) => Int -> IO (IOVector n e)
newZeroIOVector n = do
    x <- newIOVector_ n
    setZeroIOVector x
    return x

newConstantIOVector :: (Elem e) => Int -> e -> IO (IOVector n e)
newConstantIOVector n e = do
    x <- newIOVector_ n
    setConstantIOVector e x
    return x
    
setZeroIOVector :: (Elem e) => IOVector n e -> IO ()    
setZeroIOVector x 
    | strideIOVector x == 1 = unsafeWithElemPtrIOVector x 0 $ 
                                    flip clearArray (dimIOVector x)
    | otherwise               = setConstantIOVector 0 x

setConstantIOVector :: (Elem e) => e -> IOVector n e -> IO ()    
setConstantIOVector e x 
    | isConjIOVector x  = setConstantIOVector (conj e) (conjIOVector x)
    | otherwise = unsafeWithElemPtrIOVector x 0 $ go (dimIOVector x)
  where
    go !n !ptr | n <= 0 = return ()
               | otherwise = let ptr' = ptr `advancePtr` (strideIOVector x)
                                 n'   = n - 1
                             in poke ptr e >> 
                                go n' ptr'

unsafeWriteElemIOVector :: (Elem e) => IOVector n e -> Int -> e -> IO ()
unsafeWriteElemIOVector x i e =
    let e' = if isConjIOVector x then conj e else e
    in withForeignPtr (storageIOVector x) $ \ptr -> 
           pokeElemOff ptr (indexOfIOVector x i) e'
                    
canModifyElemIOVector :: IOVector n e -> Int -> IO Bool
canModifyElemIOVector _ _ = return True
{-# INLINE canModifyElemIOVector #-}

modifyWithIOVector :: (Elem e) => (e -> e) -> IOVector n e -> IO ()
modifyWithIOVector f x
    | isConjIOVector x = modifyWithIOVector (conj . f . conj) (conjIOVector x)
    | otherwise = withIOVectorPtr x $ \ptr ->
                      go (dimIOVector x) ptr
  where
    go !n !ptr | n <= 0 = return ()
               | otherwise = do
                   peek ptr >>= poke ptr . f
                   go (n-1) (ptr `advancePtr` incX)

    incX = strideIOVector x

---------------------------  Swapping Vectors --------------------------------

unsafeSwapIOVector :: (Elem e) => IOVector n e -> IOVector n e -> IO ()
unsafeSwapIOVector x y
    | isConjIOVector x && isConjIOVector y =
        unsafeSwapIOVector (conjIOVector x) (conjIOVector y)
    | otherwise =
        forM_ [0..(dimIOVector x - 1)] $ \i -> do
            tmp <- unsafeReadElemIOVector x i
            unsafeReadElemIOVector y i >>= unsafeWriteElemIOVector x i
            unsafeWriteElemIOVector y i tmp

unsafeSwapIOVectorBLAS1 :: (BLAS1 e) => 
    IOVector n e -> IOVector n e -> IO ()
unsafeSwapIOVectorBLAS1 x y
    | isConjIOVector x && isConjIOVector y =
        unsafeSwapIOVectorBLAS1 (conj x) (conj y)
    | isConjIOVector x || isConjIOVector y =
        forM_ [0..(dimIOVector x - 1)] $ \i -> do
            tmp <- unsafeReadElemIOVector x i
            unsafeReadElemIOVector y i >>= unsafeWriteElemIOVector x i
            unsafeWriteElemIOVector y i tmp
    | otherwise =
        call2 BLAS.swap x y

unsafeSwapIOVectorBLAS1Real :: (BLAS1 e) => 
    IOVector n e -> IOVector n e -> IO ()
unsafeSwapIOVectorBLAS1Real = call2 BLAS.swap

unsafeSwapIOVectorComplexDouble :: IOVector n (Complex Double) -> IOVector n (Complex Double) -> IO ()
unsafeSwapIOVectorComplexDouble = unsafeSwapIOVectorBLAS1
unsafeSwapIOVectorDouble :: IOVector n Double -> IOVector n Double -> IO ()
unsafeSwapIOVectorDouble = unsafeSwapIOVectorBLAS1Real
{-# RULES "unsafeSwapIOVector/Double" unsafeSwapIOVector = unsafeSwapIOVectorDouble #-}
{-# RULES "unsafeSwapIOVector/ComplexDouble" unsafeSwapIOVector = unsafeSwapIOVectorComplexDouble #-}


----------------------------- BLAS Operations --------------------------------

-- | Gets the sum of the absolute values of the vector entries.
getSumAbs :: (ReadVector x e m, BLAS1 e) => x n e -> m Double
getSumAbs = call BLAS.asum
    
-- | Gets the 2-norm of a vector.
getNorm2 :: (ReadVector x e m, BLAS1 e) => x n e -> m Double
getNorm2 = call BLAS.nrm2

-- | Gets the index and norm of the element with maximum magnitude.  This is 
-- undefined if any of the elements are @NaN@.  It will throw an exception if 
-- the dimension of the vector is 0.
getWhichMaxAbs :: (ReadVector x e m, BLAS1 e) => x n e -> m (Int, e)
getWhichMaxAbs x =
    case (dim x) of
        0 -> fail $ "getWhichMaxAbs of an empty vector"
        _ -> do
            i <- call BLAS.iamax x
            e <- unsafeReadElem x i
            return (i,e)


-- | Computes the dot product of two vectors.
getDot :: (ReadVector x e m, ReadVector y e m, BLAS1 e) => 
    x n e -> y n e -> m e
getDot x y = checkVecVecOp "getDot" (dim x) (dim y) $ unsafeGetDot x y
{-# INLINE getDot #-}

unsafeGetDot :: (ReadVector x e m, ReadVector y e m, BLAS1 e) => 
    x n e -> y n e -> m e
unsafeGetDot x y =
    case (isConj x, isConj y) of
        (False, False) -> call2 BLAS.dotc x y
        (True , False) -> call2 BLAS.dotu x y
        (False, True ) -> call2 BLAS.dotu x y >>= return . conj
        (True , True)  -> call2 BLAS.dotc x y >>= return . conj
{-# INLINE unsafeGetDot #-}

doConjIOVector :: (BLAS1 e) => IOVector n e -> IO ()
doConjIOVector = call BLAS.conj

scaleByIOVector :: (BLAS1 e) => e -> IOVector n e -> IO ()
scaleByIOVector 1 _ = return ()
scaleByIOVector k x | isConjIOVector x= scaleByIOVector (conj k) (conj x)
                    | otherwise       = call (flip BLAS.scal k) x
                    
shiftByIOVector :: (Elem e) => e -> IOVector n e -> IO ()
shiftByIOVector k x | isConjIOVector x  = shiftByIOVector (conj k) (conj x)
                    | otherwise = modifyWithIOVector (k+) x

unsafeAxpyVector :: (ReadVector x e m, WriteVector y e m, BLAS1 e) => 
    e -> x n e -> y n e -> m ()
unsafeAxpyVector alpha x y
    | isConj y =
        unsafeAxpyVector (conj alpha) (conj x) (conj y)
    | isConj x =
        call2 (flip BLAS.acxpy alpha) x y
    | otherwise =
        call2 (flip BLAS.axpy alpha) x y

unsafeMulVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
unsafeMulVector y x
    | isConj y =
        unsafeMulVector (conj y) (conj x)
    | isConj x =
        call2 BLAS.cmul x y
        --call2 (flip (BLAS.tbmv BLAS.colMajor BLAS.upper BLAS.conjTrans BLAS.nonUnit) 0) x y    
    | otherwise =
        call2 BLAS.mul x y
        --call2 (flip (BLAS.tbmv BLAS.colMajor BLAS.upper BLAS.noTrans BLAS.nonUnit) 0) x y

unsafeDivVector :: (WriteVector y e m, ReadVector x e m, BLAS1 e) => 
    y n e -> x n e -> m ()
unsafeDivVector y x
    | isConj y =
        unsafeDivVector (conj y) (conj x)
    | isConj x =
        call2 BLAS.cdiv x y
        --call2 (flip (BLAS.tbsv BLAS.colMajor BLAS.upper BLAS.conjTrans BLAS.nonUnit) 0) x y
    | otherwise =
        call2 BLAS.div x y
        --call2 (flip (BLAS.tbsv BLAS.colMajor BLAS.upper BLAS.noTrans BLAS.nonUnit) 0) x y


call :: (ReadVector x e m) => 
        (Int -> Ptr e -> Int -> IO a) 
    ->  x n e -> m a
call f x =
    let n    = dim x
        incX = stride x
    in unsafeIOToM $
           withVectorPtr x $ \pX -> 
               f n pX incX
{-# INLINE call #-}

call2 :: (ReadVector x e m, ReadVector y f m) => 
       (Int -> Ptr e -> Int -> Ptr f -> Int -> IO a) 
    -> x n e -> y n' f -> m a
call2 f x y =
    let n    = dim x
        incX = stride x
        incY = stride y
    in unsafeIOToM $
           withVectorPtr x $ \pX ->
               withVectorPtr y $ \pY ->
                   f n pX incX pY incY
{-# INLINE call2 #-}    

instance BaseTensor IOVector Int e where
    shape  = shapeIOVector
    bounds = boundsIOVector

instance (Elem e) => ReadTensor IOVector Int e IO where
    getSize        = getSizeIOVector

    getAssocs      = getAssocsIOVector
    getIndices     = getIndicesIOVector
    getElems       = getElemsIOVector

    getAssocs'     = getAssocsIOVector'
    getIndices'    = getIndicesIOVector'
    getElems'      = getElemsIOVector'
    
    unsafeReadElem = unsafeReadElemIOVector

instance (Elem e) => WriteTensor IOVector Int e IO where
    newZero         = newZeroIOVector
    newConstant     = newConstantIOVector

    setConstant     = setConstantIOVector
    setZero         = setZeroIOVector

    modifyWith      = modifyWithIOVector
    unsafeWriteElem = unsafeWriteElemIOVector
    canModifyElem   = canModifyElemIOVector
    
    unsafeSwap      = unsafeSwapIOVector
    
instance (Elem e) => BaseVector IOVector e where
    stride                    = strideIOVector
    isConj                    = isConjIOVector
    conjVector                = conjIOVector
    unsafeSubvectorWithStride = unsafeSubvectorWithStrideIOVector
    withVectorPtr             = withIOVectorPtr
    vectorViewArray           = DV
    arrayFromVector (DV f o n s c) = (f,o,n,s,c)

instance (Elem e) => ReadVector IOVector e IO where
    
instance (Elem e) => WriteVector IOVector e IO where
    newVector_ = newIOVector_

instance (BLAS1 e) => ReadNumeric IOVector Int e IO where

instance (BLAS1 e) => WriteNumeric IOVector Int e IO where
    doConj  = doConjIOVector
    scaleBy = scaleByIOVector
    shiftBy = shiftByIOVector

instance (BLAS1 e, ReadVector x e IO) => CopyTensor x IOVector Int e IO where
    newCopyTensor    = newCopyVector    
    unsafeCopyTensor = unsafeCopyVector

instance (BLAS1 e, ReadVector x e IO) => Numeric2 x IOVector Int e IO where
    unsafeAxpy = unsafeAxpyVector
    unsafeMul  = unsafeMulVector
    unsafeDiv  = unsafeDivVector

instance (BLAS1 e, ReadVector x e IO, ReadVector y e IO) => Numeric3 x y IOVector Int e IO where
    unsafeDoAdd = unsafeDoVectorOp2 $ flip $ unsafeAxpyVector 1
    unsafeDoSub = unsafeDoVectorOp2 $ flip $ unsafeAxpyVector (-1)
    unsafeDoMul = unsafeDoVectorOp2 $ unsafeMulVector
    unsafeDoDiv = unsafeDoVectorOp2 $ unsafeDivVector
        