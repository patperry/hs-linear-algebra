{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.IOBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.IOBase
    where

import Control.Monad
import Foreign
import System.IO.Unsafe

import BLAS.Internal ( clearArray )
import Data.Elem.BLAS ( Complex, Elem, BLAS1, conjugate )
import qualified Data.Elem.BLAS as BLAS

import Data.Tensor.Class
import Data.Tensor.Class.MTensor


data IOVector n e = 
      IOVector {-# UNPACK #-} !(ForeignPtr e) -- memory owner
               {-# UNPACK #-} !(Ptr e)        -- pointer to the first element
               {-# UNPACK #-} !Int            -- the length of the vector
               {-# UNPACK #-} !Int            -- the stride (in elements, not bytes) between elements.
               {-# UNPACK #-} !Bool           -- indicates whether or not the vector is conjugated

dimIOVector :: IOVector n e -> Int
dimIOVector (IOVector _ _ n _ _) = n
{-# INLINE dimIOVector #-}

strideIOVector :: IOVector n e -> Int
strideIOVector (IOVector _ _ _ incX _) = incX
{-# INLINE strideIOVector #-}

isConjIOVector :: IOVector n e -> Bool
isConjIOVector (IOVector _ _ _ _ c) = c
{-# INLINE isConjIOVector #-}

conjIOVector :: IOVector n e -> IOVector n e
conjIOVector (IOVector f p n incX c) = (IOVector f p n incX (not c))
{-# INLINE conjIOVector #-}

unsafeSubvectorViewWithStrideIOVector :: (Elem e) =>
    Int -> IOVector n e -> Int -> Int -> IOVector n' e
unsafeSubvectorViewWithStrideIOVector s' (IOVector f p _ inc c) o' n' =
    IOVector f (p `advancePtr` (inc*o')) n' (inc*s') c
{-# INLINE unsafeSubvectorViewWithStrideIOVector #-}

withIOVectorPtr :: IOVector n e -> (Ptr e -> IO a) -> IO a
withIOVectorPtr (IOVector f p _ _ _) g = do
    a <- g p
    touchForeignPtr f
    return a
{-# INLINE withIOVectorPtr #-}

newIOVector_ :: (Elem e) => Int -> IO (IOVector n e)
newIOVector_ n
    | n < 0 = 
        fail $ "Tried to create a vector with `" ++ show n ++ "' elements."
    | otherwise = do
        arr <- mallocForeignPtrArray n
        return $ IOVector arr (unsafeForeignPtrToPtr arr) n 1 False

newCopyIOVector :: (BLAS1 e) => IOVector n e -> IO (IOVector n e)
newCopyIOVector (IOVector f p n incX c) = do
    (IOVector f' p' _ _ _) <- newIOVector_ n
    BLAS.copy n p incX p' 1
    touchForeignPtr f
    touchForeignPtr f'
    return (IOVector f' p' n 1 c)

shapeIOVector :: IOVector n e -> Int
shapeIOVector = dimIOVector
{-# INLINE shapeIOVector #-}

boundsIOVector :: IOVector n e -> (Int,Int)
boundsIOVector x = (0, dimIOVector x - 1)
{-# INLINE boundsIOVector #-}

sizeIOVector :: IOVector n e -> Int
sizeIOVector = dimIOVector
{-# INLINE sizeIOVector #-}

getSizeIOVector :: IOVector n e -> IO Int
getSizeIOVector = return . sizeIOVector
{-# INLINE getSizeIOVector #-}

getMaxSizeIOVector :: IOVector n e -> IO Int
getMaxSizeIOVector = getSizeIOVector
{-# INLINE getMaxSizeIOVector #-}

indicesIOVector :: IOVector n e -> [Int]
indicesIOVector x = [ 0..n-1 ] where n = dimIOVector x
{-# INLINE indicesIOVector #-}

getIndicesIOVector :: IOVector n e -> IO [Int]
getIndicesIOVector = return . indicesIOVector
{-# INLINE getIndicesIOVector #-}

getIndicesIOVector' :: IOVector n e -> IO [Int]
getIndicesIOVector' = getIndicesIOVector
{-# INLINE getIndicesIOVector' #-}

getElemsIOVector :: (Elem e) => IOVector n e -> IO [e]
getElemsIOVector (IOVector f p n incX True) = do
    es <- getElemsIOVector (IOVector f p n incX False)
    return $ map conjugate es
getElemsIOVector (IOVector f p n incX False) =
    let end = p `advancePtr` (n*incX)
        go p' | p' == end = do
                   touchForeignPtr f
                   return []
              | otherwise = unsafeInterleaveIO $ do
                  e   <- peek p'
                  es  <- go (p' `advancePtr` incX)
                  return (e:es)
    in go p
{-# SPECIALIZE INLINE getElemsIOVector :: IOVector n Double -> IO [Double] #-}
{-# SPECIALIZE INLINE getElemsIOVector :: IOVector n (Complex Double) -> IO [Complex Double] #-}

getElemsIOVector' :: (Elem e) => IOVector n e -> IO [e]
getElemsIOVector' (IOVector f p n incX True) = do
    es <- getElemsIOVector' (IOVector f p n incX False)
    return $ map conjugate es    
getElemsIOVector' (IOVector f p n incX False) =
    let end = p `advancePtr` (-incX)
        go p' es | p' == end = do
                      touchForeignPtr f
                      return es
                 | otherwise = do
                      e <- peek p'
                      go (p' `advancePtr` (-incX)) (e:es)
    in go (p `advancePtr` ((n-1)*incX)) []
{-# SPECIALIZE INLINE getElemsIOVector' :: IOVector n Double -> IO [Double] #-}
{-# SPECIALIZE INLINE getElemsIOVector' :: IOVector n (Complex Double) -> IO [Complex Double] #-}

getAssocsIOVector :: (Elem e) => IOVector n e -> IO [(Int,e)]
getAssocsIOVector x = liftM2 zip (getIndicesIOVector x) (getElemsIOVector x)
{-# INLINE getAssocsIOVector #-}

getAssocsIOVector' :: (Elem e) => IOVector n e -> IO [(Int,e)]
getAssocsIOVector' x = liftM2 zip (getIndicesIOVector' x) (getElemsIOVector' x)
{-# INLINE getAssocsIOVector' #-}

unsafeReadElemIOVector :: (Elem e) => IOVector n e -> Int -> IO e
unsafeReadElemIOVector (IOVector f p n incX c) i
    | c = liftM conjugate $ unsafeReadElemIOVector (IOVector f p n incX False) i
    | otherwise = do
        e <- peekElemOff p (i*incX)
        touchForeignPtr f
        return e
{-# SPECIALIZE INLINE unsafeReadElemIOVector :: IOVector n Double -> Int -> IO (Double) #-}
{-# SPECIALIZE INLINE unsafeReadElemIOVector :: IOVector n (Complex Double) -> Int -> IO (Complex Double) #-}

canModifyElemIOVector :: IOVector n e -> Int -> IO Bool
canModifyElemIOVector _ _ = return True
{-# INLINE canModifyElemIOVector #-}

unsafeWriteElemIOVector :: (Elem e) => IOVector n e -> Int -> e -> IO ()
unsafeWriteElemIOVector (IOVector f p _ incX c) i e =
    let e' = if c then conjugate e else e
    in do
        pokeElemOff p (i*incX) e'
        touchForeignPtr f
{-# SPECIALIZE INLINE unsafeWriteElemIOVector :: IOVector n Double -> Int -> Double -> IO () #-}
{-# SPECIALIZE INLINE unsafeWriteElemIOVector :: IOVector n (Complex Double) -> Int -> Complex Double -> IO () #-}

unsafeModifyElemIOVector :: (Elem e) => IOVector n e -> Int -> (e -> e) -> IO ()
unsafeModifyElemIOVector (IOVector f p _ incX c) i g =
    let g' = if c then conjugate . g . conjugate else g
        p' = p `advancePtr` (i*incX)
    in do
        e <- peek p'
        poke p' (g' e)
        touchForeignPtr f
{-# SPECIALIZE INLINE unsafeModifyElemIOVector :: IOVector n Double -> Int -> (Double -> Double) -> IO () #-}
{-# SPECIALIZE INLINE unsafeModifyElemIOVector :: IOVector n (Complex Double) -> Int -> (Complex Double -> Complex Double) -> IO () #-}

unsafeSwapElemsIOVector :: (Elem e) => IOVector n e -> Int -> Int -> IO ()
unsafeSwapElemsIOVector (IOVector f p _ incX _) i1 i2 =
    let p1 = p `advancePtr` (i1*incX)
        p2 = p `advancePtr` (i2*incX)
    in do
        e1 <- peek p1
        e2 <- peek p2
        poke p2 e1
        poke p1 e2
        touchForeignPtr f
{-# SPECIALIZE INLINE unsafeSwapElemsIOVector :: IOVector n Double -> Int -> Int -> IO () #-}
{-# SPECIALIZE INLINE unsafeSwapElemsIOVector :: IOVector n (Complex Double) -> Int -> Int -> IO () #-}
        
                    
modifyWithIOVector :: (Elem e) => (e -> e) -> IOVector n e -> IO ()
modifyWithIOVector g (IOVector f p n incX c) =
    let g'  = if c then (conjugate . g . conjugate) else g
        end = p `advancePtr` (n*incX)
        go p' | p' == end = touchForeignPtr f
              | otherwise = do
                   e <- peek p'
                   poke p' (g' e)
                   go (p' `advancePtr` incX)
    in go p
{-# SPECIALIZE INLINE modifyWithIOVector :: (Double -> Double) -> IOVector n Double -> IO () #-}
{-# SPECIALIZE INLINE modifyWithIOVector :: (Complex Double -> Complex Double) -> IOVector n (Complex Double) -> IO () #-}

-- | Set every element in the vector to zero.
setZeroIOVector :: (Elem e) => IOVector n e -> IO ()
setZeroIOVector x@(IOVector f p n incX _)
    | incX == 1 = clearArray p n >> touchForeignPtr f
    | otherwise = setConstantIOVector 0 x
{-# INLINE setZeroIOVector #-}

-- | Set every element in the vector to a constant.
setConstantIOVector :: (Elem e) => e -> IOVector n e -> IO ()
setConstantIOVector 0 x | strideIOVector x == 1 = setZeroIOVector x
setConstantIOVector e (IOVector f p n incX c) =
    let e'   = if c then conjugate e else e
        end = p `advancePtr` (n*incX)
        go p' | p' == end = touchForeignPtr f
              | otherwise = do
                   poke p' e'
                   go (p' `advancePtr` incX)
    in go p
{-# INLINE setConstantIOVector #-}


doConjIOVector :: (BLAS1 e) => IOVector n e -> IO ()
doConjIOVector (IOVector f p n incX _) =
    BLAS.vconj n p incX >> touchForeignPtr f
{-# INLINE doConjIOVector #-}

scaleByIOVector :: (BLAS1 e) => e -> IOVector n e -> IO ()
scaleByIOVector 1 _ = return ()
scaleByIOVector k (IOVector f p n incX c) =
    let k' = if c then conjugate k else k
    in BLAS.scal n k' p incX >> touchForeignPtr f
{-# INLINE scaleByIOVector #-}
                    
shiftByIOVector :: (Elem e) => e -> IOVector n e -> IO ()                    
shiftByIOVector k x | isConjIOVector x = 
                        shiftByIOVector (conjugate k) (conjIOVector x)
                    | otherwise = 
                        modifyWithIOVector (k+) x
{-# INLINE shiftByIOVector #-}

instance Shaped IOVector Int e where
    shape = shapeIOVector
    {-# INLINE shape #-}
    bounds = boundsIOVector
    {-# INLINE bounds #-}

instance (Elem e) => ReadTensor IOVector Int e IO where
    getSize = getSizeIOVector
    {-# INLINE getSize #-}
    unsafeReadElem = unsafeReadElemIOVector
    {-# INLINE unsafeReadElem #-}
    getIndices = getIndicesIOVector
    {-# INLINE getIndices #-}
    getIndices' = getIndicesIOVector'
    {-# INLINE getIndices' #-}
    getElems = getElemsIOVector
    {-# INLINE getElems #-}
    getElems' = getElemsIOVector'
    {-# INLINE getElems' #-}
    getAssocs = getAssocsIOVector
    {-# INLINE getAssocs #-}
    getAssocs' = getAssocsIOVector'
    {-# INLINE getAssocs' #-}

instance (BLAS1 e) => WriteTensor IOVector Int e IO where
    getMaxSize = getMaxSizeIOVector
    {-# INLINE getMaxSize #-}
    setZero = setZeroIOVector
    {-# INLINE setZero #-}
    setConstant = setConstantIOVector
    {-# INLINE setConstant #-}
    canModifyElem = canModifyElemIOVector
    {-# INLINE canModifyElem #-}
    unsafeWriteElem = unsafeWriteElemIOVector
    {-# INLINE unsafeWriteElem #-}
    unsafeModifyElem = unsafeModifyElemIOVector
    {-# INLINE unsafeModifyElem #-}
    modifyWith = modifyWithIOVector
    {-# INLINE modifyWith #-}
    doConj = doConjIOVector
    {-# INLINE doConj #-}
    scaleBy = scaleByIOVector
    {-# INLINE scaleBy #-}
    shiftBy = shiftByIOVector
    {-# INLINE shiftBy #-}
