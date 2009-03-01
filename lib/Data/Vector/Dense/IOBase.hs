{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances, 
        GADTs #-}
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
import Text.Printf

import BLAS.Internal ( clearArray )
import BLAS.Types( ConjEnum(..), flipConj )
import Data.Elem.BLAS ( Complex, Elem, BLAS1, conjugate )
import qualified Data.Elem.BLAS.Base as BLAS
import qualified Data.Elem.BLAS.Level1 as BLAS

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

-- | Dense vectors in the 'IO' monad.  The type arguments are as follows:
--
--     * @n@: a phantom type for the dimension of the vector
--
--     * @e@: the element type of the vector.  Only certain element types
--       are supported.
--
data IOVector e = 
      Elem e =>
      IOVector {-# UNPACK #-} !ConjEnum       
               {-# UNPACK #-} !Int            -- the length of the vector
               {-# UNPACK #-} !(ForeignPtr e) -- memory owner
               {-# UNPACK #-} !(Ptr e)        -- pointer to the first element
               {-# UNPACK #-} !Int            -- the stride (in elements, not bytes) between elements.

-- | View an array in memory as a vector.
vectorViewArray :: (Elem e)
                => ForeignPtr e 
                -> Int          -- ^ offset
                -> Int          -- ^ length
                -> IOVector e
vectorViewArray = vectorViewArrayWithStride 1
{-# INLINE vectorViewArray #-}

-- | View an array in memory as a vector, with the given stride.
vectorViewArrayWithStride :: (Elem e)
                          => Int          -- ^ stride
                          -> ForeignPtr e
                          -> Int          -- ^ offset
                          -> Int          -- ^ length
                          -> IOVector e
vectorViewArrayWithStride s f o n =
    let p = unsafeForeignPtrToPtr f `advancePtr` o
    in IOVector NoConj n f p s
{-# INLINE vectorViewArrayWithStride #-}
                          
dimIOVector :: IOVector e -> Int
dimIOVector (IOVector _ n _ _ _) = n
{-# INLINE dimIOVector #-}

strideIOVector :: IOVector e -> Int
strideIOVector (IOVector _ _ _ _ incX) = incX
{-# INLINE strideIOVector #-}

conjEnumIOVector :: IOVector e -> ConjEnum
conjEnumIOVector (IOVector c _ _ _ _) = c
{-# INLINE conjEnumIOVector #-}

isConjIOVector :: IOVector e -> Bool
isConjIOVector x = conjEnumIOVector x == Conj
{-# INLINE isConjIOVector #-}

conjIOVector :: IOVector e -> IOVector e
conjIOVector (IOVector c n f p incX) = (IOVector (flipConj c) n f p incX)
{-# INLINE conjIOVector #-}

unsafeSubvectorViewWithStrideIOVector :: 
    Int -> IOVector e -> Int -> Int -> IOVector e
unsafeSubvectorViewWithStrideIOVector s' (IOVector c _ f p inc) o' n' =
    IOVector c n' f (p `advancePtr` (inc*o')) (inc*s')
{-# INLINE unsafeSubvectorViewWithStrideIOVector #-}

-- | Execute an 'IO' action with a pointer to the first element in the
-- vector.
withIOVector :: IOVector e -> (Ptr e -> IO a) -> IO a
withIOVector (IOVector _ _ f p _) g = do
    a <- g p
    touchForeignPtr f
    return a
{-# INLINE withIOVector #-}

newIOVector_ :: (Elem e) => Int -> IO (IOVector e)
newIOVector_ n
    | n < 0 = 
        fail $ "Tried to create a vector with `" ++ show n ++ "' elements."
    | otherwise = do
        arr <- mallocForeignPtrArray n
        return $ IOVector NoConj n arr (unsafeForeignPtrToPtr arr) 1

newIOVector :: (Elem e) => Int -> [(Int,e)] -> IO (IOVector e)
newIOVector n ies = 
    let go p ((i,e):ies') = do
            when (i < 0 || i >= n) $ error $
                printf "newVector %d [ ..., (%d,_), ... ]: invalid index" n i
            pokeElemOff p i e
            go p ies'
        go _ [] = return ()
    in do
        x  <- newZeroIOVector n
        withIOVector x $ \p -> go p ies
        return x

unsafeNewIOVector :: (Elem e) => Int -> [(Int,e)] -> IO (IOVector e)
unsafeNewIOVector n ies =
    let go p ((i,e):ies') = do
            pokeElemOff p i e
            go p ies'
        go _ [] = return ()
    in do
        x  <- newZeroIOVector n
        withIOVector x $ \p -> go p ies
        return x

newListIOVector :: (Elem e) => Int -> [e] -> IO (IOVector e)
newListIOVector n es = do
    x  <- newIOVector_ n
    withIOVector x $ \p -> do
        pokeArray p $ take n $ es ++ (repeat 0)
    return x

newZeroIOVector :: (Elem e) => Int -> IO (IOVector e)
newZeroIOVector n = do
    x  <- newIOVector_ n
    setZeroIOVector x
    return x

newConstantIOVector :: (Elem e) => Int -> e -> IO (IOVector e)
newConstantIOVector n e = do
    x  <- newIOVector_ n
    setConstantIOVector e x
    return x

newBasisIOVector :: (Elem e) => Int -> Int -> IO (IOVector e)
newBasisIOVector n i = do
    x  <- newZeroIOVector n
    setBasisIOVector i x
    return x

setBasisIOVector :: Int -> IOVector e -> IO ()
setBasisIOVector i x@(IOVector _ _ _ _ _)
    | i < 0 || i >= n =
        error $ printf "setBasisIOVector %d <vector of dim %d>: invalid index" 
                       i n
    | otherwise = do
        setZeroIOVector x
        unsafeWriteElemIOVector x i 1 
  where
    n = dimIOVector x

newCopyIOVector :: IOVector e -> IO (IOVector e)
newCopyIOVector (IOVector c n f p incX) = do
    (IOVector _ _ f' p' _) <- newIOVector_ n
    BLAS.copy c c n p incX p' 1
    touchForeignPtr f
    touchForeignPtr f'
    return (IOVector c n f' p' 1)

newCopyIOVector' :: IOVector e -> IO (IOVector e)
newCopyIOVector' x@(IOVector _ _ _ _ _) = do
    y  <- newIOVector_ (dimIOVector x)
    unsafeCopyIOVector y x
    return y

unsafeCopyIOVector :: IOVector e -> IOVector e -> IO ()
unsafeCopyIOVector (IOVector cy n fy py incy) 
                   (IOVector cx _ fx px incx) = do
    io <- BLAS.copy cx cy n px incx py incy
    touchForeignPtr fx
    touchForeignPtr fy
    return io
{-# INLINE unsafeCopyIOVector #-}

unsafeSwapIOVector :: IOVector e -> IOVector e -> IO ()
unsafeSwapIOVector (IOVector cx n fx px incx) 
                   (IOVector cy _ fy py incy) = do
    io <- BLAS.swap cx cy n px incx py incy
    touchForeignPtr fx
    touchForeignPtr fy
    return io
{-# INLINE unsafeSwapIOVector #-}

shapeIOVector :: IOVector e -> Int
shapeIOVector = dimIOVector
{-# INLINE shapeIOVector #-}

boundsIOVector :: IOVector e -> (Int,Int)
boundsIOVector x = (0, dimIOVector x - 1)
{-# INLINE boundsIOVector #-}

sizeIOVector :: IOVector e -> Int
sizeIOVector = dimIOVector
{-# INLINE sizeIOVector #-}

getSizeIOVector :: IOVector e -> IO Int
getSizeIOVector = return . sizeIOVector
{-# INLINE getSizeIOVector #-}

getMaxSizeIOVector :: IOVector e -> IO Int
getMaxSizeIOVector = getSizeIOVector
{-# INLINE getMaxSizeIOVector #-}

indicesIOVector :: IOVector e -> [Int]
indicesIOVector x = [ 0..n-1 ] where n = dimIOVector x
{-# INLINE indicesIOVector #-}

getIndicesIOVector :: IOVector e -> IO [Int]
getIndicesIOVector = return . indicesIOVector
{-# INLINE getIndicesIOVector #-}

getIndicesIOVector' :: IOVector e -> IO [Int]
getIndicesIOVector' = getIndicesIOVector
{-# INLINE getIndicesIOVector' #-}

getElemsIOVector :: IOVector e -> IO [e]
getElemsIOVector (IOVector Conj n f p incX) = do
    es <- getElemsIOVector (IOVector NoConj n f p incX)
    return $ map conjugate es
getElemsIOVector (IOVector NoConj n f p incX) =
    let end = p `advancePtr` (n*incX)
        go p' | p' == end = do
                   touchForeignPtr f
                   return []
              | otherwise = unsafeInterleaveIO $ do
                  e   <- peek p'
                  es  <- go (p' `advancePtr` incX)
                  return (e:es)
    in go p
{-# SPECIALIZE INLINE getElemsIOVector :: IOVector Double -> IO [Double] #-}
{-# SPECIALIZE INLINE getElemsIOVector :: IOVector (Complex Double) -> IO [Complex Double] #-}

getElemsIOVector' :: IOVector e -> IO [e]
getElemsIOVector' (IOVector Conj n f p incX) = do
    es <- getElemsIOVector' (IOVector NoConj n f p incX)
    return $ map conjugate es    
getElemsIOVector' (IOVector NoConj n f p incX) =
    let end = p `advancePtr` (-incX)
        go p' es | p' == end = do
                      touchForeignPtr f
                      return es
                 | otherwise = do
                      e <- peek p'
                      go (p' `advancePtr` (-incX)) (e:es)
    in go (p `advancePtr` ((n-1)*incX)) []
{-# SPECIALIZE INLINE getElemsIOVector' :: IOVector Double -> IO [Double] #-}
{-# SPECIALIZE INLINE getElemsIOVector' :: IOVector (Complex Double) -> IO [Complex Double] #-}

getAssocsIOVector :: IOVector e -> IO [(Int,e)]
getAssocsIOVector x = liftM2 zip (getIndicesIOVector x) (getElemsIOVector x)
{-# INLINE getAssocsIOVector #-}

getAssocsIOVector' :: IOVector e -> IO [(Int,e)]
getAssocsIOVector' x = liftM2 zip (getIndicesIOVector' x) (getElemsIOVector' x)
{-# INLINE getAssocsIOVector' #-}

unsafeReadElemIOVector :: IOVector e -> Int -> IO e
unsafeReadElemIOVector (IOVector Conj   n f p incX) i = 
    liftM conjugate $ unsafeReadElemIOVector (IOVector NoConj n f p incX) i
unsafeReadElemIOVector (IOVector NoConj _ f p incX) i = do
    e <- peekElemOff p (i*incX)
    touchForeignPtr f
    return e
{-# SPECIALIZE INLINE unsafeReadElemIOVector :: IOVector Double -> Int -> IO (Double) #-}
{-# SPECIALIZE INLINE unsafeReadElemIOVector :: IOVector (Complex Double) -> Int -> IO (Complex Double) #-}

canModifyElemIOVector :: IOVector e -> Int -> IO Bool
canModifyElemIOVector _ _ = return True
{-# INLINE canModifyElemIOVector #-}

unsafeWriteElemIOVector :: IOVector e -> Int -> e -> IO ()
unsafeWriteElemIOVector (IOVector c _ f p incX) i e =
    let e' = if c == Conj then conjugate e else e
    in do
        io <- pokeElemOff p (i*incX) e'
        touchForeignPtr f
        return io
{-# SPECIALIZE INLINE unsafeWriteElemIOVector :: IOVector Double -> Int -> Double -> IO () #-}
{-# SPECIALIZE INLINE unsafeWriteElemIOVector :: IOVector (Complex Double) -> Int -> Complex Double -> IO () #-}

unsafeModifyElemIOVector :: IOVector e -> Int -> (e -> e) -> IO ()
unsafeModifyElemIOVector (IOVector c _ f p incX) i g =
    let g' = if c == Conj then conjugate . g . conjugate else g
        p' = p `advancePtr` (i*incX)
    in do
        e  <- peek p'
        io <- poke p' (g' e)
        touchForeignPtr f
        return io
{-# SPECIALIZE INLINE unsafeModifyElemIOVector :: IOVector Double -> Int -> (Double -> Double) -> IO () #-}
{-# SPECIALIZE INLINE unsafeModifyElemIOVector :: IOVector (Complex Double) -> Int -> (Complex Double -> Complex Double) -> IO () #-}

unsafeSwapElemsIOVector :: IOVector e -> Int -> Int -> IO ()
unsafeSwapElemsIOVector (IOVector _ _ f p incX) i1 i2 =
    let p1 = p `advancePtr` (i1*incX)
        p2 = p `advancePtr` (i2*incX)
    in do
        e1  <- peek p1
        e2  <- peek p2
        poke p2 e1
        poke p1 e2
        touchForeignPtr f
{-# SPECIALIZE INLINE unsafeSwapElemsIOVector :: IOVector Double -> Int -> Int -> IO () #-}
{-# SPECIALIZE INLINE unsafeSwapElemsIOVector :: IOVector (Complex Double) -> Int -> Int -> IO () #-}
                            
modifyWithIOVector :: (e -> e) -> IOVector e -> IO ()
modifyWithIOVector g (IOVector c n f p incX) =
    let g'  = if c == Conj then (conjugate . g . conjugate) else g
        end = p `advancePtr` (n*incX)
        go p' | p' == end = touchForeignPtr f
              | otherwise = do
                   e  <- peek p'
                   poke p' (g' e)
                   go (p' `advancePtr` incX)
    in go p
{-# SPECIALIZE INLINE modifyWithIOVector :: (Double -> Double) -> IOVector Double -> IO () #-}
{-# SPECIALIZE INLINE modifyWithIOVector :: (Complex Double -> Complex Double) -> IOVector (Complex Double) -> IO () #-}

setZeroIOVector :: IOVector e -> IO ()
setZeroIOVector x@(IOVector _ n f p incX)
    | incX == 1 = do
        io <- clearArray p n
        touchForeignPtr f
        return io
    | otherwise = 
        setConstantIOVector 0 x
{-# INLINE setZeroIOVector #-}

setConstantIOVector :: e -> IOVector e -> IO ()
setConstantIOVector e x@(IOVector _ _ _ _ _) 
    | e == 0 && strideIOVector x == 1 = setZeroIOVector x
setConstantIOVector e (IOVector c n f p incX) =
    let e'   = if c == Conj then conjugate e else e
        end = p `advancePtr` (n*incX)
        go p' | p' == end = touchForeignPtr f
              | otherwise = do
                   poke p' e'
                   go (p' `advancePtr` incX)
    in go p
{-# INLINE setConstantIOVector #-}

doConjIOVector :: IOVector e -> IO ()
doConjIOVector (IOVector _ n f p incX) = do
    io <- BLAS.vconj n p incX
    touchForeignPtr f
    return io
{-# INLINE doConjIOVector #-}

scaleByIOVector :: (BLAS1 e) => e -> IOVector e -> IO ()
scaleByIOVector 1 _ = return ()
scaleByIOVector k (IOVector c n f p incX) =
    let k' = if c == Conj then conjugate k else k
    in do
        io <- BLAS.scal n k' p incX
        touchForeignPtr f
        return io
{-# INLINE scaleByIOVector #-}
                    
shiftByIOVector :: e -> IOVector e -> IO ()                    
shiftByIOVector k x@(IOVector _ _ _ _ _) 
    | isConjIOVector x = 
         shiftByIOVector (conjugate k) (conjIOVector x)
    | otherwise = 
         modifyWithIOVector (k+) x
{-# INLINE shiftByIOVector #-}

instance Shaped IOVector Int where
    shape = shapeIOVector
    {-# INLINE shape #-}
    bounds = boundsIOVector
    {-# INLINE bounds #-}

instance ReadTensor IOVector Int IO where
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

instance WriteTensor IOVector Int IO where
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
