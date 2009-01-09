{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances, TypeFamilies, 
        Rank2Types #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.IOBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.IOBase
    where

import Control.Monad
import Foreign
import System.IO.Unsafe

import BLAS.Internal( diagLen )

import Data.Elem.BLAS( Complex, Elem, BLAS1, conjugate )
import qualified Data.Elem.BLAS.Level1 as BLAS

import Data.Matrix.Class

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Vector.Dense.IOBase


-- | Dense matrix in the 'IO' monad.  The type arguments are as follows:
--
--     * @np@: a phantom type for the shape of the matrix.  Most functions
--       will demand that this be specified as a pair.  When writing a function
--       signature, you should always prefer @IOMatrix (n,p) e@ to
--       @IOMatrix np e@.
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
data IOMatrix np e =
      IOMatrix {-# UNPACK #-} !(ForeignPtr e) -- a pointer to the storage region
               {-# UNPACK #-} !(Ptr e)        -- a pointer to the first element
               {-# UNPACK #-} !Int            -- the number of rows in the matrix
               {-# UNPACK #-} !Int            -- the number of colunps in the matrix
               {-# UNPACK #-} !Int            -- the leading dimension size of the matrix
               {-# UNPACK #-} !Bool           -- indicates whether or not the matrix is transposed and conjugated

-- | View an array in memory as a matrix.
matrixViewArray :: (Elem e)
                => ForeignPtr e
                -> Int          -- ^ offset
                -> (Int,Int)    -- ^ shape
                -> IOMatrix (n,p) e
matrixViewArray f o (m,n) = matrixViewArrayWithLda m f o (m,n)
{-# INLINE matrixViewArray #-}

-- | View an array in memory as a matrix, with the given leading dimension
-- size.
matrixViewArrayWithLda :: (Elem e)
                       => Int          -- ^ leading dimension size
                       -> ForeignPtr e
                       -> Int          -- ^ offset
                       -> (Int,Int)    -- ^ shape
                       -> IOMatrix (n,p) e
matrixViewArrayWithLda l f o (m,n) =
    let p = unsafeForeignPtrToPtr f `advancePtr` o
    in IOMatrix f p m n l False
{-# INLINE matrixViewArrayWithLda #-}


numRowsIOMatrix :: IOMatrix np e -> Int
numRowsIOMatrix (IOMatrix _ _ m _ _ _) = m
{-# INLINE numRowsIOMatrix #-}

numColsIOMatrix :: IOMatrix np e -> Int
numColsIOMatrix (IOMatrix _ _ _ n _ _) = n
{-# INLINE numColsIOMatrix #-}

ldaMatrixIOMatrix :: IOMatrix np e -> Int
ldaMatrixIOMatrix (IOMatrix _ _ _ _ l _) = l
{-# INLINE ldaMatrixIOMatrix #-}

isHermIOMatrix :: IOMatrix np e -> Bool
isHermIOMatrix (IOMatrix _ _ _ _ _ h) = h
{-# INLINE isHermIOMatrix #-}

hermIOMatrix :: IOMatrix np e -> IOMatrix nm e
hermIOMatrix (IOMatrix f p m n l h) = (IOMatrix f p n m l (not h))
{-# INLINE hermIOMatrix #-}

unsafeSubmatrixViewIOMatrix :: (Elem e) =>
    IOMatrix np e -> (Int,Int) -> (Int,Int) -> IOMatrix np' e
unsafeSubmatrixViewIOMatrix (IOMatrix f p _ _ l h) (i,j) (m',n') =
    let o = if h then i*l+j else i+j*l
        p' = p `advancePtr` o
    in IOMatrix f p' m' n' l h
{-# INLINE unsafeSubmatrixViewIOMatrix #-}

unsafeRowViewIOMatrix :: (Elem e) => IOMatrix np e -> Int -> IOVector p e
unsafeRowViewIOMatrix (IOMatrix f p _ n l h) i =
    let (o,s) = if h then (i*l,1) else (i,l)
        p'    = p `advancePtr` o
    in IOVector f p' n s h
{-# INLINE unsafeRowViewIOMatrix #-}

unsafeColViewIOMatrix :: (Elem e) => IOMatrix np e -> Int -> IOVector n e
unsafeColViewIOMatrix (IOMatrix f p m _ l h) j =
    let (o,s) = if h then (j,l) else (j*l,1)
        p'    = p `advancePtr` o
    in IOVector f p' m s h
{-# INLINE unsafeColViewIOMatrix #-}

unsafeDiagViewIOMatrix :: (Elem e) => IOMatrix np e -> Int -> IOVector k e
unsafeDiagViewIOMatrix (IOMatrix f p m n l h) i =
    let o = if i >= 0 
                then if h then  i   else  i*l
                else if h then -i*l else -i
        p' = p `advancePtr` o
        k  = diagLen (m,n) i
        s  = l+1
    in IOVector f p' k s h
{-# INLINE unsafeDiagViewIOMatrix #-}

maybeViewVectorAsRowIOMatrix :: (Elem e) => IOVector p e -> Maybe (IOMatrix p1 e)
maybeViewVectorAsRowIOMatrix (IOVector f p n s c)
    | c && (s == 1) =
        Just $ IOMatrix f p 1 n (max 1 n) True
    | not c =
        Just $ IOMatrix f p 1 n s         False
    | otherwise =
        Nothing
{-# INLINE maybeViewVectorAsRowIOMatrix #-}

maybeViewVectorAsColIOMatrix :: (Elem e) => IOVector n e -> Maybe (IOMatrix n1 e)
maybeViewVectorAsColIOMatrix (IOVector f p n s c)
    | c =
        Just $ IOMatrix f p n 1 s         True
    | s == 1 =
        Just $ IOMatrix f p n 1 (max 1 n) False
    | otherwise =
        Nothing
{-# INLINE maybeViewVectorAsColIOMatrix #-}

maybeMatrixToVectorViewIOMatrix :: (Elem e) => IOMatrix np e -> Maybe (IOVector k e)
maybeMatrixToVectorViewIOMatrix (IOMatrix f p m n l h)
    | h = 
        liftM conjIOVector $ maybeMatrixToVectorViewIOMatrix (IOMatrix f p n m l False)
    | l == m =
        Just $ IOVector f p (m*n) 1 False
    | m == 1 =
        Just $ IOVector f p n     l False
    | otherwise =
        Nothing
{-# INLINE maybeMatrixToVectorViewIOMatrix #-}

liftIOMatrix :: (Elem e) => (forall n. IOVector n e -> IO ()) -> IOMatrix np e -> IO ()
liftIOMatrix g (IOMatrix f p m n l h)
    | h && (l == n) =
        g (IOVector f p (m*n) 1 True)
    | (not h) && (l == m) =
        g (IOVector f p (m*n) 1 False)
    | otherwise =
        let (m',n') = if h then (n,m) else (m,n)
            end     = p `advancePtr` (n'*l)
            go p' | p' == end = return ()
                  | otherwise = do
                      g (IOVector f p' m' 1 h)
                      go (p' `advancePtr` l)
        in go p

liftIOMatrix2 :: (Elem e, Elem f) =>
    (forall k. IOVector k e -> IOVector k f -> IO ()) ->
        IOMatrix np e -> IOMatrix np f -> IO ()
liftIOMatrix2 f a b =
    if isHermIOMatrix a == isHermIOMatrix b
        then case (maybeMatrixToVectorViewIOMatrix a, maybeMatrixToVectorViewIOMatrix b) of
                 ((Just x), (Just y)) -> f x y
                 _                    -> elementwise
        else elementwise
  where
    
    elementwise =
        let vecsA = if isHermIOMatrix a then rowViews
                                        else colViews
            vecsB = if isHermIOMatrix a then rowViews
                                        else colViews
            xs = vecsA a
            ys = vecsB b
        in zipWithM_ f xs ys
        
    rowViews c = [ unsafeRowViewIOMatrix c i | i <- [ 0..numRows c ] ]
    colViews c = [ unsafeColViewIOMatrix c j | j <- [ 0..numCols c ] ]

-- | Perform an 'IO' action with a pointer to the first element of the
-- matrix.
withIOMatrix :: IOMatrix (n,p) e -> (Ptr e -> IO a) -> IO a
withIOMatrix (IOMatrix f p _ _ _ _) g = do
    a <- g p
    touchForeignPtr f
    return a
{-# INLINE withIOMatrix #-}

-- | Create a new matrix of given shape, but do not initialize the elements.
newIOMatrix_ :: (Elem e) => (Int,Int) -> IO (IOMatrix np e)
newIOMatrix_ (m,n) 
    | m < 0 || n < 0 =
        fail $ 
            "Tried to create a matrix with shape `" ++ show (m,n) ++ "'"
    | otherwise =  do
        f <- mallocForeignPtrArray (m*n)
        return (IOMatrix f (unsafeForeignPtrToPtr f) m n (max 1 m) False)
{-# INLINE newIOMatrix_ #-}

newCopyIOMatrix :: (BLAS1 e) => IOMatrix np e -> IO (IOMatrix np e)
newCopyIOMatrix (IOMatrix f p m n l h) = 
    let (m',n') = if h then (n,m) else (m,n)
        l'      = max 1 m'
    in do
        (IOMatrix f' p' _ _ _ _) <- newIOMatrix_ (m',n')
        if l == m'
            then do
                BLAS.copy (m*n) p 1 p' 1
            else 
                let go src dst i | i == n'   = return ()
                                 | otherwise = do
                        BLAS.copy m' src 1 dst 1
                        go (src `advancePtr` l) (dst `advancePtr` l') (i+1)
                in go p p' 0
        touchForeignPtr f
        touchForeignPtr f'
        return (IOMatrix f' p' m n l' h)
{-# INLINE newCopyIOMatrix #-}

shapeIOMatrix :: IOMatrix np e -> (Int,Int)
shapeIOMatrix (IOMatrix _ _ m n _ _) = (m,n)
{-# INLINE shapeIOMatrix #-}

boundsIOMatrix :: IOMatrix np e -> ((Int,Int), (Int,Int))
boundsIOMatrix a = ((0,0), (m-1,n-1)) where (m,n) = shapeIOMatrix a
{-# INLINE boundsIOMatrix #-}

sizeIOMatrix :: IOMatrix np e -> Int
sizeIOMatrix (IOMatrix _ _ m n _ _) = m*n
{-# INLINE sizeIOMatrix #-}

getSizeIOMatrix :: IOMatrix np e -> IO Int
getSizeIOMatrix = return . sizeIOMatrix
{-# INLINE getSizeIOMatrix #-}

getMaxSizeIOMatrix :: IOMatrix np e -> IO Int
getMaxSizeIOMatrix = getSizeIOMatrix
{-# INLINE getMaxSizeIOMatrix #-}

indicesIOMatrix :: IOMatrix np e -> [(Int,Int)]
indicesIOMatrix (IOMatrix _ _ m n _ h)
    | h         = [ (i,j) | i <- [ 0..m-1 ], j <- [ 0..n-1 ] ]
    | otherwise = [ (i,j) | j <- [ 0..n-1 ], i <- [ 0..m-1 ] ]
{-# INLINE indicesIOMatrix #-}

getIndicesIOMatrix :: IOMatrix np e -> IO [(Int,Int)]
getIndicesIOMatrix = return . indicesIOMatrix
{-# INLINE getIndicesIOMatrix #-}

getIndicesIOMatrix' :: IOMatrix np e -> IO [(Int,Int)]
getIndicesIOMatrix' = getIndicesIOMatrix
{-# INLINE getIndicesIOMatrix' #-}

getElemsIOMatrix :: (Elem e) => IOMatrix np e -> IO [e]
getElemsIOMatrix (IOMatrix f p m n l h) 
    | h         = liftM (map conjugate) $ 
                      getElemsIOMatrix (IOMatrix f p n m l False)
    | l == m    = getElemsIOVector (IOVector f p (m*n) 1 False)
    | otherwise =
        let end   = p `advancePtr` (n*l)
            go p' | p' == end = return []
                  | otherwise = unsafeInterleaveIO $ do
                        c  <- getElemsIOVector (IOVector f p' m 1 False)
                        cs <- go (p' `advancePtr` l)
                        return (c ++ cs)
        in go p
{-# SPECIALIZE INLINE getElemsIOMatrix :: IOMatrix np Double -> IO [Double] #-}
{-# SPECIALIZE INLINE getElemsIOMatrix :: IOMatrix np (Complex Double) -> IO [Complex Double] #-}

getElemsIOMatrix' :: (Elem e) => IOMatrix np e -> IO [e]
getElemsIOMatrix' (IOMatrix f p m n l h) 
    | h         = liftM (map conjugate) $ 
                      getElemsIOMatrix' (IOMatrix f p n m l False)
    | l == m    = getElemsIOVector' (IOVector f p (m*n) 1 False)
    | otherwise =
        let end   = p `advancePtr` (n*l)
            go p' | p' == end = return []
                  | otherwise = do
                        c  <- getElemsIOVector' (IOVector f p' m 1 False)
                        cs <- go (p' `advancePtr` l)
                        return (c ++ cs)
        in go p
{-# SPECIALIZE INLINE getElemsIOMatrix' :: IOMatrix np Double -> IO [Double] #-}
{-# SPECIALIZE INLINE getElemsIOMatrix' :: IOMatrix np (Complex Double) -> IO [Complex Double] #-}

getAssocsIOMatrix :: (Elem e) => IOMatrix np e -> IO [((Int,Int),e)]
getAssocsIOMatrix a = do
    is <- getIndicesIOMatrix a
    es <- getElemsIOMatrix a
    return $ zip is es
{-# INLINE getAssocsIOMatrix #-}

getAssocsIOMatrix' :: (Elem e) => IOMatrix np e -> IO [((Int,Int),e)]
getAssocsIOMatrix' a = do
    is <- getIndicesIOMatrix' a
    es <- getElemsIOMatrix' a
    return $ zip is es
{-# INLINE getAssocsIOMatrix' #-}

unsafeReadElemIOMatrix :: (Elem e) => IOMatrix np e -> (Int,Int) -> IO e
unsafeReadElemIOMatrix (IOMatrix f p _ _ l h) (i,j)
    | h = do
        e <- liftM conjugate $ peekElemOff p (i*l+j)
        touchForeignPtr f
        return e
    | otherwise = do
        e <- peekElemOff p (i+j*l)
        touchForeignPtr f
        return e
{-# SPECIALIZE INLINE unsafeReadElemIOMatrix :: IOMatrix n Double -> (Int,Int) -> IO (Double) #-}
{-# SPECIALIZE INLINE unsafeReadElemIOMatrix :: IOMatrix n (Complex Double) -> (Int,Int) -> IO (Complex Double) #-}        

canModifyElemIOMatrix :: IOMatrix np e -> (Int,Int) -> IO Bool
canModifyElemIOMatrix _ _ = return True
{-# INLINE canModifyElemIOMatrix #-}

unsafeWriteElemIOMatrix :: (Elem e) => 
    IOMatrix np e -> (Int,Int) -> e -> IO ()
unsafeWriteElemIOMatrix (IOMatrix f p _ _ l h) (i,j) e
    | h = do
        pokeElemOff p (i*l+j) (conjugate e)
        touchForeignPtr f
    | otherwise = do
        pokeElemOff p (i+j*l) e
        touchForeignPtr f
{-# SPECIALIZE INLINE unsafeWriteElemIOMatrix :: IOMatrix n Double -> (Int,Int) -> Double -> IO () #-}
{-# SPECIALIZE INLINE unsafeWriteElemIOMatrix :: IOMatrix n (Complex Double) -> (Int,Int) -> Complex Double -> IO () #-}        

unsafeModifyElemIOMatrix :: (Elem e) => 
    IOMatrix n e -> (Int,Int) -> (e -> e) -> IO ()
unsafeModifyElemIOMatrix (IOMatrix f p _ _ l h) (i,j) g =
    let g' = if h then conjugate . g . conjugate else g
        p' = if h then p `advancePtr` (i*l+j) else p `advancePtr` (i+j*l)
    in do
        e <- peek p'
        poke p' (g' e)
        touchForeignPtr f
{-# SPECIALIZE INLINE unsafeModifyElemIOMatrix :: IOMatrix n Double -> (Int,Int) -> (Double -> Double) -> IO () #-}
{-# SPECIALIZE INLINE unsafeModifyElemIOMatrix :: IOMatrix n (Complex Double) -> (Int,Int) -> (Complex Double -> Complex Double) -> IO () #-}

unsafeSwapElemsIOMatrix :: (Elem e) => 
    IOMatrix n e -> (Int,Int) -> (Int,Int) -> IO ()
unsafeSwapElemsIOMatrix (IOMatrix f p _ _ l h) (i1,j1) (i2,j2) =
    let (p1,p2) = if h then (p `advancePtr` (i1*l+j1), p `advancePtr` (i2*l+j2))
                       else (p `advancePtr` (i1+j1*l), p `advancePtr` (i2+j2*l))
    in do
        e1 <- peek p1
        e2 <- peek p2
        poke p2 e1
        poke p1 e2
        touchForeignPtr f
{-# SPECIALIZE INLINE unsafeSwapElemsIOMatrix :: IOMatrix n Double -> (Int,Int) -> (Int,Int) -> IO () #-}
{-# SPECIALIZE INLINE unsafeSwapElemsIOMatrix :: IOMatrix n (Complex Double) -> (Int,Int) -> (Int,Int) -> IO () #-}

modifyWithIOMatrix :: (Elem e) => (e -> e) -> IOMatrix np e -> IO ()
modifyWithIOMatrix g = liftIOMatrix (modifyWithIOVector g)
{-# INLINE modifyWithIOMatrix #-}

setZeroIOMatrix :: (Elem e) => IOMatrix np e -> IO ()
setZeroIOMatrix = liftIOMatrix setZeroIOVector
{-# INLINE setZeroIOMatrix #-}

setConstantIOMatrix :: (Elem e) => e -> IOMatrix np e -> IO ()
setConstantIOMatrix k = liftIOMatrix (setConstantIOVector k)
{-# INLINE setConstantIOMatrix #-}

doConjIOMatrix :: (BLAS1 e) => IOMatrix np e -> IO ()
doConjIOMatrix = liftIOMatrix doConjIOVector
{-# INLINE doConjIOMatrix #-}

scaleByIOMatrix :: (BLAS1 e) => e -> IOMatrix np e -> IO ()
scaleByIOMatrix k = liftIOMatrix (scaleByIOVector k)
{-# INLINE scaleByIOMatrix #-}
                    
shiftByIOMatrix :: (Elem e) => e -> IOMatrix np e -> IO ()                    
shiftByIOMatrix k = liftIOMatrix (shiftByIOVector k)
{-# INLINE shiftByIOMatrix #-}

instance Shaped IOMatrix (Int,Int) e where
    shape = shapeIOMatrix
    {-# INLINE shape #-}
    bounds = boundsIOMatrix
    {-# INLINE bounds #-}

instance (Elem e) => ReadTensor IOMatrix (Int,Int) e IO where
    getSize = getSizeIOMatrix
    {-# INLINE getSize #-}
    unsafeReadElem = unsafeReadElemIOMatrix
    {-# INLINE unsafeReadElem #-}
    getIndices = getIndicesIOMatrix
    {-# INLINE getIndices #-}
    getIndices' = getIndicesIOMatrix'
    {-# INLINE getIndices' #-}
    getElems = getElemsIOMatrix
    {-# INLINE getElems #-}
    getElems' = getElemsIOMatrix'
    {-# INLINE getElems' #-}
    getAssocs = getAssocsIOMatrix
    {-# INLINE getAssocs #-}
    getAssocs' = getAssocsIOMatrix'
    {-# INLINE getAssocs' #-}

instance (BLAS1 e) => WriteTensor IOMatrix (Int,Int) e IO where
    getMaxSize = getMaxSizeIOMatrix
    {-# INLINE getMaxSize #-}
    setZero = setZeroIOMatrix
    {-# INLINE setZero #-}
    setConstant = setConstantIOMatrix
    {-# INLINE setConstant #-}
    canModifyElem = canModifyElemIOMatrix
    {-# INLINE canModifyElem #-}
    unsafeWriteElem = unsafeWriteElemIOMatrix
    {-# INLINE unsafeWriteElem #-}
    unsafeModifyElem = unsafeModifyElemIOMatrix
    {-# INLINE unsafeModifyElem #-}
    modifyWith = modifyWithIOMatrix
    {-# INLINE modifyWith #-}
    doConj = doConjIOMatrix
    {-# INLINE doConj #-}
    scaleBy = scaleByIOMatrix
    {-# INLINE scaleBy #-}
    shiftBy = shiftByIOMatrix
    {-# INLINE shiftBy #-}

instance HasVectorView IOMatrix where
    type VectorView IOMatrix = IOVector

instance MatrixShaped IOMatrix e where
    herm = hermIOMatrix
    {-# INLINE herm #-}
