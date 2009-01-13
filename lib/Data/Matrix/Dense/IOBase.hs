{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances, FlexibleContexts,
        TypeFamilies, Rank2Types, GADTs #-}
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
import Text.Printf

import BLAS.Internal( diagLen )

import Data.Elem.BLAS
import qualified Data.Elem.BLAS.Base   as BLAS
import qualified Data.Elem.BLAS.Level1 as BLAS

import Data.Matrix.Class

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Vector.Dense.IOBase
import Data.Vector.Dense.Base


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
      Elem e =>
      IOMatrix { transEnumIOMatrix :: {-# UNPACK #-} !TransEnum
               , numRowsIOMatrix   :: {-# UNPACK #-} !Int
               , numColsIOMatrix   :: {-# UNPACK #-} !Int
               , fptrIOMatrix      :: {-# UNPACK #-} !(ForeignPtr e)
               , ptrIOMatrix       :: {-# UNPACK #-} !(Ptr e)
               , ldaIOMatrix       :: {-# UNPACK #-} !Int
               }

-- | View an array in memory as a matrix.
matrixViewArray :: (Elem e)
                => ForeignPtr e
                -> Int          -- ^ offset
                -> (Int,Int)    -- ^ shape
                -> IOMatrix (n,p) e
matrixViewArray f o (m,n) = matrixViewArrayWithLda (max 1 m) f o (m,n)
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
    in IOMatrix NoTrans m n f p l
{-# INLINE matrixViewArrayWithLda #-}

isHermIOMatrix :: IOMatrix np e -> Bool
isHermIOMatrix = (ConjTrans ==) . transEnumIOMatrix
{-# INLINE isHermIOMatrix #-}

hermIOMatrix :: IOMatrix np e -> IOMatrix nm e
hermIOMatrix (IOMatrix h m n f p l) = (IOMatrix (flipTrans h) n m f p l)
{-# INLINE hermIOMatrix #-}

unsafeSubmatrixViewIOMatrix :: 
    IOMatrix np e -> (Int,Int) -> (Int,Int) -> IOMatrix np' e
unsafeSubmatrixViewIOMatrix (IOMatrix h _ _ f p l) (i,j) (m',n') =
    let o = if h == ConjTrans then i*l+j else i+j*l
        p' = p `advancePtr` o
    in IOMatrix h m' n' f p' l
{-# INLINE unsafeSubmatrixViewIOMatrix #-}

unsafeRowViewIOMatrix :: IOMatrix np e -> Int -> IOVector p e
unsafeRowViewIOMatrix (IOMatrix h _ n f p l) i =
    let (c,o,s) = if h == ConjTrans then (Conj,i*l,1) else (NoConj,i,l)
        p'      = p `advancePtr` o
    in IOVector c n f p' s
{-# INLINE unsafeRowViewIOMatrix #-}

unsafeColViewIOMatrix :: IOMatrix np e -> Int -> IOVector n e
unsafeColViewIOMatrix (IOMatrix h m _ f p l) j =
    let (c,o,s) = if h == ConjTrans then (Conj,j,l) else (NoConj,j*l,1)
        p'      = p `advancePtr` o
    in IOVector c m f p' s
{-# INLINE unsafeColViewIOMatrix #-}

unsafeDiagViewIOMatrix :: IOMatrix np e -> Int -> IOVector k e
unsafeDiagViewIOMatrix (IOMatrix h m n f p l) i =
    let o = if i >= 0 
                then if h == ConjTrans then  i   else  i*l
                else if h == ConjTrans then -i*l else -i
        c  = if h == ConjTrans then Conj else NoConj
        p' = p `advancePtr` o
        k  = diagLen (m,n) i
        s  = l+1
    in IOVector c k f p' s
{-# INLINE unsafeDiagViewIOMatrix #-}

maybeViewVectorAsRowIOMatrix :: IOVector p e -> Maybe (IOMatrix p1 e)
maybeViewVectorAsRowIOMatrix (IOVector c n f p s)
    | c == Conj && s == 1 =
        Just $ IOMatrix ConjTrans 1 n f p (max 1 n)
    | c == NoConj =
        Just $ IOMatrix NoTrans   1 n f p s
    | otherwise =
        Nothing
{-# INLINE maybeViewVectorAsRowIOMatrix #-}

maybeViewVectorAsColIOMatrix :: IOVector n e -> Maybe (IOMatrix n1 e)
maybeViewVectorAsColIOMatrix (IOVector c n f p s)
    | c == Conj =
        Just $ IOMatrix ConjTrans n 1 f p s        
    | s == 1 =
        Just $ IOMatrix NoTrans   n 1 f p (max 1 n)
    | otherwise =
        Nothing
{-# INLINE maybeViewVectorAsColIOMatrix #-}

maybeViewIOMatrixAsVector :: IOMatrix np e -> Maybe (IOVector k e)
maybeViewIOMatrixAsVector (IOMatrix h m n f p l)
    | h == ConjTrans = Nothing
    | l /= m         = Nothing
    | otherwise      = Just $ IOVector NoConj (m*n) f p 1
{-# INLINE maybeViewIOMatrixAsVector #-}

maybeViewVectorAsIOMatrix :: (Int,Int) -> IOVector k e -> Maybe (IOMatrix np e)
maybeViewVectorAsIOMatrix (m,n) (IOVector c k f p inc)
    | m*n /= k =
        error $ "maybeViewVectorAsMatrix " ++ show (m,n)
              ++ " <vector of dim " ++ show k ++ ">: vector dimension"
              ++ " must equal product of specified dimensions"
    | c == Conj = Nothing
    | inc /= 1  = Nothing
    | otherwise = Just $ IOMatrix NoTrans m n f p m
{-# INLINE maybeViewVectorAsIOMatrix #-}

liftIOMatrix :: (forall n. IOVector n e -> IO ()) -> IOMatrix np e -> IO ()
liftIOMatrix g (IOMatrix h m n f p l)
    | h == ConjTrans && l == n = do
        g (IOVector Conj   (m*n) f p 1)
    | h == NoTrans   && l == m = do
        g (IOVector NoConj (m*n) f p 1)
    | otherwise =
        let (c,m',n') = if h == ConjTrans then (Conj,n,m) else (NoConj,m,n)
            end     = p `advancePtr` (n'*l)
            go p' | p' == end = return ()
                  | otherwise = do
                      g (IOVector c m' f p' 1)
                      go (p' `advancePtr` l)
        in go p

liftIOMatrix2 :: (forall k. IOVector k e -> IOVector k f -> IO ())
              -> IOMatrix np e -> IOMatrix np f -> IO ()
liftIOMatrix2 f a b =
    if isHermIOMatrix a == isHermIOMatrix b
        then case (maybeViewIOMatrixAsVector a, maybeViewIOMatrixAsVector b) of
                 ((Just x), (Just y)) -> f x y
                 _                    -> elementwise
        else elementwise
  where
    
    elementwise =
        let vecsA = if isHermIOMatrix a then rowViews
                                        else colViews
            vecsB = if isHermIOMatrix a then rowViews
                                        else colViews
            go (x:xs) (y:ys) = do f x y
                                  go xs ys
            go []     []     = return ()
            go _      _      = error $ printf 
                ("liftMatrix2 <matrix of shape %s> <matrix of shape %s>:"
                ++ " shape mismatch") (show $ shape a) (show $ shape b)
        in go (vecsA a) (vecsB b)
        
    rowViews c = [ unsafeRowViewIOMatrix c i | i <- [ 0..numRows c - 1] ]
    colViews c = [ unsafeColViewIOMatrix c j | j <- [ 0..numCols c - 1] ]
{-# INLINE liftIOMatrix2 #-}

-- | Perform an 'IO' action with a pointer to the first element of the
-- matrix.
withIOMatrix :: IOMatrix (n,p) e -> (Ptr e -> IO a) -> IO a
withIOMatrix (IOMatrix _ _ _ f p _) g = do
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
        return $ IOMatrix NoTrans m n f (unsafeForeignPtrToPtr f) (max 1 m)

newIOMatrix :: (Elem e) 
            => (Int,Int) -> [((Int,Int), e)] -> IO (IOMatrix (n,p) e)
newIOMatrix (m,n) ies =
    let go p (((i,j),e):ies') = do
            when (i < 0 || i >= m || j < 0 || j >= n) $ fail $
                error $ ( printf 
                    "newMatrix (%d,%d) [ ..., ((%d,%d),_), ... ]: invalid index"
                    m n i j )
            pokeElemOff p (i+j*m) e
            go p ies'
        go _ [] = return ()
    in do
        a  <- newZeroIOMatrix (m,n)
        withIOMatrix a $ \p -> go p ies
        return a
    
unsafeNewIOMatrix :: (Elem e) 
                  => (Int,Int) -> [((Int,Int), e)] -> IO (IOMatrix (n,p) e)
unsafeNewIOMatrix (m,n) ies =
    let go p (((i,j),e):ies') = do
            pokeElemOff p (i+j*m) e
            go p ies'
        go _ [] = return ()
    in do
        a  <- newZeroIOMatrix (m,n)
        withIOMatrix a $ \p -> go p ies
        return a

newListIOMatrix :: (Elem e) => (Int,Int) -> [e] -> IO (IOMatrix (n,p) e)
newListIOMatrix (m,n) es = do
    a  <- newZeroIOMatrix (m,n)
    withIOMatrix a $ flip pokeArray (take (m*n) es)
    return a

newIdentityIOMatrix :: (Elem e) => (Int,Int) -> IO (IOMatrix np e)
newIdentityIOMatrix mn = do
    a  <- newIOMatrix_ mn
    setIdentityIOMatrix a
    return a

setIdentityIOMatrix :: IOMatrix np e -> IO ()
setIdentityIOMatrix a@(IOMatrix _ _ _ _ _ _) = do
    setZeroIOMatrix a
    setConstantIOVector 1 (unsafeDiagViewIOMatrix a 0)

newColsIOMatrix :: (ReadVector x IO, Elem e)
                => (Int,Int) -> [x n e] -> IO (IOMatrix (n,p) e)
newColsIOMatrix (m,n) cs = 
    let go _ j (_:_)   | j == n = 
            error $ printf 
                "newColsMatrix (%d,%d) <list of length at least %d>: too many columns"
                m n (n+1)
        go _ j []      | j < n =
            error $ printf
                "newColsMatrix (%d,%d) <list of length %d>: too few columns"
                m n j
        go _ _ []      = return ()
        go a j (c:cs') = do
             unsafeCopyVector (unsafeColViewIOMatrix a j) c
             go a (j+1) cs'
    in do
        a  <- newZeroIOMatrix (m,n)
        go a 0 cs
        return a

newRowsIOMatrix :: (ReadVector x IO, Elem e) 
                => (Int,Int) -> [x p e] -> IO (IOMatrix (n,p) e)
newRowsIOMatrix (m,n) rs = 
    let go _ i (_:_)   | i == m = 
            error $ printf 
                "newRowsMatrix (%d,%d) <list of length at least %d>: too many rows"
                m n (m+1)
        go _ i []      | i < m =
            error $ printf
                "newRowsMatrix (%d,%d) <list of length %d>: too few rows"
                m n i
        go _ _ []      = return ()
        go a i (r:rs') = do
             unsafeCopyVector (unsafeRowViewIOMatrix a i) r
             go a (i+1) rs'
    in do
        a  <- newZeroIOMatrix (m,n)
        go a 0 rs
        return a

newColIOMatrix :: (ReadVector x IO, Elem e)
               => x n e -> IO (IOMatrix (n,one) e)
newColIOMatrix x = newColsIOMatrix (dim x,1) [x]

newRowIOMatrix :: (ReadVector x IO, Elem e)
               => x p e -> IO (IOMatrix (one,p) e)
newRowIOMatrix x = newRowsIOMatrix (1,dim x) [x]

newZeroIOMatrix :: (Elem e) => (Int,Int) -> IO (IOMatrix (n,p) e)
newZeroIOMatrix mn = do
    a  <- newIOMatrix_ mn
    setZeroIOMatrix a
    return a

newConstantIOMatrix :: (Elem e) => (Int,Int) -> e -> IO (IOMatrix (n,p) e)
newConstantIOMatrix mn e = do
    a  <- newIOMatrix_ mn
    setConstantIOMatrix e a
    return a

newCopyIOMatrix :: IOMatrix np e -> IO (IOMatrix np e)
newCopyIOMatrix a@(IOMatrix h m n _ _ _) =
    let (m',n') = if h == ConjTrans then (n,m) else (m,n)
    in do
        b <- newIOMatrix_ (m',n')
        let b' = if h == ConjTrans then hermIOMatrix b else b
        unsafeCopyIOMatrix b' a
        return b'
{-# INLINE newCopyIOMatrix #-}

newCopyIOMatrix' :: IOMatrix np e -> IO (IOMatrix np e)
newCopyIOMatrix' a@(IOMatrix _ m n _ _ _) = do
    b  <- newIOMatrix_ (m,n)
    unsafeCopyIOMatrix b a
    return b
{-# INLINE newCopyIOMatrix' #-}

unsafeCopyIOMatrix :: IOMatrix np e -> IOMatrix np e -> IO ()
unsafeCopyIOMatrix = liftIOMatrix2 unsafeCopyIOVector
{-# INLINE unsafeCopyIOMatrix #-}

shapeIOMatrix :: IOMatrix np e -> (Int,Int)
shapeIOMatrix (IOMatrix _ m n _ _ _) = (m,n)
{-# INLINE shapeIOMatrix #-}

boundsIOMatrix :: IOMatrix np e -> ((Int,Int), (Int,Int))
boundsIOMatrix a = ((0,0), (m-1,n-1)) where (m,n) = shapeIOMatrix a
{-# INLINE boundsIOMatrix #-}

sizeIOMatrix :: IOMatrix np e -> Int
sizeIOMatrix (IOMatrix _ m n _ _ _) = m*n
{-# INLINE sizeIOMatrix #-}

getSizeIOMatrix :: IOMatrix np e -> IO Int
getSizeIOMatrix = return . sizeIOMatrix
{-# INLINE getSizeIOMatrix #-}

getMaxSizeIOMatrix :: IOMatrix np e -> IO Int
getMaxSizeIOMatrix = getSizeIOMatrix
{-# INLINE getMaxSizeIOMatrix #-}

indicesIOMatrix :: IOMatrix np e -> [(Int,Int)]
indicesIOMatrix (IOMatrix h m n _ _ _)
    | h == ConjTrans = [ (i,j) | i <- [ 0..m-1 ], j <- [ 0..n-1 ] ]
    | otherwise      = [ (i,j) | j <- [ 0..n-1 ], i <- [ 0..m-1 ] ]
{-# INLINE indicesIOMatrix #-}

getIndicesIOMatrix :: IOMatrix np e -> IO [(Int,Int)]
getIndicesIOMatrix = return . indicesIOMatrix
{-# INLINE getIndicesIOMatrix #-}

getIndicesIOMatrix' :: IOMatrix np e -> IO [(Int,Int)]
getIndicesIOMatrix' = getIndicesIOMatrix
{-# INLINE getIndicesIOMatrix' #-}

getElemsIOMatrix :: IOMatrix np e -> IO [e]
getElemsIOMatrix (IOMatrix h m n f p l)
    | h == ConjTrans = 
        liftM (map conjugate) $ 
            getElemsIOMatrix (IOMatrix NoTrans n m f p l)
    | l == m = 
        getElemsIOVector (IOVector NoConj (m*n) f p 1)
    | otherwise =
        let end   = p `advancePtr` (n*l)
            go p' | p' == end = return []
                  | otherwise = unsafeInterleaveIO $ do
                        c  <- getElemsIOVector (IOVector NoConj m f p' 1)
                        cs <- go (p' `advancePtr` l)
                        return (c ++ cs)
        in go p
{-# SPECIALIZE INLINE getElemsIOMatrix :: IOMatrix np Double -> IO [Double] #-}
{-# SPECIALIZE INLINE getElemsIOMatrix :: IOMatrix np (Complex Double) -> IO [Complex Double] #-}

getElemsIOMatrix' :: IOMatrix np e -> IO [e]
getElemsIOMatrix' (IOMatrix h m n f p l) 
    | h == ConjTrans = 
        liftM (map conjugate) $ 
            getElemsIOMatrix' (IOMatrix NoTrans n m f p l)
    | l == m    = 
        getElemsIOVector' (IOVector NoConj (m*n) f p 1)
    | otherwise =
        let end   = p `advancePtr` (n*l)
            go p' | p' == end = return []
                  | otherwise = do
                        c  <- getElemsIOVector' (IOVector NoConj m f p' 1)
                        cs <- go (p' `advancePtr` l)
                        return (c ++ cs)
        in go p
{-# SPECIALIZE INLINE getElemsIOMatrix' :: IOMatrix np Double -> IO [Double] #-}
{-# SPECIALIZE INLINE getElemsIOMatrix' :: IOMatrix np (Complex Double) -> IO [Complex Double] #-}

getAssocsIOMatrix :: IOMatrix np e -> IO [((Int,Int),e)]
getAssocsIOMatrix a = do
    is <- getIndicesIOMatrix a
    es <- getElemsIOMatrix a
    return $ zip is es
{-# INLINE getAssocsIOMatrix #-}

getAssocsIOMatrix' :: IOMatrix np e -> IO [((Int,Int),e)]
getAssocsIOMatrix' a = do
    is <- getIndicesIOMatrix' a
    es <- getElemsIOMatrix' a
    return $ zip is es
{-# INLINE getAssocsIOMatrix' #-}

unsafeReadElemIOMatrix :: IOMatrix np e -> (Int,Int) -> IO e
unsafeReadElemIOMatrix (IOMatrix h _ _ f p l) (i,j)
    | h == ConjTrans = do
        e <- peekElemOff p (i*l+j)
        touchForeignPtr f
        return $! conjugate e
    | otherwise = do
        e <- peekElemOff p (i+j*l)
        touchForeignPtr f
        return e
{-# SPECIALIZE INLINE unsafeReadElemIOMatrix :: IOMatrix n Double -> (Int,Int) -> IO (Double) #-}
{-# SPECIALIZE INLINE unsafeReadElemIOMatrix :: IOMatrix n (Complex Double) -> (Int,Int) -> IO (Complex Double) #-}        

canModifyElemIOMatrix :: IOMatrix np e -> (Int,Int) -> IO Bool
canModifyElemIOMatrix _ _ = return True
{-# INLINE canModifyElemIOMatrix #-}

unsafeWriteElemIOMatrix :: IOMatrix np e -> (Int,Int) -> e -> IO ()
unsafeWriteElemIOMatrix (IOMatrix h _ _ f p l) (i,j) e
    | h == ConjTrans = do
        io <- pokeElemOff p (i*l+j) (conjugate e)
        touchForeignPtr f
        return io
    | otherwise = do
        io <- pokeElemOff p (i+j*l) e
        touchForeignPtr f
        return io
{-# SPECIALIZE INLINE unsafeWriteElemIOMatrix :: IOMatrix n Double -> (Int,Int) -> Double -> IO () #-}
{-# SPECIALIZE INLINE unsafeWriteElemIOMatrix :: IOMatrix n (Complex Double) -> (Int,Int) -> Complex Double -> IO () #-}        

unsafeModifyElemIOMatrix :: IOMatrix n e -> (Int,Int) -> (e -> e) -> IO ()
unsafeModifyElemIOMatrix (IOMatrix h _ _ f p l) (i,j) g =
    let g' = if h == ConjTrans then conjugate . g . conjugate else g
        p' = if h == ConjTrans then p `advancePtr` (i*l+j) 
                               else p `advancePtr` (i+j*l)
    in do
        e  <- peek p'
        io <- poke p' (g' e)
        touchForeignPtr f
        return io
{-# SPECIALIZE INLINE unsafeModifyElemIOMatrix :: IOMatrix n Double -> (Int,Int) -> (Double -> Double) -> IO () #-}
{-# SPECIALIZE INLINE unsafeModifyElemIOMatrix :: IOMatrix n (Complex Double) -> (Int,Int) -> (Complex Double -> Complex Double) -> IO () #-}

unsafeSwapElemsIOMatrix :: IOMatrix n e -> (Int,Int) -> (Int,Int) -> IO ()
unsafeSwapElemsIOMatrix (IOMatrix h _ _ f p l) (i1,j1) (i2,j2) =
    let (p1,p2) = 
            if h == ConjTrans 
                then (p `advancePtr` (i1*l+j1), p `advancePtr` (i2*l+j2))
                else (p `advancePtr` (i1+j1*l), p `advancePtr` (i2+j2*l))
    in do
        e1  <- peek p1
        e2  <- peek p2
        poke p2 e1
        poke p1 e2
        touchForeignPtr f
{-# SPECIALIZE INLINE unsafeSwapElemsIOMatrix :: IOMatrix n Double -> (Int,Int) -> (Int,Int) -> IO () #-}
{-# SPECIALIZE INLINE unsafeSwapElemsIOMatrix :: IOMatrix n (Complex Double) -> (Int,Int) -> (Int,Int) -> IO () #-}

modifyWithIOMatrix :: (e -> e) -> IOMatrix np e -> IO ()
modifyWithIOMatrix g = liftIOMatrix (modifyWithIOVector g)
{-# INLINE modifyWithIOMatrix #-}

setZeroIOMatrix :: IOMatrix np e -> IO ()
setZeroIOMatrix = liftIOMatrix setZeroIOVector

setConstantIOMatrix :: e -> IOMatrix np e -> IO ()
setConstantIOMatrix k = liftIOMatrix (setConstantIOVector k)

getConjIOMatrix :: (BLAS1 e) => IOMatrix np e -> IO (IOMatrix np e)
getConjIOMatrix a = do
    b  <- newCopyIOMatrix a 
    doConjIOMatrix b
    return b

doConjIOMatrix :: (BLAS1 e) => IOMatrix np e -> IO ()
doConjIOMatrix = liftIOMatrix doConjIOVector

getScaledIOMatrix :: (BLAS1 e) => e -> IOMatrix np e -> IO (IOMatrix np e)
getScaledIOMatrix e a = do
    b  <- newCopyIOMatrix a
    scaleByIOMatrix e b
    return b

scaleByIOMatrix :: (BLAS1 e) => e -> IOMatrix np e -> IO ()
scaleByIOMatrix k = liftIOMatrix (scaleByIOVector k)

getShiftedIOMatrix :: (BLAS1 e) => e -> IOMatrix np e -> IO (IOMatrix np e)
getShiftedIOMatrix e a = do
    b  <- newCopyIOMatrix a
    shiftByIOMatrix e b
    return b
                    
shiftByIOMatrix :: (BLAS1 e) => e -> IOMatrix np e -> IO ()                    
shiftByIOMatrix k = liftIOMatrix (shiftByIOVector k)


instance Shaped IOMatrix (Int,Int) where
    shape = shapeIOMatrix
    {-# INLINE shape #-}
    bounds = boundsIOMatrix
    {-# INLINE bounds #-}

instance ReadTensor IOMatrix (Int,Int) IO where
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

instance WriteTensor IOMatrix (Int,Int) IO where
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

instance MatrixShaped IOMatrix where
    herm = hermIOMatrix
    {-# INLINE herm #-}
