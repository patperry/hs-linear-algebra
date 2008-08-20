{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses, UndecidableInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Internal (
    -- * IOMatrix matrix data types
    IOMatrix,

    module BLAS.Tensor,
    module BLAS.Numeric,
    module BLAS.Matrix,
    module Data.Matrix.Dense.Class,
    ) where

import Control.Monad
import Data.Ix( range )
import Foreign
import System.IO.Unsafe( unsafeInterleaveIO )
import Unsafe.Coerce

import BLAS.Tensor  hiding ( ITensor(..) )
import BLAS.Numeric hiding ( INumeric(..) )
import BLAS.Matrix  hiding ( BaseMatrix, IRowCol(..) )
import qualified BLAS.Matrix as BLAS

import BLAS.Elem( Elem, BLAS1 )
import BLAS.Internal( diagStart, diagLen )

import Data.Vector.Dense.Internal

import Data.Matrix.Dense.Class


-- | The mutable dense matrix data type.  It can either store elements in 
-- column-major order, or provide a view into another matrix.  The view 
-- transposes and conjugates the underlying matrix.
data IOMatrix mn e =
      DM { storageIOMatrix :: {-# UNPACK #-} !(ForeignPtr e) -- ^ a pointer to the storage region
         , offsetIOMatrix  :: {-# UNPACK #-} !Int            -- ^ an offset (in elements, not bytes) to the first element in the matrix. 
         , numRowsIOMatrix :: {-# UNPACK #-} !Int            -- ^ the number of rows in the matrix
         , numColsIOMatrix :: {-# UNPACK #-} !Int            -- ^ the number of columns in the matrix
         , ldaIOMatrix     :: {-# UNPACK #-} !Int            -- ^ the leading dimension size of the matrix
         , isHermIOMatrix  :: {-# UNPACK #-} !Bool           -- ^ indicates whether or not the matrix is transposed and conjugated
         }

coerceIOMatrix :: IOMatrix mn e -> IOMatrix mn' e
coerceIOMatrix = unsafeCoerce


-- | Create a new matrix of given shape, but do not initialize the elements.
newIOMatrix_ :: (Elem e) => (Int,Int) -> IO (IOMatrix mn e)
newIOMatrix_ (m,n) 
    | m < 0 || n < 0 =
        ioError $ userError $ 
            "Tried to create a matrix with shape `" ++ show (m,n) ++ "'"
    | otherwise = do
        f <- mallocForeignPtrArray (m*n)
        return $ DM f 0 m n (max 1 m) False

hermIOMatrix :: IOMatrix (m,n) e -> IOMatrix (n,m) e
hermIOMatrix a = a{ isHermIOMatrix=(not . isHermIOMatrix) a }

unsafeSubmatrixIOMatrix :: 
    IOMatrix mn e -> (Int,Int) -> (Int,Int) -> IOMatrix mn' e
unsafeSubmatrixIOMatrix a (i,j) (m,n)
    | isHermIOMatrix a  = 
        coerceIOMatrix $ hermIOMatrix $ 
            unsafeSubmatrixIOMatrix (hermIOMatrix $ coerceIOMatrix a) (j,i) (n,m)
    | otherwise =
        let f = storageIOMatrix a
            o = indexOfIOMatrix a (i,j)
            l = ldaIOMatrix a
        in DM f o m n l False

unsafeRowViewIOMatrix :: (Elem e) => IOMatrix (m,n) e -> Int -> IOVector n e
unsafeRowViewIOMatrix a i
    | isHermIOMatrix a =
        conj $ unsafeColViewIOMatrix (herm a) i
    | otherwise =
        let f = storageIOMatrix a
            o = indexOfIOMatrix a (i,0)
            n = numColsIOMatrix a
            s = ldaIOMatrix a
            c = False
        in vectorViewArray f o n s c

unsafeColViewIOMatrix :: (Elem e) => IOMatrix (m,n) e -> Int -> IOVector m e
unsafeColViewIOMatrix a j 
    | isHermIOMatrix a =
        conj $ unsafeRowViewIOMatrix (hermIOMatrix a) j
    | otherwise =
        let f = storageIOMatrix a
            o = indexOfIOMatrix a (0,j)
            m = numRowsIOMatrix a
            s = 1
            c = False
        in vectorViewArray f o m s c

unsafeDiagViewIOMatrix :: (Elem e) => IOMatrix (m,n) e -> Int -> IOVector k e
unsafeDiagViewIOMatrix a i 
    | isHermIOMatrix a = 
        conj $ unsafeDiagViewIOMatrix (hermIOMatrix a) (negate i)
    | otherwise =            
        let f = storageIOMatrix a
            o = indexOfIOMatrix a (diagStart i)
            n = diagLen (shapeIOMatrix a) i
            s = ldaIOMatrix a + 1
            c = False
        in vectorViewArray f o n s c

colViewsIOMatrix :: (Elem e) => IOMatrix (m,n) e -> [IOVector m e]
colViewsIOMatrix a = [ unsafeColViewIOMatrix a i | i <- [0..(numCols a - 1)] ]

    
withIOMatrixPtr :: (Elem e) => IOMatrix mn e -> (Ptr e -> IO a) -> IO a
withIOMatrixPtr a f =
    withForeignPtr (storageIOMatrix a) $ \ptr ->
        let ptr' = ptr `advancePtr` (offsetIOMatrix a)
        in f ptr'


indexOfIOMatrix :: IOMatrix mn e -> (Int,Int) -> Int
indexOfIOMatrix a (i,j) = 
    let (i',j') = case isHermIOMatrix a of
                        True  -> (j,i)
                        False -> (i,j)
        o = offsetIOMatrix a
        l = ldaIOMatrix a
    in o + i' + j'*l
{-# INLINE indexOfIOMatrix #-}
    
    
shapeIOMatrix :: IOMatrix mn e -> (Int,Int)
shapeIOMatrix a = (numRowsIOMatrix a, numColsIOMatrix a)

boundsIOMatrix :: IOMatrix mn e -> ((Int,Int), (Int,Int))
boundsIOMatrix a = case shapeIOMatrix a of (m,n) -> ((0,0),(m-1,n-1))

getSizeIOMatrix :: IOMatrix mn e -> IO Int
getSizeIOMatrix a = return (m*n) where (m,n) = shapeIOMatrix a

indicesIOMatrix :: IOMatrix mn e -> [(Int,Int)]
indicesIOMatrix a 
    | isHermIOMatrix a = [ (i,j) | i <- range (0,m-1), j <- range (0,n-1) ]
    | otherwise        = [ (i,j) | j <- range (0,n-1), i <- range (0,m-1) ]
  where (m,n) = shapeIOMatrix a

getIndicesIOMatrix :: IOMatrix mn e -> IO [(Int,Int)]
getIndicesIOMatrix = return . indicesIOMatrix
{-# INLINE getIndicesIOMatrix #-}

getElemsIOMatrix :: (Elem e) => IOMatrix mn e -> IO [e]
getElemsIOMatrix a
    | isHermIOMatrix a = getElemsIOMatrix (hermIOMatrix $ coerceIOMatrix a) >>= 
                             return . map conj
    | otherwise = 
        liftM concat $
            unsafeInterleaveIO $ 
                mapM getElems (colViewsIOMatrix $ coerceIOMatrix a)

getAssocsIOMatrix :: (Elem e) => IOMatrix mn e -> IO [((Int,Int),e)]
getAssocsIOMatrix a = do
    is <- getIndicesIOMatrix a
    es <- getElemsIOMatrix a
    return $ zip is es
    
getIndicesIOMatrix' :: IOMatrix mn e -> IO [(Int,Int)]
getIndicesIOMatrix' = getIndicesIOMatrix
{-# INLINE getIndicesIOMatrix' #-}

getElemsIOMatrix' :: (Elem e) => IOMatrix mn e -> IO [e]
getElemsIOMatrix' a
    | isHermIOMatrix a = getElemsIOMatrix' (hermIOMatrix $ coerceIOMatrix a) >>= 
                             return . map conj
    | otherwise = 
        liftM concat $
            mapM getElems' (colViewsIOMatrix $ coerceIOMatrix a)

getAssocsIOMatrix' :: (Elem e) => IOMatrix mn e -> IO [((Int,Int),e)]
getAssocsIOMatrix' a = do
    is <- getIndicesIOMatrix' a
    es <- getElemsIOMatrix' a
    return $ zip is es

unsafeReadElemIOMatrix :: (Elem e) => IOMatrix mn e -> (Int,Int) -> IO e
unsafeReadElemIOMatrix a (i,j)
    | isHerm a  = unsafeReadElem (hermIOMatrix $ coerceIOMatrix a) (j,i) >>= 
                      return . conj
    | otherwise = withForeignPtr (storageIOMatrix a) $ \ptr ->
                      peekElemOff ptr (indexOfIOMatrix a (i,j))
{-# INLINE unsafeReadElemIOMatrix #-}

newZeroIOMatrix :: (Elem e) => (Int,Int) -> IO (IOMatrix mn e)
newZeroIOMatrix mn = do
    a <- newMatrix_ mn
    setZero a
    return a

newConstantIOMatrix :: (Elem e) => (Int,Int) -> e -> IO (IOMatrix mn e)
newConstantIOMatrix mn e = do
    a <- newMatrix_ mn
    setConstant e a
    return a

setZeroIOMatrix :: (Elem e) => IOMatrix mn e -> IO ()    
setZeroIOMatrix = liftMatrix setZero

setConstantIOMatrix :: (Elem e) => e -> IOMatrix mn e -> IO ()
setConstantIOMatrix e = liftMatrix (setConstant e)

unsafeWriteElemIOMatrix :: (Elem e) => 
    IOMatrix mn e -> (Int,Int) -> e -> IO ()
unsafeWriteElemIOMatrix a (i,j) e
    | isHermIOMatrix a = unsafeWriteElem a' (j,i) $ conj e
    | otherwise        = withForeignPtr (storageIOMatrix a) $ \ptr ->
                             pokeElemOff ptr (indexOfIOMatrix a (i,j)) e
  where
    a' = (herm . coerceIOMatrix) a

modifyWithIOMatrix :: (Elem e) => (e -> e) -> IOMatrix mn e -> IO ()
modifyWithIOMatrix f = liftMatrix (modifyWith f)

canModifyElemIOMatrix :: IOMatrix mn e -> (Int,Int) -> IO Bool
canModifyElemIOMatrix _ _ = return True
{-# INLINE canModifyElemIOMatrix #-}

unsafeSwapIOMatrix :: (Elem e) => IOMatrix mn e -> IOMatrix mn e -> IO ()
unsafeSwapIOMatrix = liftMatrix2 unsafeSwap


instance BaseTensor IOMatrix (Int,Int) e where
    shape  = shapeIOMatrix
    bounds = boundsIOMatrix

instance (Elem e) => ReadTensor IOMatrix (Int,Int) e IO where
    getSize        = getSizeIOMatrix

    getAssocs      = getAssocsIOMatrix
    getIndices     = getIndicesIOMatrix
    getElems       = getElemsIOMatrix

    getAssocs'     = getAssocsIOMatrix'
    getIndices'    = getIndicesIOMatrix'
    getElems'      = getElemsIOMatrix'
    
    unsafeReadElem = unsafeReadElemIOMatrix

instance (Elem e) => WriteTensor IOMatrix (Int,Int) e IO where
    newZero         = newZeroIOMatrix
    newConstant     = newConstantIOMatrix

    setConstant     = setConstantIOMatrix
    setZero         = setZeroIOMatrix

    modifyWith      = modifyWithIOMatrix
    unsafeWriteElem = unsafeWriteElemIOMatrix
    canModifyElem   = canModifyElemIOMatrix
    
    unsafeSwap      = unsafeSwapIOMatrix

instance (Elem e) => BLAS.BaseMatrix IOMatrix e where
    herm            = hermIOMatrix
    
instance (Elem e) => BaseMatrix IOMatrix e where
    lda             = ldaIOMatrix
    isHerm          = isHermIOMatrix
    unsafeSubmatrix = unsafeSubmatrixIOMatrix
    withMatrixPtr   = withIOMatrixPtr
    matrixViewArray f o (m,n) = DM f o m n
    arrayFromMatrix (DM f o m n l h) = (f,o,(m,n),l,h)

instance (Elem e) => ReadMatrix IOMatrix IOVector e IO where
    
instance (Elem e) => WriteMatrix IOMatrix IOVector e IO where
    newMatrix_ = newIOMatrix_

instance (Elem e) => RowColView IOMatrix IOVector e where
    unsafeRowView = unsafeRowViewIOMatrix
    unsafeColView = unsafeColViewIOMatrix

instance (BLAS1 e) => RowColRead IOMatrix IOVector e IO where
    unsafeGetRow a i = newCopyVector (unsafeRowView a i)
    unsafeGetCol a j = newCopyVector (unsafeColView a j)

instance (Elem e) => DiagView IOMatrix IOVector e where
    unsafeDiagView = unsafeDiagViewIOMatrix

instance (BLAS1 e) => DiagRead IOMatrix IOVector e IO where
    unsafeGetDiag a i = newCopyVector (unsafeDiagView a i)

instance (BLAS1 e) => ReadNumeric IOMatrix (Int,Int) e IO where

instance (BLAS1 e) => WriteNumeric IOMatrix (Int,Int) e IO where
    doConj    = liftMatrix doConj
    scaleBy k = liftMatrix (scaleBy k)
    shiftBy k = liftMatrix (shiftBy k)

instance (BLAS1 e, ReadMatrix a x e IO) => CopyTensor a IOMatrix (Int,Int) e IO where
    newCopyTensor    = newCopyMatrix
    unsafeCopyTensor = unsafeCopyMatrix

instance (BLAS1 e, ReadMatrix a x e IO) => Numeric2 a IOMatrix (Int,Int) e IO where
    unsafeAxpy k = liftMatrix2 (unsafeAxpy k)
    unsafeMul    = liftMatrix2 unsafeMul
    unsafeDiv    = liftMatrix2 unsafeDiv

instance (BLAS1 e, ReadMatrix a x e IO, ReadMatrix b y e IO) => Numeric3 a b IOMatrix (Int,Int) e IO where
    unsafeDoAdd = unsafeDoMatrixOp2 $ flip $ unsafeAxpy 1
    unsafeDoSub = unsafeDoMatrixOp2 $ flip $ unsafeAxpy (-1)
    unsafeDoMul = unsafeDoMatrixOp2 $ unsafeMul
    unsafeDoDiv = unsafeDoMatrixOp2 $ unsafeDiv
