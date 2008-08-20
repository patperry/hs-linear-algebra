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

    -- * converting to and from vectors
    maybeFromRow,
    maybeFromCol,
    maybeToVector,
    liftMatrix,
    liftMatrix2,

    module BLAS.Tensor,
    module BLAS.Numeric,
    module BLAS.Matrix,
    module Data.Matrix.Dense.Class,
    
    -- * Internal operations
    unsafeDoMatrixOp2

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
--import BLAS.Internal( inlinePerformIO, checkedRow, checkedCol, checkedDiag,
--    checkedSubmatrix, diagStart, diagLen )

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

-- | Create a matrix view of a row vector.  This will fail if the
-- vector is conjugated and the stride is not @1@.
maybeFromRow :: (Elem e, BaseMatrix a e, BaseVector x e) => 
    x m e -> Maybe (a (one,m) e)
maybeFromRow x
    | c && s == 1 =
        Just $ matrixViewArray f o (n,1) (max 1 n) True
    | not c =
        Just $ matrixViewArray f o (1,n) s         False
    | otherwise =
        Nothing
  where
    (f,o,n,s,c) = arrayFromVector x

-- | Possibly create a matrix view of a column vector.  This will fail
-- if the stride of the vector is not @1@ and the vector is not conjugated.
maybeFromCol :: (Elem e, BaseMatrix a e, BaseVector x e) => 
    x n e -> Maybe (a (n,one) e)
maybeFromCol x
    | c = maybeFromRow (conj x) >>= return . herm
    | s == 1 =
        Just $ matrixViewArray f o (n,1) (max 1 n) False
    | otherwise =
        Nothing
  where
    (f,o,n,s,c) = arrayFromVector x

maybeToVector :: (Elem e, BaseMatrix a e, RowColView a x e) => 
    a mn e -> Maybe (x k e)
maybeToVector a
    | h = 
        maybeToVector a' >>= return . conj
    | ld == m =
        Just $ vectorViewArray f o (m*n) 1  False
    | m == 1 =
        Just $ vectorViewArray f o n     ld False
    | otherwise =
        Nothing
  where
    a' = (coerceMatrix . herm . coerceMatrix) a
    (f,o,(m,n),ld,h) = arrayFromMatrix a

-- | Take a unary elementwise vector operation and apply it to the elements
-- of a matrix.
liftMatrix :: (Elem e, Monad m, BaseMatrix a e, RowColView a x e) => 
    (x k e -> m ()) -> a mn e -> m ()
liftMatrix f a =
    case maybeToVector a of
        Just x -> f x
        _ -> 
            let xs  = case isHerm a of
                          True  -> rowViews (coerceMatrix a)
                          False -> colViews (coerceMatrix a)
            in mapM_ f xs

-- | Take a binary elementwise vector operation and apply it to the elements
-- of a pair of matrices.
liftMatrix2 :: (Elem e, Monad m, BaseMatrix a e, BaseMatrix b e,
                    RowColView a x e, RowColView b y e) => 
    (x k e -> y k e -> m ()) ->
        a mn e -> b mn e -> m ()
liftMatrix2 f a b =
    if isHerm a == isHerm b
        then case (maybeToVector a, maybeToVector b) of
                 ((Just x), (Just y)) -> f x y
                 _                    -> elementwise
        else elementwise
  where
    elementwise =             
        let vecsA = if isHerm a then rowViews . coerceMatrix
                                else colViews . coerceMatrix
            vecsB = if isHerm a then rowViews . coerceMatrix
                                else colViews . coerceMatrix
            xs = vecsA a
            ys = vecsB b
        in zipWithM_ f xs ys


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

newCopyMatrix :: (BLAS1 e, ReadMatrix a x e m, WriteMatrix b y e m) => 
    a mn e -> m (b mn e)
newCopyMatrix a 
    | isHerm a =
        newCopyMatrix ((herm . coerceMatrix) a) >>= 
            return . coerceMatrix . herm
    | otherwise = do
        a' <- newMatrix_ (shape a)
        unsafeCopyMatrix a' a
        return a'

unsafeCopyMatrix :: (BLAS1 e, ReadMatrix a x e m, WriteMatrix b y e m) => 
    b mn e -> a mn e -> m ()
unsafeCopyMatrix = liftMatrix2 unsafeCopyVector

unsafeSwapIOMatrix :: (Elem e) => IOMatrix mn e -> IOMatrix mn e -> IO ()
unsafeSwapIOMatrix = liftMatrix2 unsafeSwap

{-
-- | @diag a 0@ gets a vector view of the main diagonal of @a@.  @diag a k@ for 
-- @k@ positive gets a view of the @k@th superdiagonal.  For @k@ negative, it
-- gets a view of the @(-k)@th subdiagonal.
diag :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t k e
diag a = checkedDiag (shape a) (unsafeDiag a)

-- | Same as 'diag', but does not do any bounds checking.
unsafeDiag :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t k e
unsafeDiag a i 
    | isHerm a = 
        conj $ unsafeDiag (herm a) (negate i)
    | otherwise =            
        let f = storageOf a
            o = indexOf a (diagStart i)
            n = diagLen (shape a) i
            s = ldaOf a + 1
            c = False
        in V.fromForeignPtr f o n s c





-}

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

instance (Elem e, WriteVector x e IO) => RowColRead IOMatrix x e IO where
    unsafeGetRow = undefined
    unsafeGetCol = undefined

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

unsafeDoMatrixOp2 :: (BLAS1 e, ReadMatrix a x e m, ReadMatrix b y e m, WriteMatrix c z e m) =>
    (c n e -> b n e -> m ()) -> a n e -> b n e -> c n e -> m ()
unsafeDoMatrixOp2 f a b c = do
    unsafeCopyMatrix c a
    f c b
