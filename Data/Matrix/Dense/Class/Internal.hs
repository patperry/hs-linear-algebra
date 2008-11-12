{-# LANGUAGE TypeFamilies, FlexibleInstances, MultiParamTypeClasses, FunctionalDependencies
    , FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Internal (

    -- * Matrix types
    IOMatrix,
    STMatrix,
    unsafeIOMatrixToSTMatrix,
    unsafeSTMatrixToIOMatrix,

    -- * Matrix type classes
    BaseMatrix_(..),
    BaseMatrix,
    ReadMatrix,
    WriteMatrix,
    
    -- * Basic matrix properties
    ldaOfMatrix,
    isHermMatrix,

    -- * Coercing the matrix shape
    coerceMatrix,

    -- * Converting to and from vectors
    maybeFromRow,
    maybeFromCol,
    maybeToVector,
    
    -- * Lifting vector operations
    liftMatrix,
    liftMatrix2,

    -- * BaseTensor functions
    shapeMatrix,
    boundsMatrix,
    
    -- * BaseMatrix functions
    hermMatrix,

    -- * ReadTensor functions
    getSizeMatrix,
    getAssocsMatrix,
    getIndicesMatrix,
    getElemsMatrix,
    getAssocsMatrix',
    getIndicesMatrix',
    getElemsMatrix',
    unsafeReadElemMatrix,

    -- * WriteTensor functions
    newMatrix_,
    newZeroMatrix,
    setZeroMatrix,
    newConstantMatrix,
    setConstantMatrix,
    modifyWithMatrix,
    canModifyElemMatrix,
    unsafeWriteElemMatrix,
    
    -- * CopyTensor functions
    newCopyMatrix,
    unsafeCopyMatrix,
    
    -- * SwapTensor functions
    unsafeSwapMatrix,

    -- * Vector views
    rowViews,
    colViews,
    unsafeRowView,
    unsafeColView,
    unsafeDiagView,
    unsafeGetRowMatrix,
    unsafeGetColMatrix,
    
    -- * Numeric functions
    doConjMatrix,
    scaleByMatrix,
    shiftByMatrix,
    
    -- * Numeric2 functions
    unsafeAxpyMatrix,
    unsafeMulMatrix,
    unsafeDivMatrix,
    
    -- * Numeric3 functions
    unsafeDoAddMatrix,
    unsafeDoSubMatrix,
    unsafeDoMulMatrix,
    unsafeDoDivMatrix,
    
    -- * ReadApply functions
    gemv,
    gemm,
    
    -- * Utility functions
    withMatrixPtr,
    indexOfMatrix,
    indicesMatrix,
    unsafeDoMatrixOp2,
    
    ) where

import Control.Monad
import Control.Monad.ST
import Data.Ix
import Foreign
import Unsafe.Coerce

import BLAS.Elem
import BLAS.C.Types
import qualified BLAS.C.Level2 as BLAS
import qualified BLAS.C.Level3 as BLAS
import BLAS.Internal( diagStart, diagLen )
import BLAS.UnsafeIOToM

import BLAS.Tensor

import Data.Vector.Dense.Class.Internal( IOVector, STVector,
    BaseVector(..), ReadVector, WriteVector, 
    newCopyVector, unsafeCopyVector, unsafeSwapVector, 
    doConjVector, scaleByVector, shiftByVector, unsafeAxpyVector, 
    unsafeMulVector, unsafeDivVector, withVectorPtr, dim, stride, isConj )

import BLAS.Matrix.Base hiding ( BaseMatrix )
import qualified BLAS.Matrix.Base as BLAS


class (BLAS.BaseMatrix a) => BaseMatrix_ a where
    type VectorView a :: * -> * -> *
    matrixViewArray :: ForeignPtr e -> Ptr e -> Int -> Int -> Int -> Bool -> a mn e
    arrayFromMatrix :: a mn e -> (ForeignPtr e, Ptr e, Int, Int, Int, Bool)

class (BaseMatrix_ a, BaseVector (VectorView a)) => BaseMatrix a

class (BaseMatrix a, UnsafeIOToM m, ReadTensor a (Int,Int) m
      , ReadVector (VectorView a) m) =>
    ReadMatrix a m where

class (ReadMatrix a m, WriteTensor a (Int,Int) m
      , WriteVector (VectorView a) m) =>
    WriteMatrix a m | m -> a where

------------------------- Basic Matrix Properties ---------------------------

size1 :: (BaseMatrix a) => a mn e -> Int
size1 a = let (_,_,m,_,_,_) = arrayFromMatrix a in m
{-# INLINE size1  #-}

size2 :: (BaseMatrix a) => a mn e -> Int
size2 a = let (_,_,_,n,_,_) = arrayFromMatrix a in n
{-# INLINE size2 #-}

ldaOfMatrix :: (BaseMatrix a) => a mn e -> Int
ldaOfMatrix a = let (_,_,_,_,l,_) = arrayFromMatrix a in l
{-# INLINE ldaOfMatrix #-}

isHermMatrix :: (BaseMatrix a) => a mn e -> Bool
isHermMatrix a = let (_,_,_,_,_,h) = arrayFromMatrix a in h
{-# INLINE isHermMatrix #-}

-- | Cast the shape type of the matrix.
coerceMatrix :: (BaseMatrix a) => a mn e -> a mn' e
coerceMatrix = unsafeCoerce
{-# INLINE coerceMatrix #-}

----------------------- Converting to/from Vectors --------------------------

-- | Create a matrix view of a row vector.  This will fail if the
-- vector is conjugated and the stride is not @1@.
maybeFromRow :: (BaseMatrix a) => 
    VectorView a n e -> Maybe (a (one,n) e)
maybeFromRow x
    | c && s == 1 =
        Just $ matrixViewArray f p n 1 (max 1 n) True
    | not c =
        Just $ matrixViewArray f p 1 n s         False
    | otherwise =
        Nothing
  where
    (f,p,n,s,c) = arrayFromVector x

-- | Possibly create a matrix view of a column vector.  This will fail
-- if the stride of the vector is not @1@ and the vector is not conjugated.
maybeFromCol :: (BaseMatrix a) => 
    VectorView a n e -> Maybe (a (n,one) e)
maybeFromCol x
    | c =
        Just $ herm $ matrixViewArray f p 1 n s False
    | s == 1 =
        Just $ matrixViewArray f p n 1 (max 1 n) False
    | otherwise =
        Nothing
  where
    (f,p,n,s,c) = arrayFromVector x

maybeToVector :: (BaseMatrix a) => 
    a mn e -> Maybe (VectorView a k e)
maybeToVector a
    | h = 
        maybeToVector a' >>= return . conj
    | ld == m =
        Just $ vectorViewArray f p (m*n) 1  False
    | m == 1 =
        Just $ vectorViewArray f p n     ld False
    | otherwise =
        Nothing
  where
    a' = (coerceMatrix . herm . coerceMatrix) a
    (f,p,m,n,ld,h) = arrayFromMatrix a


----------------------- Lifting vector operations ---------------------------

-- | Take a unary elementwise vector operation and apply it to the elements
-- of a matrix.
liftMatrix :: (ReadMatrix a m, Storable e) =>
    (VectorView a k e -> m ()) -> a mn e -> m ()
liftMatrix f a =
    case maybeToVector a of
        Just x -> f x
        _ -> 
            let xs  = case isHermMatrix a of
                          True  -> rowViews (coerceMatrix a)
                          False -> colViews (coerceMatrix a)
            in mapM_ f xs

-- | Take a binary elementwise vector operation and apply it to the elements
-- of a pair of matrices.
liftMatrix2 :: (Monad m, BaseMatrix a, BaseMatrix b, Storable e) =>
    (VectorView a k e -> VectorView b k e -> m ()) ->
        a mn e -> b mn e -> m ()
liftMatrix2 f a b =
    if isHermMatrix a == isHermMatrix b
        then case (maybeToVector a, maybeToVector b) of
                 ((Just x), (Just y)) -> f x y
                 _                    -> elementwise
        else elementwise
  where
    elementwise =             
        let vecsA = if isHermMatrix a then rowViews . coerceMatrix
                                      else colViews . coerceMatrix
            vecsB = if isHermMatrix a then rowViews . coerceMatrix
                                      else colViews . coerceMatrix
            xs = vecsA a
            ys = vecsB b
        in zipWithM_ f xs ys


-------------------------- BaseTensor functions -----------------------------

shapeMatrix :: (BaseMatrix a) => a mn e -> (Int,Int)
shapeMatrix a | isHermMatrix a  = (size2 a, size1 a)
              | otherwise       = (size1 a, size2 a)
{-# INLINE shapeMatrix #-}

boundsMatrix :: (BaseMatrix a) => a mn e -> ((Int,Int), (Int,Int))
boundsMatrix a = ((0,0), (m-1,n-1)) where (m,n) = shapeMatrix a
{-# INLINE boundsMatrix #-}


-------------------------- BaseMatrix functions -----------------------------

hermMatrix :: (BaseMatrix a) => a (m,n) e -> a (n,m) e
hermMatrix a = let (f,p,m,n,l,h) = arrayFromMatrix a
               in matrixViewArray f p m n l (not h)
{-# INLINE hermMatrix #-}


-------------------------- ReadTensor functions -----------------------------

getSizeMatrix :: (ReadMatrix a m) => a mn e -> m Int
getSizeMatrix a = return (m*n) where (m,n) = shape a

getIndicesMatrix :: (ReadMatrix a m) => a mn e -> m [(Int,Int)]
getIndicesMatrix = return . indicesMatrix
{-# INLINE getIndicesMatrix #-}

getElemsMatrix :: (ReadMatrix a m, Elem e) => a mn e -> m [e]
getElemsMatrix a
    | isHermMatrix a = getElemsMatrix (herm $ coerceMatrix a) >>= 
                           return . map conj
    | otherwise = 
        liftM concat $
            unsafeInterleaveM $ 
                mapM getElems (colViews $ coerceMatrix a)

getAssocsMatrix :: (ReadMatrix a m, Elem e) => a mn e -> m [((Int,Int),e)]
getAssocsMatrix a = do
    is <- getIndicesMatrix a
    es <- getElemsMatrix a
    return $ zip is es
    
getIndicesMatrix' :: (ReadMatrix a m) => a mn e -> m [(Int,Int)]
getIndicesMatrix' = getIndicesMatrix
{-# INLINE getIndicesMatrix' #-}

getElemsMatrix' :: (ReadMatrix a m, Elem e) => a mn e -> m [e]
getElemsMatrix' a
    | isHermMatrix a = getElemsMatrix' (herm $ coerceMatrix a) >>= 
                           return . map conj
    | otherwise = 
        liftM concat $
            mapM getElems' (colViews $ coerceMatrix a)

getAssocsMatrix' :: (ReadMatrix a m, Elem e) => a mn e -> m [((Int,Int),e)]
getAssocsMatrix' a = do
    is <- getIndicesMatrix' a
    es <- getElemsMatrix' a
    return $ zip is es

unsafeReadElemMatrix :: (ReadMatrix a m, Elem e) => a mn e -> (Int,Int) -> m e
unsafeReadElemMatrix a (i,j)
    | isHermMatrix a = unsafeReadElem (herm $ coerceMatrix a) (j,i) >>= 
                           return . conj
    | otherwise = unsafeIOToM $
                      withMatrixPtr a $ \ptr ->
                          peekElemOff ptr (indexOfMatrix a (i,j))
{-# INLINE unsafeReadElemMatrix #-}


------------------------- WriteTensor functions -----------------------------

-- | Create a new matrix of given shape, but do not initialize the elements.
newMatrix_ :: (WriteMatrix a m, Elem e) => (Int,Int) -> m (a mn e)
newMatrix_ (m,n) 
    | m < 0 || n < 0 =
        fail $ 
            "Tried to create a matrix with shape `" ++ show (m,n) ++ "'"
    | otherwise = unsafeIOToM $ do
        f <- mallocForeignPtrArray (m*n)
        return $ matrixViewArray f (unsafeForeignPtrToPtr f) m n (max 1 m) False

-- | Create a zero matrix of the specified shape.
newZeroMatrix :: (WriteMatrix a m, Elem e) => (Int,Int) -> m (a mn e)
newZeroMatrix mn = do
    a <- newMatrix_ mn
    setZero a
    return a

-- | Create a constant matrix of the specified shape.
newConstantMatrix :: (WriteMatrix a m, Elem e) => (Int,Int) -> e -> m (a mn e)
newConstantMatrix mn e = do
    a <- newMatrix_ mn
    setConstant e a
    return a

setZeroMatrix :: (WriteMatrix a m, Elem e) => a mn e -> m ()    
setZeroMatrix = liftMatrix setZero

setConstantMatrix :: (WriteMatrix a m, Elem e) => e -> a mn e -> m ()
setConstantMatrix e = liftMatrix (setConstant e)

unsafeWriteElemMatrix :: (WriteMatrix a m, Elem e) => 
    a mn e -> (Int,Int) -> e -> m ()
unsafeWriteElemMatrix a (i,j) e
    | isHermMatrix a  = unsafeWriteElem a' (j,i) $ conj e
    | otherwise       = unsafeIOToM $
                            withMatrixPtr a $ \ptr ->
                                pokeElemOff ptr (indexOfMatrix a (i,j)) e
  where
    a' = (herm . coerceMatrix) a

modifyWithMatrix :: (WriteMatrix a m, Elem e) => (e -> e) -> a mn e -> m ()
modifyWithMatrix f = liftMatrix (modifyWith f)

canModifyElemMatrix :: (WriteMatrix a m) => a mn e -> (Int,Int) -> m Bool
canModifyElemMatrix _ _ = return True
{-# INLINE canModifyElemMatrix #-}


------------------------- CopyTensor functions ------------------------------

newCopyMatrix :: (BLAS1 e, ReadMatrix a m, WriteMatrix b m) => 
    a mn e -> m (b mn e)
newCopyMatrix a 
    | isHermMatrix a =
        newCopyMatrix ((herm . coerceMatrix) a) >>= 
            return . coerceMatrix . herm
    | otherwise = do
        a' <- newMatrix_ (shape a)
        unsafeCopyMatrix a' a
        return a'

unsafeCopyMatrix :: (BLAS1 e, WriteMatrix b m,  ReadMatrix a m) => 
    b mn e -> a mn e -> m ()
unsafeCopyMatrix = liftMatrix2 unsafeCopyVector


------------------------- SwapTensor functions ------------------------------

unsafeSwapMatrix :: (WriteMatrix a m, BLAS1 e) => a mn e -> a mn e -> m ()
unsafeSwapMatrix = liftMatrix2 unsafeSwapVector


------------------------------ Vector views ---------------------------------

unsafeRowView :: (BaseMatrix a, Storable e) => 
    a (k,l) e -> Int -> VectorView a l e
unsafeRowView a i
    | isHermMatrix a =
        conj $ unsafeColView (herm a) i
    | otherwise =
        let (fp,p,_,n,ld,_) = arrayFromMatrix a
            o  = indexOfMatrix a (i,0)
            p' = p `advancePtr` o
            s  = ld
            c  = False
        in vectorViewArray fp p' n s c

unsafeColView :: (BaseMatrix a, Storable e) => 
    a (k,l) e -> Int -> VectorView a k e
unsafeColView a j 
    | isHermMatrix a =
        conj $ unsafeRowView (herm a) j
    | otherwise =
        let (fp,p,m,_,_,_) = arrayFromMatrix a
            o  = indexOfMatrix a (0,j)
            p' = p `advancePtr` o
            s  = 1
            c  = False
        in vectorViewArray fp p' m s c

unsafeDiagView :: (BaseMatrix a, Storable e) => 
    a mn e -> Int -> VectorView a k e
unsafeDiagView a i 
    | isHermMatrix a = 
        conj $ unsafeDiagView (herm $ coerceMatrix a) (negate i)
    | otherwise =            
        let (fp,p,m,n,ld,_) = arrayFromMatrix a
            o  = indexOfMatrix a (diagStart i)
            p' = p `advancePtr` o
            n' = diagLen (m,n) i
            s  = ld + 1
            c  = False
        in vectorViewArray fp p' n' s c

-- | Get a list of vector views of the rows of the matrix.
rowViews :: (BaseMatrix a, Storable e) => a (m,n) e -> [VectorView a n e]
rowViews a = [ unsafeRowView a i | i <- [0..numRows a - 1] ]

-- | Get a list of vector views of the columns of the matrix.
colViews :: (BaseMatrix a, Storable e) => a (m,n) e -> [VectorView a m e]
colViews a = [ unsafeColView a j | j <- [0..numCols a - 1] ]

-- | Same as 'getRow' but not range-checked.
unsafeGetRowMatrix :: (ReadMatrix a m, WriteVector y m, BLAS1 e) => 
    a (k,l) e -> Int -> m (y l e)    
unsafeGetRowMatrix a i = newCopyVector (unsafeRowView a i)

-- | Same as 'getCol' but not range-checked.
unsafeGetColMatrix :: (ReadMatrix a m, WriteVector y m, BLAS1 e) => 
    a (k,l) e -> Int -> m (y k e)
unsafeGetColMatrix a j = newCopyVector (unsafeColView a j)


--------------------------- Numeric functions -------------------------------

doConjMatrix :: (WriteMatrix a m, BLAS1 e) => a mn e -> m ()
doConjMatrix = liftMatrix doConjVector

scaleByMatrix :: (WriteMatrix a m, BLAS1 e) => e -> a mn e -> m ()
scaleByMatrix k = liftMatrix (scaleByVector k)

shiftByMatrix :: (WriteMatrix a m, BLAS1 e) => e -> a mn e -> m ()
shiftByMatrix k = liftMatrix (shiftByVector k)


-------------------------- Numeric2 functions -------------------------------

unsafeAxpyMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    e -> a mn e -> b mn e -> m ()
unsafeAxpyMatrix alpha = liftMatrix2 (unsafeAxpyVector alpha)

unsafeMulMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b mn e -> a mn e -> m ()
unsafeMulMatrix = liftMatrix2 unsafeMulVector

unsafeDivMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b mn e -> a mn e -> m ()
unsafeDivMatrix = liftMatrix2 unsafeDivVector


-------------------------- Numeric3 functions -------------------------------

unsafeDoAddMatrix :: (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a mn e -> b mn e -> c mn e -> m ()
unsafeDoAddMatrix = unsafeDoMatrixOp2 $ flip $ unsafeAxpyMatrix 1

unsafeDoSubMatrix :: (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a mn e -> b mn e -> c mn e -> m ()
unsafeDoSubMatrix = unsafeDoMatrixOp2 $ flip $ unsafeAxpyMatrix (-1)

unsafeDoMulMatrix :: (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a mn e -> b mn e -> c mn e -> m ()
unsafeDoMulMatrix = unsafeDoMatrixOp2 $ unsafeMulMatrix

unsafeDoDivMatrix :: (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a mn e -> b mn e -> c mn e -> m ()
unsafeDoDivMatrix = unsafeDoMatrixOp2 $ unsafeDivMatrix


-------------------------- ReadApply functions -------------------------------

-- | @gemv alpha a x beta y@ replaces @y := alpha a * x + beta y@.
gemv :: (ReadMatrix a m, ReadVector x m, WriteVector y m, BLAS3 e) => 
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
gemv alpha a x beta y
    | numRows a == 0 || numCols a == 0 =
        scaleBy beta y
        
    | isConj y && (isConj x || stride x == 1) =
        let order  = colMajor
            transA = if isConj x then noTrans else conjTrans
            transB = blasTransOf (herm a)
            m      = 1
            n      = dim y
            k      = dim x
            ldA    = stride x
            ldB    = ldaOfMatrix a
            ldC    = stride y
            alpha' = conj alpha
            beta'  = conj beta
        in unsafeIOToM $
               withVectorPtr x $ \pA ->
               withMatrixPtr a $ \pB ->
               withVectorPtr y $ \pC ->
                   BLAS.gemm order transA transB m n k alpha' pA ldA pB ldB beta' pC ldC
    
    | (isConj y && otherwise) || isConj x = do
        doConj y
        gemv alpha a x beta (conj y)
        doConj y
        
    | otherwise =
        let order  = colMajor
            transA = blasTransOf a
            (m,n)  = case (isHermMatrix a) of
                         False -> shape a
                         True  -> (flipShape . shape) a
            ldA    = ldaOfMatrix a
            incX   = stride x
            incY   = stride y
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withVectorPtr x $ \pX ->
               withVectorPtr y $ \pY -> do
                   BLAS.gemv order transA m n alpha pA ldA pX incX beta pY incY

-- | @gemm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gemm :: (BLAS3 e, ReadMatrix a m, ReadMatrix b m, WriteMatrix c m) => 
    e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
gemm alpha a b beta c
    | numRows a == 0 || numCols a == 0 || numCols b == 0 = 
        scaleBy beta c
    | isHermMatrix c = gemm (conj alpha) (herm b) (herm a) (conj beta) (herm c)
    | otherwise =
        let order  = colMajor
            transA = blasTransOf a
            transB = blasTransOf b
            (m,n)  = shape c
            k      = numCols a
            ldA    = ldaOfMatrix a
            ldB    = ldaOfMatrix b
            ldC    = ldaOfMatrix c
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withMatrixPtr b $ \pB ->
               withMatrixPtr c $ \pC ->
                   BLAS.gemm order transA transB m n k alpha pA ldA pB ldB beta pC ldC


--------------------------- Utility functions -------------------------------

blasTransOf :: (BaseMatrix a) => a mn e -> CBLASTrans
blasTransOf a = 
    case (isHermMatrix a) of
          False -> noTrans
          True  -> conjTrans

flipShape :: (Int,Int) -> (Int,Int)
flipShape (m,n) = (n,m)

withMatrixPtr :: (BaseMatrix a) =>
    a mn e -> (Ptr e -> IO b) -> IO b
withMatrixPtr a f =
    let (fp,p,_,_,_,_) = arrayFromMatrix a
    in do
        b <- f p
        touchForeignPtr fp
        return b

indexOfMatrix :: (BaseMatrix a) => a mn e -> (Int,Int) -> Int
indexOfMatrix a (i,j) = 
    let (i',j') = case isHermMatrix a of
                        True  -> (j,i)
                        False -> (i,j)
        l = ldaOfMatrix a
    in i' + j'*l
{-# INLINE indexOfMatrix #-}

indicesMatrix :: (BaseMatrix a) => a mn e -> [(Int,Int)]
indicesMatrix a 
    | isHermMatrix a = [ (i,j) | i <- range (0,m-1), j <- range (0,n-1) ]
    | otherwise      = [ (i,j) | j <- range (0,n-1), i <- range (0,m-1) ]
  where (m,n) = shape a

unsafeDoMatrixOp2 :: (BLAS1 e, ReadMatrix a m, ReadMatrix b m, WriteMatrix c m) =>
    (c n e -> b n e -> m ()) -> a n e -> b n e -> c n e -> m ()
unsafeDoMatrixOp2 f a b c = do
    unsafeCopyMatrix c a
    f c b


------------------------------------ Instances ------------------------------

-- | The mutable dense matrix data type.  It can either store elements in 
-- column-major order, or provide a view into another matrix.  The view 
-- transposes and conjugates the underlying matrix.
data IOMatrix mn e =
      DM {-# UNPACK #-} !(ForeignPtr e) -- a pointer to the storage region
         {-# UNPACK #-} !(Ptr e)        -- a pointer to the first element
         {-# UNPACK #-} !Int            -- the number of rows in the matrix
         {-# UNPACK #-} !Int            -- the number of columns in the matrix
         {-# UNPACK #-} !Int            -- the leading dimension size of the matrix
         {-# UNPACK #-} !Bool           -- indicates whether or not the matrix is transposed and conjugated

newtype STMatrix s n e = ST (IOMatrix n e)

unsafeIOMatrixToSTMatrix :: IOMatrix n e -> STMatrix s n e
unsafeIOMatrixToSTMatrix = ST
{-# INLINE unsafeIOMatrixToSTMatrix #-}

unsafeSTMatrixToIOMatrix :: STMatrix s n e -> IOMatrix n e
unsafeSTMatrixToIOMatrix (ST x) = x
{-# INLINE unsafeSTMatrixToIOMatrix #-}


instance BaseMatrix_ IOMatrix where
    type VectorView IOMatrix = IOVector
    matrixViewArray f p m n = DM f p m n
    arrayFromMatrix (DM f p m n l h) = (f,p,m,n,l,h)

instance BaseMatrix_ (STMatrix s) where
    type VectorView (STMatrix s) = STVector s
    matrixViewArray f p m n l h = ST $ DM f p m n l h
    arrayFromMatrix (ST (DM f p m n l h)) = (f,p,m,n,l,h)

instance BaseMatrix IOMatrix
instance BaseMatrix (STMatrix s)

instance BaseTensor IOMatrix (Int,Int) where
    shape  = shapeMatrix
    bounds = boundsMatrix

instance BaseTensor (STMatrix s) (Int,Int) where
    shape  = shapeMatrix
    bounds = boundsMatrix

instance ReadTensor IOMatrix (Int,Int) IO where
    getSize        = getSizeMatrix
    getAssocs      = getAssocsMatrix
    getIndices     = getIndicesMatrix
    getElems       = getElemsMatrix
    getAssocs'     = getAssocsMatrix'
    getIndices'    = getIndicesMatrix'
    getElems'      = getElemsMatrix'
    unsafeReadElem = unsafeReadElemMatrix

instance ReadTensor (STMatrix s) (Int,Int) (ST s) where
    getSize        = getSizeMatrix
    getAssocs      = getAssocsMatrix
    getIndices     = getIndicesMatrix
    getElems       = getElemsMatrix
    getAssocs'     = getAssocsMatrix'
    getIndices'    = getIndicesMatrix'
    getElems'      = getElemsMatrix'
    unsafeReadElem = unsafeReadElemMatrix

instance WriteTensor IOMatrix (Int,Int) IO where
    setConstant     = setConstantMatrix
    setZero         = setZeroMatrix
    modifyWith      = modifyWithMatrix
    unsafeWriteElem = unsafeWriteElemMatrix
    canModifyElem   = canModifyElemMatrix
    doConj          = doConjMatrix
    scaleBy         = scaleByMatrix
    shiftBy         = shiftByMatrix


instance WriteTensor (STMatrix s) (Int,Int) (ST s) where
    setConstant     = setConstantMatrix
    setZero         = setZeroMatrix
    modifyWith      = modifyWithMatrix
    unsafeWriteElem = unsafeWriteElemMatrix
    canModifyElem   = canModifyElemMatrix
    doConj          = doConjMatrix
    scaleBy         = scaleByMatrix
    shiftBy         = shiftByMatrix

instance BLAS.BaseMatrix IOMatrix where
    herm = hermMatrix

instance BLAS.BaseMatrix (STMatrix s) where
    herm = hermMatrix
    
instance ReadMatrix IOMatrix IO where
instance ReadMatrix (STMatrix s) (ST s) where    
    
instance WriteMatrix IOMatrix IO where
instance WriteMatrix (STMatrix s) (ST s) where
