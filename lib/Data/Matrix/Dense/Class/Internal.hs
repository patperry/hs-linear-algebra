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

    -- * Matrix views
    unsafeSubmatrixView,

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
    
    -- * MMatrix functions
    MMatrix(..),
    gemv,
    gemm,
    hemv,
    hemm,
    hemv',
    hemm',
    unsafeDoSApplyAddTriMatrix,
    unsafeDoSApplyAddMatTriMatrix,
    trmv,
    trmm,
    unsafeDoSSolveTriMatrix,
    unsafeDoSSolveMatTriMatrix,
    trsv,
    trsm,
    
    -- * MSolve functions
    MSolve(..),
    
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
import BLAS.Types

import Data.Matrix.Herm
import Data.Matrix.Tri.Internal
import Data.Vector.Dense.Class.Internal( IOVector, STVector,
    BaseVector(..), ReadVector, WriteVector, newVector_,
    newCopyVector, unsafeCopyVector, unsafeSwapVector, 
    doConjVector, scaleByVector, shiftByVector, unsafeAxpyVector, 
    unsafeMulVector, unsafeDivVector, withVectorPtr, dim, stride, isConj,
    coerceVector )
import Data.Vector.Dense.Class.Views( unsafeSubvectorView )
import Data.Vector.Dense.Class.Special( newBasisVector )

import BLAS.Matrix.Shaped
import BLAS.Matrix.HasVectorView


class (Storable e, MatrixShaped a e, HasVectorView a) => BaseMatrix_ a e where
    matrixViewArray :: ForeignPtr e -> Ptr e -> Int -> Int -> Int -> Bool -> a mn e
    arrayFromMatrix :: a mn e -> (ForeignPtr e, Ptr e, Int, Int, Int, Bool)

class (BaseMatrix_ a e, BaseVector (VectorView a) e) => BaseMatrix a e

class (BaseMatrix a e, UnsafeIOToM m, BLAS3 e, ReadTensor a (Int,Int) e m
      , MMatrix a e m, MMatrix (Herm a) e m, MMatrix (Tri a) e m
      , MSolve (Tri a) e m
      , ReadVector (VectorView a) e m) =>
    ReadMatrix a e m where

class (ReadMatrix a e m, WriteTensor a (Int,Int) e m
      , WriteVector (VectorView a) e m) =>
    WriteMatrix a e m | m -> a where

------------------------- Basic Matrix Properties ---------------------------

size1 :: (BaseMatrix a e) => a mn e -> Int
size1 a = let (_,_,m,_,_,_) = arrayFromMatrix a in m
{-# INLINE size1  #-}

size2 :: (BaseMatrix a e) => a mn e -> Int
size2 a = let (_,_,_,n,_,_) = arrayFromMatrix a in n
{-# INLINE size2 #-}

ldaOfMatrix :: (BaseMatrix a e) => a mn e -> Int
ldaOfMatrix a = let (_,_,_,_,l,_) = arrayFromMatrix a in l
{-# INLINE ldaOfMatrix #-}

isHermMatrix :: (BaseMatrix a e) => a mn e -> Bool
isHermMatrix a = let (_,_,_,_,_,h) = arrayFromMatrix a in h
{-# INLINE isHermMatrix #-}

-- | Cast the shape type of the matrix.
coerceMatrix :: (BaseMatrix a e) => a mn e -> a mn' e
coerceMatrix = unsafeCoerce
{-# INLINE coerceMatrix #-}

----------------------- Converting to/from Vectors --------------------------

-- | Create a matrix view of a row vector.  This will fail if the
-- vector is conjugated and the stride is not @1@.
maybeFromRow :: (BaseMatrix a e) => 
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
maybeFromCol :: (BaseMatrix a e) => 
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

maybeToVector :: (BaseMatrix a e) => 
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
liftMatrix :: (ReadMatrix a e m) =>
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
liftMatrix2 :: (Monad m, BaseMatrix a e, BaseMatrix b f) =>
    (VectorView a k e -> VectorView b k f -> m ()) ->
        a mn e -> b mn f -> m ()
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

shapeMatrix :: (BaseMatrix a e) => a mn e -> (Int,Int)
shapeMatrix a | isHermMatrix a  = (size2 a, size1 a)
              | otherwise       = (size1 a, size2 a)
{-# INLINE shapeMatrix #-}

boundsMatrix :: (BaseMatrix a e) => a mn e -> ((Int,Int), (Int,Int))
boundsMatrix a = ((0,0), (m-1,n-1)) where (m,n) = shapeMatrix a
{-# INLINE boundsMatrix #-}


-------------------------- BaseMatrix functions -----------------------------

hermMatrix :: (BaseMatrix a e) => a (m,n) e -> a (n,m) e
hermMatrix a = let (f,p,m,n,l,h) = arrayFromMatrix a
               in matrixViewArray f p m n l (not h)
{-# INLINE hermMatrix #-}


-------------------------- ReadTensor functions -----------------------------

getSizeMatrix :: (ReadMatrix a e m) => a mn e -> m Int
getSizeMatrix a = return (m*n) where (m,n) = shape a

getIndicesMatrix :: (ReadMatrix a e m) => a mn e -> m [(Int,Int)]
getIndicesMatrix = return . indicesMatrix
{-# INLINE getIndicesMatrix #-}

getElemsMatrix :: (ReadMatrix a e m) => a mn e -> m [e]
getElemsMatrix a
    | isHermMatrix a = getElemsMatrix (herm $ coerceMatrix a) >>= 
                           return . map conj
    | otherwise = 
        liftM concat $
            unsafeInterleaveM $ 
                mapM getElems (colViews $ coerceMatrix a)

getAssocsMatrix :: (ReadMatrix a e m) => a mn e -> m [((Int,Int),e)]
getAssocsMatrix a = do
    is <- getIndicesMatrix a
    es <- getElemsMatrix a
    return $ zip is es
    
getIndicesMatrix' :: (ReadMatrix a e m) => a mn e -> m [(Int,Int)]
getIndicesMatrix' = getIndicesMatrix
{-# INLINE getIndicesMatrix' #-}

getElemsMatrix' :: (ReadMatrix a e m) => a mn e -> m [e]
getElemsMatrix' a
    | isHermMatrix a = getElemsMatrix' (herm $ coerceMatrix a) >>= 
                           return . map conj
    | otherwise = 
        liftM concat $
            mapM getElems' (colViews $ coerceMatrix a)

getAssocsMatrix' :: (ReadMatrix a e m) => a mn e -> m [((Int,Int),e)]
getAssocsMatrix' a = do
    is <- getIndicesMatrix' a
    es <- getElemsMatrix' a
    return $ zip is es

unsafeReadElemMatrix :: (ReadMatrix a e m) => a mn e -> (Int,Int) -> m e
unsafeReadElemMatrix a (i,j)
    | isHermMatrix a = unsafeReadElem (herm $ coerceMatrix a) (j,i) >>= 
                           return . conj
    | otherwise = unsafeIOToM $
                      withMatrixPtr a $ \ptr ->
                          peekElemOff ptr (indexOfMatrix a (i,j))
{-# INLINE unsafeReadElemMatrix #-}


------------------------- WriteTensor functions -----------------------------

-- | Create a new matrix of given shape, but do not initialize the elements.
newMatrix_ :: (WriteMatrix a e m) => (Int,Int) -> m (a mn e)
newMatrix_ (m,n) 
    | m < 0 || n < 0 =
        fail $ 
            "Tried to create a matrix with shape `" ++ show (m,n) ++ "'"
    | otherwise = unsafeIOToM $ do
        f <- mallocForeignPtrArray (m*n)
        return $ matrixViewArray f (unsafeForeignPtrToPtr f) m n (max 1 m) False

-- | Create a zero matrix of the specified shape.
newZeroMatrix :: (WriteMatrix a e m) => (Int,Int) -> m (a mn e)
newZeroMatrix mn = do
    a <- newMatrix_ mn
    setZero a
    return a

-- | Create a constant matrix of the specified shape.
newConstantMatrix :: (WriteMatrix a e m) => (Int,Int) -> e -> m (a mn e)
newConstantMatrix mn e = do
    a <- newMatrix_ mn
    setConstant e a
    return a

setZeroMatrix :: (WriteMatrix a e m) => a mn e -> m ()    
setZeroMatrix = liftMatrix setZero

setConstantMatrix :: (WriteMatrix a e m) => e -> a mn e -> m ()
setConstantMatrix e = liftMatrix (setConstant e)

unsafeWriteElemMatrix :: (WriteMatrix a e m) => 
    a mn e -> (Int,Int) -> e -> m ()
unsafeWriteElemMatrix a (i,j) e
    | isHermMatrix a  = unsafeWriteElem a' (j,i) $ conj e
    | otherwise       = unsafeIOToM $
                            withMatrixPtr a $ \ptr ->
                                pokeElemOff ptr (indexOfMatrix a (i,j)) e
  where
    a' = (herm . coerceMatrix) a

modifyWithMatrix :: (WriteMatrix a e m) => (e -> e) -> a mn e -> m ()
modifyWithMatrix f = liftMatrix (modifyWith f)

canModifyElemMatrix :: (WriteMatrix a e m) => a mn e -> (Int,Int) -> m Bool
canModifyElemMatrix _ _ = return True
{-# INLINE canModifyElemMatrix #-}


------------------------- CopyTensor functions ------------------------------

newCopyMatrix :: (ReadMatrix a e m, WriteMatrix b e m) => 
    a mn e -> m (b mn e)
newCopyMatrix a 
    | isHermMatrix a =
        newCopyMatrix ((herm . coerceMatrix) a) >>= 
            return . coerceMatrix . herm
    | otherwise = do
        a' <- newMatrix_ (shape a)
        unsafeCopyMatrix a' a
        return a'

unsafeCopyMatrix :: (WriteMatrix b e m, ReadMatrix a e m) => 
    b mn e -> a mn e -> m ()
unsafeCopyMatrix = liftMatrix2 unsafeCopyVector


------------------------- SwapTensor functions ------------------------------

unsafeSwapMatrix :: (WriteMatrix a e m) => a mn e -> a mn e -> m ()
unsafeSwapMatrix = liftMatrix2 unsafeSwapVector


------------------------------ Vector views ---------------------------------

-- | Same as 'submatrixView' but indices are not range-checked.
unsafeSubmatrixView :: (BaseMatrix a e) => 
    a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
unsafeSubmatrixView a (i,j) (m,n)
    | isHermMatrix a  = 
        coerceMatrix $ herm $ 
            unsafeSubmatrixView (herm $ coerceMatrix a) (j,i) (n,m)
    | otherwise =
        let (fp,p,_,_,ld,_) = arrayFromMatrix a
            o  = indexOfMatrix a (i,j)
            p' = p `advancePtr` o
        in matrixViewArray fp p' m n ld False

unsafeRowView :: (BaseMatrix a e) => 
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

unsafeColView :: (BaseMatrix a e) => 
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

unsafeDiagView :: (BaseMatrix a e) => 
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
rowViews :: (BaseMatrix a e) => a (m,n) e -> [VectorView a n e]
rowViews a = [ unsafeRowView a i | i <- [0..numRows a - 1] ]

-- | Get a list of vector views of the columns of the matrix.
colViews :: (BaseMatrix a e) => a (m,n) e -> [VectorView a m e]
colViews a = [ unsafeColView a j | j <- [0..numCols a - 1] ]

-- | Same as 'getRow' but not range-checked.
unsafeGetRowMatrix :: (ReadMatrix a e m, WriteVector y e m) => 
    a (k,l) e -> Int -> m (y l e)    
unsafeGetRowMatrix a i = newCopyVector (unsafeRowView a i)

-- | Same as 'getCol' but not range-checked.
unsafeGetColMatrix :: (ReadMatrix a e m, WriteVector y e m) => 
    a (k,l) e -> Int -> m (y k e)
unsafeGetColMatrix a j = newCopyVector (unsafeColView a j)


--------------------------- Numeric functions -------------------------------

doConjMatrix :: (WriteMatrix a e m) => a mn e -> m ()
doConjMatrix = liftMatrix doConjVector

scaleByMatrix :: (WriteMatrix a e m) => e -> a mn e -> m ()
scaleByMatrix k = liftMatrix (scaleByVector k)

shiftByMatrix :: (WriteMatrix a e m) => e -> a mn e -> m ()
shiftByMatrix k = liftMatrix (shiftByVector k)


-------------------------- Numeric2 functions -------------------------------

unsafeAxpyMatrix :: (ReadMatrix a e m, WriteMatrix b e m) =>
    e -> a mn e -> b mn e -> m ()
unsafeAxpyMatrix alpha = liftMatrix2 (unsafeAxpyVector alpha)

unsafeMulMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b mn e -> a mn e -> m ()
unsafeMulMatrix = liftMatrix2 unsafeMulVector

unsafeDivMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b mn e -> a mn e -> m ()
unsafeDivMatrix = liftMatrix2 unsafeDivVector


-------------------------- Numeric3 functions -------------------------------

unsafeDoAddMatrix :: (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    a mn e -> b mn e -> c mn e -> m ()
unsafeDoAddMatrix = unsafeDoMatrixOp2 $ flip $ unsafeAxpyMatrix 1

unsafeDoSubMatrix :: (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    a mn e -> b mn e -> c mn e -> m ()
unsafeDoSubMatrix = unsafeDoMatrixOp2 $ flip $ unsafeAxpyMatrix (-1)

unsafeDoMulMatrix :: (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    a mn e -> b mn e -> c mn e -> m ()
unsafeDoMulMatrix = unsafeDoMatrixOp2 $ unsafeMulMatrix

unsafeDoDivMatrix :: (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    a mn e -> b mn e -> c mn e -> m ()
unsafeDoDivMatrix = unsafeDoMatrixOp2 $ unsafeDivMatrix



--------------------------- Utility functions -------------------------------

blasTransOf :: (BaseMatrix a e) => a mn e -> CBLASTrans
blasTransOf a = 
    case (isHermMatrix a) of
          False -> noTrans
          True  -> conjTrans

flipShape :: (Int,Int) -> (Int,Int)
flipShape (m,n) = (n,m)

withMatrixPtr :: (BaseMatrix a e) =>
    a mn e -> (Ptr e -> IO b) -> IO b
withMatrixPtr a f =
    let (fp,p,_,_,_,_) = arrayFromMatrix a
    in do
        b <- f p
        touchForeignPtr fp
        return b

indexOfMatrix :: (BaseMatrix a e) => a mn e -> (Int,Int) -> Int
indexOfMatrix a (i,j) = 
    let (i',j') = case isHermMatrix a of
                        True  -> (j,i)
                        False -> (i,j)
        l = ldaOfMatrix a
    in i' + j'*l
{-# INLINE indexOfMatrix #-}

indicesMatrix :: (BaseMatrix a e) => a mn e -> [(Int,Int)]
indicesMatrix a 
    | isHermMatrix a = [ (i,j) | i <- range (0,m-1), j <- range (0,n-1) ]
    | otherwise      = [ (i,j) | j <- range (0,n-1), i <- range (0,m-1) ]
  where (m,n) = shape a

unsafeDoMatrixOp2 :: (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    (c n e -> b n e -> m ()) -> a n e -> b n e -> c n e -> m ()
unsafeDoMatrixOp2 f a b c = do
    unsafeCopyMatrix c a
    f c b


-------------------------------- MMatrix ----------------------------------

-- | Minimal complete definition: (unsafeDoSApplyAdd, unsafeDoSApplyAddMat)
class (MatrixShaped a e, Monad m) => MMatrix a e m where
    unsafeGetSApply :: (ReadVector x e m, WriteVector y e m) =>
        e -> a (k,l) e -> x l e -> m (y k e)
    unsafeGetSApply alpha a x = do
        y <- newVector_ (numRows a)
        unsafeDoSApplyAdd alpha a x 0 y
        return y

    unsafeGetSApplyMat :: (ReadMatrix b e m, WriteMatrix c e m) =>
        e -> a (r,s) e -> b (s,t) e -> m (c (r,t) e)
    unsafeGetSApplyMat alpha a b = do
        c <- newMatrix_ (numRows a, numCols b)
        unsafeDoSApplyAddMat alpha a b 0 c
        return c

    unsafeDoSApplyAdd :: (ReadVector x e m, WriteVector y e m) =>
        e -> a (k,l) e -> x l e -> e -> y k e -> m ()
    unsafeDoSApplyAdd alpha a x beta y = do
        y' <- unsafeGetSApply alpha a x
        scaleBy beta y
        unsafeAxpyVector 1 y' y

    unsafeDoSApplyAddMat :: (ReadMatrix b e m, WriteMatrix c e m) =>
        e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
    unsafeDoSApplyAddMat alpha a b beta c = do
        c' <- unsafeGetSApplyMat alpha a b
        scaleBy beta c
        unsafeAxpyMatrix 1 c' c

    unsafeDoSApply_ :: (WriteVector y e m) =>
        e -> a (n,n) e -> y n e -> m ()
    unsafeDoSApply_ alpha a x = do
        y <- newVector_ (dim x)
        unsafeDoSApplyAdd alpha a x 0 y
        unsafeCopyVector x y

    unsafeDoSApplyMat_ :: (WriteMatrix b e m) =>
        e -> a (k,k) e -> b (k,l) e -> m ()
    unsafeDoSApplyMat_ alpha a b = do
        c <- newMatrix_ (shape b)
        unsafeDoSApplyAddMat alpha a b 0 c
        unsafeCopyMatrix b c

    unsafeGetRow :: (WriteVector x e m) => a (k,l) e -> Int -> m (x l e)
    unsafeGetRow a i = do
        e <- newBasisVector (numRows a) i
        liftM conj $ unsafeGetSApply 1 (herm a) e
        
    unsafeGetCol :: (WriteVector x e m) => a (k,l) e -> Int -> m (x k e)
    unsafeGetCol a j = do
        e <- newBasisVector (numCols a) j
        unsafeGetSApply 1 a e


-- | @gemv alpha a x beta y@ replaces @y := alpha a * x + beta y@.
gemv :: (ReadMatrix a e m, ReadVector x e m, WriteVector y e m) => 
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
gemv alpha a x beta y
    | numRows a == 0 || numCols a == 0 =
        scaleBy beta y
        
    | isConj y && (isConj x || stride x == 1) =
        let transA = if isConj x then noTrans else conjTrans
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
                   BLAS.gemm transA transB m n k alpha' pA ldA pB ldB beta' pC ldC
    
    | (isConj y && otherwise) || isConj x = do
        doConj y
        gemv alpha a x beta (conj y)
        doConj y
        
    | otherwise =
        let transA = blasTransOf a
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
                   BLAS.gemv transA m n alpha pA ldA pX incX beta pY incY

-- | @gemm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gemm :: (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
gemm alpha a b beta c
    | numRows a == 0 || numCols a == 0 || numCols b == 0 = 
        scaleBy beta c
    | isHermMatrix c = gemm (conj alpha) (herm b) (herm a) (conj beta) (herm c)
    | otherwise =
        let transA = blasTransOf a
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
                   BLAS.gemm transA transB m n k alpha pA ldA pB ldB beta pC ldC

hemv :: (ReadMatrix a e m, ReadVector x e m, WriteVector y e m) => 
    e -> Herm a (k,k) e -> x k e -> e -> y k e -> m ()
hemv alpha h x beta y
    | numRows h == 0 =
        return ()
    | isConj y = do
        doConj y
        hemv alpha h x beta (conj y)
        doConj y
    | isConj x = do
        x' <- newCopyVector x
        doConj x'
        hemv alpha h (conj x') beta y
    | otherwise =
        let (u,a) = hermToBase h
            n     = numCols a
            u'    = case isHermMatrix a of
                        True  -> flipUpLo u
                        False -> u
            uploA = cblasUpLo u'
            ldA   = ldaOfMatrix a
            incX  = stride x
            incY  = stride y
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withVectorPtr x $ \pX ->
               withVectorPtr y $ \pY ->
                   BLAS.hemv uploA n alpha pA ldA pX incX beta pY incY

hemm :: (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> Herm a (k,k) e -> b (k,l) e -> e -> c (k,l) e -> m ()
hemm alpha h b beta c
    | numRows b == 0 || numCols b == 0 || numCols c == 0 = return ()
    | (isHermMatrix a) /= (isHermMatrix c) || (isHermMatrix a) /= (isHermMatrix b) =
        zipWithM_ (\x y -> hemv alpha h x beta y) (colViews b) (colViews c)
    | otherwise =
        let (m,n)   = shape c
            (side,u',m',n')
                    = if isHermMatrix a
                          then (rightSide, flipUpLo u, n, m)
                          else (leftSide,  u,          m, n)
            uploA   = cblasUpLo u'
            ldA     = ldaOfMatrix a
            ldB     = ldaOfMatrix b
            ldC     = ldaOfMatrix c
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withMatrixPtr b $ \pB ->
               withMatrixPtr c $ \pC ->
                   BLAS.hemm side uploA m' n' alpha pA ldA pB ldB beta pC ldC
    where
      (u,a) = hermToBase h

hemv' :: (ReadMatrix a e m, ReadVector x e m, WriteVector y e m) => 
    e -> Herm a (r,s) e -> x s e -> e -> y r e -> m ()
hemv' alpha a x beta y = 
    hemv alpha (coerceHerm a) x beta (coerceVector y)

hemm' :: (ReadMatrix a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> Herm a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
hemm' alpha a b beta c = 
    hemm alpha (coerceHerm a) b beta (coerceMatrix c)

unsafeDoSApplyAddTriMatrix :: (ReadMatrix a e m, MMatrix a e m, 
    ReadVector x e m, WriteVector y e m) =>
        e -> Tri a (k,l) e -> x l e -> e -> y k e -> m ()
unsafeDoSApplyAddTriMatrix alpha t x beta y =
    if beta == 0
        then unsafeDoSApplyTriMatrix alpha t x y
        else do
            y' <- newCopyVector y
            unsafeDoSApplyTriMatrix alpha t x y'
            scaleBy beta y
            unsafeAxpyVector 1 y' y

unsafeDoSApplyAddMatTriMatrix :: (ReadMatrix a e m, MMatrix a e m, 
    ReadMatrix b e m, WriteMatrix c e m) =>
        e -> Tri a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
unsafeDoSApplyAddMatTriMatrix alpha t b beta c =
    if beta == 0
        then unsafeDoSApplyMatTriMatrix alpha t b c
        else do
            c' <- newCopyMatrix c
            unsafeDoSApplyMatTriMatrix alpha t b c'
            scaleBy beta c
            unsafeAxpyMatrix 1 c' c

unsafeDoSApplyTriMatrix :: (ReadMatrix a e m, MMatrix a e m, 
    ReadVector x e m, WriteVector y e m) =>
        e -> Tri a (k,l) e -> x l e -> y k e -> m ()
unsafeDoSApplyTriMatrix alpha t x y =
    case (u, toLower d a, toUpper d a) of
        (Lower,Left t',_) -> do
            unsafeCopyVector y (coerceVector x)
            trmv alpha t' y
            
        (Lower,Right (t',r),_) -> do
            let y1 = unsafeSubvectorView y 0            (numRows t')
                y2 = unsafeSubvectorView y (numRows t') (numRows r)
            unsafeCopyVector y1 x
            trmv alpha t' y1
            unsafeDoSApplyAdd alpha r x 0 y2
            
        (Upper,_,Left t') -> do
            unsafeCopyVector (coerceVector y) x
            trmv alpha t' (coerceVector y)

        (Upper,_,Right (t',r)) ->
            let x1 = unsafeSubvectorView x 0            (numCols t')
                x2 = unsafeSubvectorView x (numCols t') (numCols r)
            in do
                unsafeCopyVector y x1
                trmv alpha t' y
                unsafeDoSApplyAdd alpha r x2 1 y
  where
    (u,d,a) = triToBase t

unsafeDoSApplyMatTriMatrix :: (ReadMatrix a e m, MMatrix a e m, 
    ReadMatrix b e m, WriteMatrix c e m) =>
        e -> Tri a (r,s) e -> b (s,t) e -> c (r,t) e -> m ()
unsafeDoSApplyMatTriMatrix alpha t b c =
    case (u, toLower d a, toUpper d a) of
        (Lower,Left t',_) -> do
            unsafeCopyMatrix c (coerceMatrix b)
            trmm alpha t' c
            
        (Lower,Right (t',r),_) -> do
            let c1 = unsafeSubmatrixView c (0,0)          (numRows t',numCols c)
                c2 = unsafeSubmatrixView c (numRows t',0) (numRows r ,numCols c)
            unsafeCopyMatrix c1 b
            trmm alpha t' c1
            unsafeDoSApplyAddMat alpha r b 0 c2
            
        (Upper,_,Left t') -> do
            unsafeCopyMatrix (coerceMatrix c) b
            trmm alpha t' (coerceMatrix c)

        (Upper,_,Right (t',r)) ->
            let b1 = unsafeSubmatrixView b (0,0)          (numCols t',numCols b)
                b2 = unsafeSubmatrixView b (numCols t',0) (numCols r ,numCols b)
            in do
                unsafeCopyMatrix c b1
                trmm alpha t' c
                unsafeDoSApplyAddMat alpha r b2 1 c
  where
    (u,d,a) = triToBase t


toLower :: (BaseMatrix a e) => Diag -> a (m,n) e 
        -> Either (Tri a (m,m) e) 
                  (Tri a (n,n) e, a (d,n) e)
toLower diag a =
    if m <= n
        then Left $  triFromBase Lower diag (unsafeSubmatrixView a (0,0) (m,m))
        else let t = triFromBase Lower diag (unsafeSubmatrixView a (0,0) (n,n))
                 r = unsafeSubmatrixView a (n,0) (d,n)
             in Right (t,r)
  where
    (m,n) = shape a
    d     = m - n
    
toUpper :: (BaseMatrix a e) => Diag -> a (m,n) e
        -> Either (Tri a (n,n) e)
                  (Tri a (m,m) e, a (m,d) e)
toUpper diag a =
    if n <= m
        then Left $  triFromBase Upper diag (unsafeSubmatrixView a (0,0) (n,n))
        else let t = triFromBase Upper diag (unsafeSubmatrixView a (0,0) (m,m))
                 r = unsafeSubmatrixView a (0,m) (m,d)
             in Right (t,r)
  where
    (m,n) = shape a
    d     = n - m

trmv :: (ReadMatrix a e m, WriteVector y e m) => 
    e -> Tri a (k,k) e -> y n e -> m ()
trmv alpha t x 
    | dim x == 0 = 
        return ()
        
    | isConj x =
        let (u,d,a) = triToBase t
            side    = rightSide
            (h,u')  = if isHermMatrix a then (NoTrans, flipUpLo u) else (ConjTrans, u)
            uploA   = cblasUpLo u'
            transA  = cblasTrans h
            diagA   = cblasDiag d
            m       = 1
            n       = dim x
            alpha'  = conj alpha
            ldA     = ldaOfMatrix a
            ldB     = stride x
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withVectorPtr x $ \pB ->
                   BLAS.trmm side uploA transA diagA m n alpha' pA ldA pB ldB

    | otherwise =
        let (u,d,a)   = triToBase t
            (transA,u') = if isHermMatrix a then (conjTrans, flipUpLo u) else (noTrans, u)
            uploA     = cblasUpLo u'
            diagA     = cblasDiag d
            n         = dim x
            ldA       = ldaOfMatrix a
            incX      = stride x
        in do
            when (alpha /= 1) $ scaleBy alpha x
            unsafeIOToM $
                withMatrixPtr a $ \pA ->
                withVectorPtr x $ \pX -> do
                   BLAS.trmv uploA transA diagA n pA ldA pX incX


trmm :: (ReadMatrix a e m, WriteMatrix b e m) => 
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
trmm _ _ b
    | numRows b == 0 || numCols b == 0 = return ()
trmm alpha t b =
    let (u,d,a)   = triToBase t
        (h,u')    = if isHermMatrix a then (ConjTrans, flipUpLo u) else (NoTrans, u)
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if isHermMatrix b
                        then (rightSide, flipTrans h, n, m, conj alpha)
                        else (leftSide , h          , m, n, alpha       )
        uploA     = cblasUpLo u'
        transA    = cblasTrans h'
        diagA     = cblasDiag d
        ldA       = ldaOfMatrix a
        ldB       = ldaOfMatrix b
    in unsafeIOToM $
           withMatrixPtr a $ \pA ->
           withMatrixPtr b $ \pB ->
               BLAS.trmm side uploA transA diagA m' n' alpha' pA ldA pB ldB

unsafeDoSSolveTriMatrix :: (ReadMatrix a e m,
    ReadVector y e m, WriteVector x e m) =>
        e -> Tri a (k,l) e -> y k e -> x l e -> m ()
unsafeDoSSolveTriMatrix alpha t y x =
    case (u, toLower d a, toUpper d a) of
        (Lower,Left t',_) -> do
            unsafeCopyVector x (coerceVector y)
            trsv alpha t' (coerceVector x)
            
        (Lower,Right (t',_),_) -> do
            let y1 = unsafeSubvectorView y 0            (numRows t')
            unsafeCopyVector x y1
            trsv alpha t' x
            
        (Upper,_,Left t') -> do
            unsafeCopyVector x (coerceVector y)
            trsv alpha t' x

        (Upper,_,Right (t',r)) ->
            let x1 = unsafeSubvectorView x 0            (numCols t')
                x2 = unsafeSubvectorView x (numCols t') (numCols r)
            in do
                unsafeCopyVector x1 y
                trsv alpha t' x1
                setZero x2
  where
    (u,d,a) = triToBase t


unsafeDoSSolveMatTriMatrix :: (ReadMatrix a e m,
    ReadMatrix c e m, WriteMatrix b e m) =>
        e -> Tri a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
unsafeDoSSolveMatTriMatrix alpha t c b =
    case (u, toLower d a, toUpper d a) of
        (Lower,Left t',_) -> do
            unsafeCopyMatrix b (coerceMatrix c)
            trsm alpha t' (coerceMatrix b)
            
        (Lower,Right (t',_),_) -> do
            let c1 = unsafeSubmatrixView c (0,0)          (numRows t',numCols c)
            unsafeCopyMatrix b c1
            trsm alpha t' b
            
        (Upper,_,Left t') -> do
            unsafeCopyMatrix (coerceMatrix b) c
            trsm alpha t' (coerceMatrix b)

        (Upper,_,Right (t',r)) ->
            let b1 = unsafeSubmatrixView b (0,0)          (numCols t',numCols b)
                b2 = unsafeSubmatrixView b (numCols t',0) (numCols r ,numCols b)
            in do
                unsafeCopyMatrix b1 c
                trsm alpha t' b1
                setZero b2
  where
    (u,d,a) = triToBase t


trsv :: (ReadMatrix a e m, WriteVector y e m) => 
    e -> Tri a (k,k) e -> y n e -> m ()
trsv alpha t x
    | dim x == 0 = return ()

    | isConj x =
        let (u,d,a) = triToBase t
            side    = rightSide
            (h,u')  = if isHermMatrix a then (NoTrans, flipUpLo u) else (ConjTrans, u)
            uploA   = cblasUpLo u'
            transA  = cblasTrans h
            diagA   = cblasDiag d
            m       = 1
            n       = dim x
            alpha'  = conj alpha
            ldA     = ldaOfMatrix a
            ldB     = stride x
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withVectorPtr x $ \pB ->
                   BLAS.trsm side uploA transA diagA m n alpha' pA ldA pB ldB

    | otherwise =
        let (u,d,a) = triToBase t
            (transA,u') = if isHermMatrix a then (conjTrans, flipUpLo u) else (noTrans, u)
            uploA     = cblasUpLo u'
            diagA     = cblasDiag d
            n         = dim x
            ldA       = ldaOfMatrix a
            incX      = stride x
        in do
            when (alpha /= 1) $ scaleBy alpha x
            unsafeIOToM $
                withMatrixPtr a $ \pA ->
                withVectorPtr x $ \pX ->
                    BLAS.trsv uploA transA diagA n pA ldA pX incX

trsm :: (ReadMatrix a e m, WriteMatrix b e m) => 
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
trsm _ _ b
    | numRows b == 0 || numCols b == 0 = return ()
trsm alpha t b =
    let (u,d,a)   = triToBase t
        (h,u')    = if isHermMatrix a then (ConjTrans, flipUpLo u) else (NoTrans, u)
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if isHermMatrix b
                        then (rightSide, flipTrans h, n, m, conj alpha)
                        else (leftSide , h          , m, n, alpha     )
        uploA     = cblasUpLo u'
        transA    = cblasTrans h'
        diagA     = cblasDiag d
        ldA       = ldaOfMatrix a
        ldB       = ldaOfMatrix b
    in unsafeIOToM $     
           withMatrixPtr a $ \pA ->
           withMatrixPtr b $ \pB -> do
               BLAS.trsm side uploA transA diagA m' n' alpha' pA ldA pB ldB


------------------------------------ MSolve ------------------------------

class (MatrixShaped a e, Monad m) => MSolve a e m where
    unsafeDoSolve :: (ReadVector y e m, WriteVector x e m) =>
        a (k,l) e -> y k e -> x l e -> m ()
    unsafeDoSolve = unsafeDoSSolve 1
    
    unsafeDoSolveMat :: (ReadMatrix c e m, WriteMatrix b e m) =>
        a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
    unsafeDoSolveMat = unsafeDoSSolveMat 1
    
    unsafeDoSSolve :: (ReadVector y e m, WriteVector x e m) =>
        e -> a (k,l) e -> y k e -> x l e -> m ()
    unsafeDoSSolve alpha a y x = do
        unsafeDoSolve a y x
        scaleBy alpha x
    
    unsafeDoSSolveMat :: (ReadMatrix c e m, WriteMatrix b e m) =>
        e -> a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
    unsafeDoSSolveMat alpha a c b = do
        unsafeDoSolveMat a c b
        scaleBy alpha b

    unsafeDoSolve_ :: (WriteVector x e m) => a (k,k) e -> x k e -> m ()
    unsafeDoSolve_ = unsafeDoSSolve_ 1

    unsafeDoSSolve_ :: (WriteVector x e m) => e -> a (k,k) e -> x k e -> m ()
    unsafeDoSSolve_ alpha a x = do
        scaleBy alpha x
        unsafeDoSolve_ a x
        
    unsafeDoSolveMat_ :: (WriteMatrix b e m) => a (k,k) e -> b (k,l) e -> m ()
    unsafeDoSolveMat_ = unsafeDoSSolveMat_ 1
        
    unsafeDoSSolveMat_ :: (WriteMatrix b e m) => e -> a (k,k) e -> b (k,l) e -> m ()         
    unsafeDoSSolveMat_ alpha a b = do
        scaleBy alpha b
        unsafeDoSolveMat_ a b


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

newtype STMatrix s mn e = ST (IOMatrix mn e)

unsafeIOMatrixToSTMatrix :: IOMatrix mn e -> STMatrix s mn e
unsafeIOMatrixToSTMatrix = ST
{-# INLINE unsafeIOMatrixToSTMatrix #-}

unsafeSTMatrixToIOMatrix :: STMatrix s mn e -> IOMatrix mn e
unsafeSTMatrixToIOMatrix (ST x) = x
{-# INLINE unsafeSTMatrixToIOMatrix #-}

instance HasVectorView IOMatrix where
    type VectorView IOMatrix = IOVector

instance HasVectorView (STMatrix s) where
    type VectorView (STMatrix s) = STVector s

instance (Storable e) => BaseMatrix_ IOMatrix e where
    matrixViewArray f p m n = DM f p m n
    arrayFromMatrix (DM f p m n l h) = (f,p,m,n,l,h)

instance (Storable e) => BaseMatrix_ (STMatrix s) e where
    matrixViewArray f p m n l h = ST $ DM f p m n l h
    arrayFromMatrix (ST (DM f p m n l h)) = (f,p,m,n,l,h)

instance (Storable e) => BaseMatrix IOMatrix e
instance (Storable e) => BaseMatrix (STMatrix s) e

instance (Storable e) => BaseTensor IOMatrix (Int,Int) e where
    shape  = shapeMatrix
    bounds = boundsMatrix

instance (Storable e) => BaseTensor (STMatrix s) (Int,Int) e where
    shape  = shapeMatrix
    bounds = boundsMatrix

instance (BLAS3 e) => ReadTensor IOMatrix (Int,Int) e IO where
    getSize        = getSizeMatrix
    getAssocs      = getAssocsMatrix
    getIndices     = getIndicesMatrix
    getElems       = getElemsMatrix
    getAssocs'     = getAssocsMatrix'
    getIndices'    = getIndicesMatrix'
    getElems'      = getElemsMatrix'
    unsafeReadElem = unsafeReadElemMatrix

instance (BLAS3 e) => ReadTensor (STMatrix s) (Int,Int) e (ST s) where
    getSize        = getSizeMatrix
    getAssocs      = getAssocsMatrix
    getIndices     = getIndicesMatrix
    getElems       = getElemsMatrix
    getAssocs'     = getAssocsMatrix'
    getIndices'    = getIndicesMatrix'
    getElems'      = getElemsMatrix'
    unsafeReadElem = unsafeReadElemMatrix

instance (BLAS3 e) => WriteTensor IOMatrix (Int,Int) e IO where
    setConstant     = setConstantMatrix
    setZero         = setZeroMatrix
    modifyWith      = modifyWithMatrix
    unsafeWriteElem = unsafeWriteElemMatrix
    canModifyElem   = canModifyElemMatrix
    doConj          = doConjMatrix
    scaleBy         = scaleByMatrix
    shiftBy         = shiftByMatrix

instance (BLAS3 e) => WriteTensor (STMatrix s) (Int,Int) e (ST s) where
    setConstant     = setConstantMatrix
    setZero         = setZeroMatrix
    modifyWith      = modifyWithMatrix
    unsafeWriteElem = unsafeWriteElemMatrix
    canModifyElem   = canModifyElemMatrix
    doConj          = doConjMatrix
    scaleBy         = scaleByMatrix
    shiftBy         = shiftByMatrix

instance (Storable e) => MatrixShaped IOMatrix e where
    herm = hermMatrix

instance (Storable e) => MatrixShaped (STMatrix s) e where
    herm = hermMatrix
    
instance (BLAS3 e) => ReadMatrix IOMatrix     e IO where
instance (BLAS3 e) => ReadMatrix (STMatrix s) e (ST s) where    
    
instance (BLAS3 e) => WriteMatrix IOMatrix     e IO where
instance (BLAS3 e) => WriteMatrix (STMatrix s) e (ST s) where

instance (BLAS3 e) => MMatrix IOMatrix e IO where
    unsafeDoSApplyAdd    = gemv
    unsafeDoSApplyAddMat = gemm
    unsafeGetRow         = unsafeGetRowMatrix
    unsafeGetCol         = unsafeGetColMatrix

instance (BLAS3 e) => MMatrix (STMatrix s) e (ST s) where
    unsafeDoSApplyAdd    = gemv
    unsafeDoSApplyAddMat = gemm
    unsafeGetRow         = unsafeGetRowMatrix
    unsafeGetCol         = unsafeGetColMatrix

instance (BLAS3 e) => MMatrix (Herm (STMatrix s)) e (ST s) where
    unsafeDoSApplyAdd    = hemv'
    unsafeDoSApplyAddMat = hemm'

instance (BLAS3 e) => MMatrix (Herm IOMatrix) e IO where
    unsafeDoSApplyAdd    = hemv'
    unsafeDoSApplyAddMat = hemm'

instance (BLAS3 e) => MMatrix (Tri IOMatrix) e IO where
    unsafeDoSApplyAdd    = unsafeDoSApplyAddTriMatrix
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatTriMatrix
    unsafeDoSApply_      = trmv
    unsafeDoSApplyMat_   = trmm

instance (BLAS3 e) => MMatrix (Tri (STMatrix s)) e (ST s) where
    unsafeDoSApplyAdd    = unsafeDoSApplyAddTriMatrix
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatTriMatrix
    unsafeDoSApply_      = trmv
    unsafeDoSApplyMat_   = trmm

instance (BLAS3 e) => MSolve (Tri IOMatrix) e IO where
    unsafeDoSSolve     = unsafeDoSSolveTriMatrix
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriMatrix
    unsafeDoSSolve_    = trsv
    unsafeDoSSolveMat_ = trsm

instance (BLAS3 e) => MSolve (Tri (STMatrix s)) e (ST s) where
    unsafeDoSSolve     = unsafeDoSSolveTriMatrix
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriMatrix
    unsafeDoSSolve_    = trsv
    unsafeDoSSolveMat_ = trsm
