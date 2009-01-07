{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances,
        RankNTypes #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Base
    where

import Control.Monad
import Foreign
import Unsafe.Coerce

import BLAS.Internal( checkBinaryOp, checkedSubmatrix, checkedDiag,
    checkedRow, checkedCol )

import Data.Elem.BLAS( Complex, Elem, BLAS1, BLAS3, conjugate )

import Data.Tensor.Class
import Data.Tensor.Class.ITensor
import Data.Tensor.Class.MTensor

import Data.Matrix.Class
import Data.Vector.Dense.Base

import Data.Matrix.Dense.IOBase

newtype Matrix np e = Matrix (IOMatrix np e)

freezeIOMatrix :: (BLAS1 e) => IOMatrix np e -> IO (Matrix np e)
freezeIOMatrix x = do
    y <- newCopyIOMatrix x
    return (Matrix y)

thawIOMatrix :: (BLAS1 e) => Matrix np e -> IO (IOMatrix np e)
thawIOMatrix (Matrix x) =
    newCopyIOMatrix x

unsafeFreezeIOMatrix :: IOMatrix np e -> IO (Matrix np e)
unsafeFreezeIOMatrix = return . Matrix

unsafeThawIOMatrix :: Matrix np e -> IO (IOMatrix np e)
unsafeThawIOMatrix (Matrix x) = return x

class (HasVectorView a, Elem e, MatrixShaped a e
      , BaseVector (VectorView a) e) => BaseMatrix a e where
          
    ldaMatrix :: a (n,p) e -> Int
    isHermMatrix :: a (n,p) e -> Bool

    -- | Cast the shape type of the matrix.
    coerceMatrix :: a np e -> a np' e
    coerceMatrix = unsafeCoerce
    {-# INLINE coerceMatrix #-}

    unsafeSubmatrixView :: a (n,p) e -> (Int,Int) -> (Int,Int) -> a (n',p') e
    unsafeDiagView :: a (n,p) e -> Int -> VectorView a k e
    unsafeRowView :: a (n,p) e -> Int -> VectorView a p e
    unsafeColView :: a (n,p) e -> Int -> VectorView a n e

    -- | Possibly create a vector view of a matrix.
    maybeToVectorView :: a (n,p) e -> Maybe (VectorView a np e)
    
    -- | Possibly create a matrix view of a row vector.  This will fail if the
    -- vector is conjugated and the stride is not @1@.
    maybeViewAsRow  :: VectorView a p e -> Maybe (a (one,p) e)
    
    -- | Possibly create a matrix view of a column vector.  This will fail
    -- if the stride of the vector is not @1@ and the vector is not conjugated.
    maybeViewAsCol  :: VectorView a n e -> Maybe (a (n,one) e)

    withMatrixPtrIO :: a np e -> (Ptr e -> IO r) -> IO r

    unsafeIOMatrixToMatrix :: IOMatrix np e -> a np e
    unsafeMatrixToIOMatrix :: a np e -> IOMatrix np e


class (BaseMatrix a e, BLAS3 e, ReadTensor a (Int,Int) e m
      -- , MMatrix a e m, MMatrix (Herm a) e m, MMatrix (Tri a) e m
      -- , MSolve (Tri a) e m
      , ReadVector (VectorView a) e m) => ReadMatrix a e m where
          
    -- | Creates a new matrix of the given shape.  The elements will be 
    -- uninitialized.
    newMatrix_ :: (Int,Int) -> m (a (n,p) e)

    -- | Execture an @IO@ action with a pointer to the first element and then
    -- convert the @IO@ action to an action in the monad @m@.
    withMatrixPtr :: a (n,p) e -> (Ptr e -> IO r) -> m r

    -- | Convert a mutable matrix to an immutable one by taking a complete
    -- copy of it.
    freezeMatrix :: a (n,p) e -> m (Matrix (n,p) e)

    -- | Convert an immutable matrix to a mutable one by taking a complete
    -- copy of it.
    thawMatrix :: Matrix (n,p) e -> m (a (n,p) e)

    unsafeFreezeMatrix :: a (n,p) e -> m (Matrix (n,p) e)
    unsafeThawMatrix :: Matrix (n,p) e -> m (a (n,p) e)

class (ReadMatrix a e m, WriteTensor a (Int,Int) e m
      , WriteVector (VectorView a) e m) =>
    WriteMatrix a e m where

-- | Creates a new matrix with the given association list.  Unspecified
-- indices will get initialized to zero.
newMatrix :: (ReadMatrix a e m) => 
    (Int,Int) -> [((Int,Int), e)] -> m (a (n,p) e)
newMatrix (m,n) ies = do
    a <- newZeroMatrix (m,n)
    withMatrixPtr a $ \p ->
        forM_ ies $ \((i,j),e) -> do
            when (i < 0 || i >= m || j < 0 || j >= n) $ fail $
                "Index `" ++ show (i,j) ++ 
                    "' is invalid for a matrix with shape `" ++ show (m,n) ++
                    "'"
            pokeElemOff p (i+j*m) e
    return a
{-# INLINE newMatrix #-}
    
unsafeNewMatrix :: (ReadMatrix a e m) => 
    (Int,Int) -> [((Int,Int), e)] -> m (a (n,p) e)
unsafeNewMatrix (m,n) ies = do
    a <- newZeroMatrix (m,n)
    withMatrixPtr a $ \p ->
        forM_ ies $ \((i,j),e) -> do
            pokeElemOff p (i+j*m) e
    return a
{-# INLINE unsafeNewMatrix #-}

-- | Create a new matrix with the given elements in column-major order.
newListMatrix :: (ReadMatrix a e m) => (Int,Int) -> [e] -> m (a (n,p) e)
newListMatrix (m,n) es = do
    a <- newZeroMatrix (m,n)
    withMatrixPtr a $ flip pokeArray (take (m*n) es)
    return a
{-# INLINE newListMatrix #-}

-- | Form a matrix from a list of column vectors.
newColsMatrix :: (ReadVector x e m, ReadMatrix a e m) => 
    (Int,Int) -> [x n e] -> m (a (n,p) e)
newColsMatrix (m,n) cs = do
    a <- newZeroMatrix (m,n)
    forM_ (zip [0..(n-1)] cs) $ \(j,c) ->
        unsafeCopyVector (unsafeColView a j) c
    return a
{-# INLINE newColsMatrix #-}

-- | Form a matrix from a list of row vectors.
newRowsMatrix :: (ReadVector x e m, ReadMatrix a e m) => 
    (Int,Int) -> [x p e] -> m (a (n,p) e)
newRowsMatrix (m,n) rs = do
    a <- newZeroMatrix (m,n)
    forM_ (zip [0..(m-1)] rs) $ \(i,r) ->
        unsafeCopyVector (unsafeRowView a i) r
    return a
{-# INLINE newRowsMatrix #-}

-- | Create a new matrix from a column vector.
newColMatrix :: (ReadVector x e m, ReadMatrix a e m) => 
    x n e -> m (a (n,one) e)
newColMatrix x = newColsMatrix (dim x,1) [x]
{-# INLINE newColMatrix #-}

-- | Create a new matrix from a row vector.
newRowMatrix :: (ReadVector x e m, ReadMatrix a e m) => 
    x p e -> m (a (one,p) e)
newRowMatrix x = newRowsMatrix (1,dim x) [x]
{-# INLINE newRowMatrix #-}

-- | Create a zero matrix of the specified shape.
newZeroMatrix :: (ReadMatrix a e m) => (Int,Int) -> m (a (n,p) e)
newZeroMatrix mn = do
    a <- newMatrix_ mn
    unsafeSetZeroMatrix a
    return a
{-# INLINE newZeroMatrix #-}

-- | Set every element in the matrix to zero.
setZeroMatrix :: (WriteMatrix a e m) => a (n,p) e -> m ()
setZeroMatrix = unsafeSetZeroMatrix
{-# INLINE setZeroMatrix #-}

unsafeSetZeroMatrix :: (ReadMatrix a e m) => a (n,p) e -> m ()
unsafeSetZeroMatrix a =
    withMatrixPtr a $ \_ ->
        setZeroIOMatrix (unsafeMatrixToIOMatrix a)
{-# INLINE unsafeSetZeroMatrix #-}

-- | Create a constant matrix of the specified shape.
newConstantMatrix :: (ReadMatrix a e m) => (Int,Int) -> e -> m (a (n,p) e)
newConstantMatrix mn e = do
    a <- newMatrix_ mn
    unsafeSetConstantMatrix e a
    return a
{-# INLINE newConstantMatrix #-}

-- | Set every element in the matrix to the given constant.
setConstantMatrix :: (WriteMatrix a e m) => e -> a (n,p) e -> m ()
setConstantMatrix = unsafeSetConstantMatrix
{-# INLINE setConstantMatrix #-}

unsafeSetConstantMatrix :: (ReadMatrix a e m) => e -> a (n,p) e -> m ()
unsafeSetConstantMatrix e a =
    withMatrixPtr a $ \_ ->
        setConstantIOMatrix e (unsafeMatrixToIOMatrix a)
{-# INLINE unsafeSetConstantMatrix #-}

-- | Create a new matrix of the given shape with ones along the diagonal, 
-- and zeros everywhere else.
newIdentityMatrix :: (ReadMatrix a e m) => (Int,Int) -> m (a (n,p) e)
newIdentityMatrix mn = do
    a <- newMatrix_ mn
    unsafeSetIdentityMatrix a
    return a
{-# INLINE newIdentityMatrix #-}

-- | Set diagonal elements to one and all other elements to zero.
setIdentityMatrix :: (WriteMatrix a e m) => a (n,p) e -> m ()
setIdentityMatrix = unsafeSetIdentityMatrix
{-# INLINE setIdentityMatrix #-}

unsafeSetIdentityMatrix :: (ReadMatrix a e m) => a (n,p) e -> m ()
unsafeSetIdentityMatrix a = do
    unsafeSetZeroMatrix a
    unsafeSetConstantVector 1 (unsafeDiagView a 0)
{-# INLINE unsafeSetIdentityMatrix #-}

-- | Get a copy of a matrix.
newCopyMatrix :: (ReadMatrix a e m, ReadMatrix b e m) => 
    a (n,p) e -> m (b (n,p) e)
newCopyMatrix a | isHermMatrix a = liftM herm $ newCopyMatrix (herm a)
                | otherwise      = do
    b <- newMatrix_ (shape a)
    unsafeCopyMatrix b a
    return b
{-# INLINE newCopyMatrix #-}

-- | Get a copy of a matrix and make sure the returned matrix is not
-- a view.  Specififially, the returned matrix will have @isHermMatrix@
-- equal to @False@.
newCopyMatrix' :: (ReadMatrix a e m, ReadMatrix b e m) =>
    a (n,p) e -> m (b (n,p) e)
newCopyMatrix' a = do
    b <- newMatrix_ (shape a)
    unsafeCopyMatrix b a
    return b
{-# INLINE newCopyMatrix' #-}

-- | @copyMatrix dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyMatrix :: (WriteMatrix b e m, ReadMatrix a e m) => 
    b (n,p) e -> a (n,p) e -> m ()
copyMatrix b a = checkBinaryOp (shape b) (shape a) $ unsafeCopyMatrix b a
{-# INLINE copyMatrix #-}

unsafeCopyMatrix :: (ReadMatrix b e m, ReadMatrix a e m) => 
    b (n,p) e -> a (n,p) e -> m ()
unsafeCopyMatrix = liftMatrix2 unsafeCopyVector
{-# INLINE unsafeCopyMatrix #-}

-- | @swapMatrix x y@ swaps the values stored in two matrices.
swapMatrix :: (WriteMatrix a e m) => 
    a (n,p) e -> a (n,p) e -> m ()
swapMatrix a b = checkBinaryOp (shape b) (shape a) $ unsafeSwapMatrix a b
{-# INLINE swapMatrix #-}

unsafeSwapMatrix :: (ReadMatrix a e m, ReadMatrix b e m) => 
    a (n,p) e -> b (n,p) e -> m ()
unsafeSwapMatrix = liftMatrix2 unsafeSwapVector
{-# INLINE unsafeSwapMatrix #-}

-- | Swap the elements in two rows of a matrix.
swapRows :: (WriteMatrix a e m) => a (n,p) e -> Int -> Int -> m ()
swapRows a i j = 
    when (i /= j) $ unsafeSwapVector (rowView a i) (rowView a j)
{-# INLINE swapRows #-}

-- | Swap the elements in two columns of a matrix.
swapCols :: (WriteMatrix a e m) => a (n,p) e -> Int -> Int -> m ()
swapCols a i j = 
    when (i /= j) $ unsafeSwapVector (colView a i) (colView a j)
{-# INLINE swapCols #-}

unsafeSwapRows :: (ReadMatrix a e m) => a (n,p) e -> Int -> Int -> m ()
unsafeSwapRows a i j = 
    when (i /= j) $ unsafeSwapVector (unsafeRowView a i) (unsafeRowView a j)
{-# INLINE unsafeSwapRows #-}

unsafeSwapCols :: (ReadMatrix a e m) => a (n,p) e -> Int -> Int -> m ()
unsafeSwapCols a i j = 
    when (i /= j) $ unsafeSwapVector (unsafeColView a i) (unsafeColView a j)
{-# INLINE unsafeSwapCols #-}

-- | @submatrixView a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrixView :: (BaseMatrix a e) => a (n,p) e -> (Int,Int) -> (Int,Int) -> a (n',p') e
submatrixView a = checkedSubmatrix (shape a) (unsafeSubmatrixView a)
{-# INLINE submatrixView #-}

-- | Divide the rows of a matrix into two blocks and return views into the
-- blocks.  The integer argument indicates how many rows should be in the
-- first block.
splitRowsAt :: (BaseMatrix a e) =>
    Int -> a (n,p) e -> (a (n1,p) e, a (n2,p) e)
splitRowsAt m1 a = ( submatrixView a (0,0)  (m1,n)
                   , submatrixView a (m1,0) (m2,n)
                   )
  where 
    (m,n) = shape a
    m2    = m - m1
{-# INLINE splitRowsAt #-}

unsafeSplitRowsAt :: (BaseMatrix a e) =>
    Int -> a (n,p) e -> (a (n1,p) e, a (n2,p) e)
unsafeSplitRowsAt m1 a = ( unsafeSubmatrixView a (0,0)  (m1,n)
                         , unsafeSubmatrixView a (m1,0) (m2,n)
                         )
  where 
    (m,n) = shape a
    m2    = m - m1
{-# INLINE unsafeSplitRowsAt #-}

-- | Divide the columns of a matrix into two blocks and return views into the
-- blocks.  The integer argument indicates how many columns should be in the
-- first block.
splitColsAt :: (BaseMatrix a e) =>
    Int -> a (n,p) e -> (a (n,p1) e, a (n,p2) e)
splitColsAt n1 a = ( submatrixView a (0,0)  (m,n1)
                   , submatrixView a (0,n1) (m,n2)
                   )
  where
    (m,n) = shape a
    n2    = n - n1
{-# INLINE splitColsAt #-}

unsafeSplitColsAt :: (BaseMatrix a e) =>
    Int -> a (n,p) e -> (a (n,p1) e, a (n,p2) e)
unsafeSplitColsAt n1 a = ( unsafeSubmatrixView a (0,0)  (m,n1)
                         , unsafeSubmatrixView a (0,n1) (m,n2)
                         )
  where
    (m,n) = shape a
    n2    = n - n1
{-# INLINE unsafeSplitColsAt #-}

-- | Get a list of vector views of the rows of the matrix.
rowViews :: (BaseMatrix a e) => a (n,p) e -> [VectorView a p e]
rowViews a = [ unsafeRowView a i | i <- [0..numRows a - 1] ]
{-# INLINE rowViews #-}

-- | Get a list of vector views of the columns of the matrix.
colViews :: (BaseMatrix a e) => a (n,p) e -> [VectorView a n e]
colViews a = [ unsafeColView a j | j <- [0..numCols a - 1] ]
{-# INLINE colViews #-}

-- | Get a vector view of the given row in a matrix.
rowView :: (BaseMatrix a e) => a (n,p) e -> Int -> VectorView a p e
rowView a = checkedRow (shape a) (unsafeRowView a)
{-# INLINE rowView #-}

unsafeGetRowMatrix :: (ReadMatrix a e m, ReadVector y e m) => 
    a (n,p) e -> Int -> m (y p e)
unsafeGetRowMatrix a i = newCopyVector (unsafeRowView a i)
{-# INLINE unsafeGetRowMatrix #-}

-- | Get a vector view of the given column in a matrix.
colView :: (BaseMatrix a e) => a (n,p) e -> Int -> VectorView a n e
colView a = checkedCol (shape a) (unsafeColView a)
{-# INLINE colView #-}

unsafeGetColMatrix :: (ReadMatrix a e m, ReadVector y e m) => 
    a (n,p) e -> Int -> m (y n e)
unsafeGetColMatrix a j = newCopyVector (unsafeColView a j)
{-# INLINE unsafeGetColMatrix #-}

-- | Get a vector view of the given diagonal in a matrix.
diagView :: (BaseMatrix a e) => a (n,p) e -> Int -> VectorView a k e
diagView a = checkedDiag (shape a) (unsafeDiagView a)
{-# INLINE diagView #-}

-- | Get the given diagonal in a matrix.  Negative indices correspond
-- to sub-diagonals.
getDiag :: (ReadMatrix a e m, WriteVector y e m) => 
    a (n,p) e -> Int -> m (y k e)
getDiag a = checkedDiag (shape a) (unsafeGetDiag a)
{-# INLINE getDiag #-}

-- | Same as 'getDiag' but not range-checked.
unsafeGetDiag :: (ReadMatrix a e m, WriteVector y e m) => 
    a (n,p) e -> Int -> m (y k e)
unsafeGetDiag a i = newCopyVector (unsafeDiagView a i)
{-# INLINE unsafeGetDiag #-}

-- | Conjugate every element of a matrix.
doConjMatrix :: (WriteMatrix a e m) => a (n,p) e -> m ()
doConjMatrix = unsafeDoConjMatrix
{-# INLINE doConjMatrix #-}

unsafeDoConjMatrix :: (ReadMatrix a e m) => a (n,p) e -> m ()
unsafeDoConjMatrix = liftMatrix unsafeDoConjVector
{-# INLINE unsafeDoConjMatrix #-}

-- | Get a new matrix with elements with the conjugates of the elements
-- of the given matrix.
getConjMatrix :: (ReadMatrix a e m, WriteMatrix b e m) =>
    a (n,p) e -> m (b (n,p) e)
getConjMatrix = getUnaryMatrixOp doConjMatrix
{-# INLINE getConjMatrix #-}

-- | Scale every element of a matrix by the given value.
scaleByMatrix :: (WriteMatrix a e m) => e -> a (n,p) e -> m ()
scaleByMatrix = unsafeScaleByMatrix
{-# INLINE scaleByMatrix #-}

unsafeScaleByMatrix :: (ReadMatrix a e m) => e -> a (n,p) e -> m ()
unsafeScaleByMatrix k = liftMatrix (unsafeScaleByVector k)
{-# INLINE unsafeScaleByMatrix #-}

-- | Get a new matrix by scaling the elements of another matrix
-- by a given value.
getScaledMatrix :: (ReadMatrix a e m, WriteMatrix b e m) =>
    e -> a (n,p) e -> m (b (n,p) e)
getScaledMatrix e = getUnaryMatrixOp (scaleByMatrix e)
{-# INLINE getScaledMatrix #-}

-- | Add a constant to every element in a matrix.
shiftByMatrix :: (WriteMatrix a e m) => e -> a (n,p) e -> m ()
shiftByMatrix = unsafeShiftByMatrix
{-# INLINE shiftByMatrix #-}

unsafeShiftByMatrix :: (WriteMatrix a e m) => e -> a (n,p) e -> m ()
unsafeShiftByMatrix k = liftMatrix (unsafeShiftByVector k)
{-# INLINE unsafeShiftByMatrix #-}

-- | Get a new matrix by shifting the elements of another matrix
-- by a given value.
getShiftedMatrix :: (ReadMatrix a e m, WriteMatrix b e m) =>
    e -> a (n,p) e -> m (b (n,p) e)
getShiftedMatrix e = getUnaryMatrixOp (shiftByMatrix e)
{-# INLINE getShiftedMatrix #-}

-- | Replace the first argument with the elementwise sum.
addMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b (n,p) e -> a (n,p) e -> m ()
addMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeAddMatrix b a
{-# INLINE addMatrix #-}

unsafeAddMatrix :: (ReadMatrix b e m, ReadMatrix a e m) =>
    b (n,p) e -> a (n,p) e -> m ()
unsafeAddMatrix b a = unsafeAxpyMatrix 1 a b
{-# INLINE unsafeAddMatrix #-}

-- | @getAddMatrix a b@ creates a new matrix equal to the sum @a+b@.  The 
-- operands must have the same shape.
getAddMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, ReadMatrix c e m) => 
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
getAddMatrix = checkMatrixOp2 unsafeGetAddMatrix
{-# INLINE getAddMatrix #-}

unsafeGetAddMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, ReadMatrix c e m) => 
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
unsafeGetAddMatrix = unsafeGetBinaryMatrixOp unsafeAddMatrix
{-# INLINE unsafeGetAddMatrix #-}

-- | Replace the first argument with the elementwise sum.
subMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b (n,p) e -> a (n,p) e -> m ()    
subMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeSubMatrix b a
{-# INLINE subMatrix #-}

unsafeSubMatrix :: (ReadMatrix b e m, ReadMatrix a e m) =>
    b (n,p) e -> a (n,p) e -> m ()
unsafeSubMatrix b a = unsafeAxpyMatrix (-1) a b
{-# INLINE unsafeSubMatrix #-}

-- | @getSubMatrix a b@ creates a new matrix equal to the difference @a-b@.  The 
-- operands must have the same shape.
getSubMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, ReadMatrix c e m) => 
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
getSubMatrix = checkMatrixOp2 unsafeGetSubMatrix
{-# INLINE getSubMatrix #-}

unsafeGetSubMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, ReadMatrix c e m) => 
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
unsafeGetSubMatrix = unsafeGetBinaryMatrixOp unsafeSubMatrix
{-# INLINE unsafeGetSubMatrix #-}

-- | @axpyMatrix a x y@ replaces @y := a x + y@.
axpyMatrix :: (ReadMatrix a e m, WriteMatrix b e m) =>
    e -> a (n,p) e -> b (n,p) e -> m ()
axpyMatrix alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyMatrix alpha x y
{-# INLINE axpyMatrix #-}

unsafeAxpyMatrix :: (ReadMatrix a e m, ReadMatrix b e m) =>
    e -> a (n,p) e -> b (n,p) e -> m ()
unsafeAxpyMatrix alpha = liftMatrix2 (unsafeAxpyVector alpha)
{-# INLINE unsafeAxpyMatrix #-}

-- | Replace the first argument with the elementwise product.
mulMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b (n,p) e -> a (n,p) e -> m ()    
mulMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeMulMatrix b a
{-# INLINE mulMatrix #-}

unsafeMulMatrix :: (ReadMatrix b e m, ReadMatrix a e m) =>
    b (n,p) e -> a (n,p) e -> m ()
unsafeMulMatrix = liftMatrix2 unsafeMulVector
{-# INLINE unsafeMulMatrix #-}

-- | @getMulMatrix a b@ creates a new matrix equal to the elementwise product 
-- @a*b@.  The operands must have the same shape.
getMulMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, ReadMatrix c e m) => 
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
getMulMatrix = checkMatrixOp2 unsafeGetMulMatrix
{-# INLINE getMulMatrix #-}

unsafeGetMulMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, ReadMatrix c e m) => 
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
unsafeGetMulMatrix = unsafeGetBinaryMatrixOp unsafeMulMatrix
{-# INLINE unsafeGetMulMatrix #-}

-- | Replace the first argument with the elementwise quotient.
divMatrix :: (WriteMatrix b e m, ReadMatrix a e m) =>
    b (n,p) e -> a (n,p) e -> m ()    
divMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeDivMatrix b a
{-# INLINE divMatrix #-}

unsafeDivMatrix :: (ReadMatrix b e m, ReadMatrix a e m) =>
    b (n,p) e -> a (n,p) e -> m ()
unsafeDivMatrix = liftMatrix2 unsafeDivVector
{-# INLINE unsafeDivMatrix #-}

-- | @getDivMatrix a b@ creates a new matrix equal to the elementwise ratio
-- @a/b@.  The operands must have the same shape.
getDivMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, ReadMatrix c e m) => 
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
getDivMatrix = checkMatrixOp2 unsafeGetDivMatrix
{-# INLINE getDivMatrix #-}

unsafeGetDivMatrix :: 
    (ReadMatrix a e m, ReadMatrix b e m, ReadMatrix c e m) => 
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
unsafeGetDivMatrix = unsafeGetBinaryMatrixOp unsafeDivMatrix
{-# INLINE unsafeGetDivMatrix #-}

instance (Elem e) => BaseMatrix IOMatrix e where
    ldaMatrix = ldaMatrixIOMatrix
    {-# INLINE ldaMatrix #-}
    isHermMatrix = isHermIOMatrix
    {-# INLINE isHermMatrix #-}
    unsafeSubmatrixView = unsafeSubmatrixViewIOMatrix
    {-# INLINE unsafeSubmatrixView #-}
    unsafeDiagView = unsafeDiagViewIOMatrix
    {-# INLINE unsafeDiagView #-}
    unsafeRowView = unsafeRowViewIOMatrix
    {-# INLINE unsafeRowView #-}
    unsafeColView = unsafeColViewIOMatrix
    {-# INLINE unsafeColView #-}
    maybeToVectorView = maybeToVectorViewIOMatrix
    {-# INLINE maybeToVectorView #-}
    maybeViewAsRow = maybeViewAsRowIOMatrix
    {-# INLINE maybeViewAsRow #-}    
    maybeViewAsCol = maybeViewAsColIOMatrix
    {-# INLINE maybeViewAsCol #-}
    withMatrixPtrIO = withIOMatrixPtr
    {-# INLINE withMatrixPtrIO #-}
    unsafeIOMatrixToMatrix = id
    {-# INLINE unsafeIOMatrixToMatrix #-}
    unsafeMatrixToIOMatrix = id
    {-# INLINE unsafeMatrixToIOMatrix #-}

instance (BLAS3 e) => ReadMatrix IOMatrix e IO where
    newMatrix_ = newIOMatrix_
    {-# INLINE newMatrix_ #-}
    withMatrixPtr = withIOMatrixPtr
    {-# INLINE withMatrixPtr #-}
    freezeMatrix = freezeIOMatrix
    {-# INLINE freezeMatrix #-}
    unsafeFreezeMatrix = unsafeFreezeIOMatrix
    {-# INLINE unsafeFreezeMatrix #-}
    thawMatrix = thawIOMatrix
    {-# INLINE thawMatrix #-}
    unsafeThawMatrix = unsafeThawIOMatrix
    {-# INLINE unsafeThawMatrix #-}
    

-- | Take a unary elementwise vector operation and apply it to the elements
-- of a matrix.
liftMatrix :: (ReadMatrix a e m) =>
    (forall k. VectorView a k e -> m ()) -> a (n,p) e -> m ()
liftMatrix f a =
    case maybeToVectorView a of
        Just x -> f x
        _ ->
            let xs = case isHermMatrix a of
                          True ->  rowViews (coerceMatrix a)
                          False -> colViews (coerceMatrix a)
            in mapM_ f xs
{-# INLINE liftMatrix #-}

-- | Take a binary elementwise vector operation and apply it to the elements
-- of a pair of matrices.
liftMatrix2 :: (ReadMatrix a e m, ReadMatrix b f m) =>
    (forall k. VectorView a k e -> VectorView b k f -> m ()) ->
        a (n,p) e -> b (n,p) f -> m ()
liftMatrix2 f a b =
    if isHermMatrix a == isHermMatrix b
        then case (maybeToVectorView a, maybeToVectorView b) of
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
{-# INLINE liftMatrix2 #-}

checkMatrixOp2 :: (BaseMatrix x e, BaseMatrix y f) => 
    (x n e -> y n f -> a) ->
        x n e -> y n f -> a
checkMatrixOp2 f x y = 
    checkBinaryOp (shape x) (shape y) $ f x y
{-# INLINE checkMatrixOp2 #-}

getUnaryMatrixOp :: (ReadMatrix a e m, ReadMatrix b e m) =>
    (b (n,p) e -> m ()) -> a (n,p) e -> m (b (n,p) e)
getUnaryMatrixOp f a = do
    b <- newCopyMatrix a
    f b
    return b
{-# INLINE getUnaryMatrixOp #-}

unsafeGetBinaryMatrixOp :: 
    (ReadMatrix c e m, ReadMatrix a e m, ReadMatrix b f m) => 
    (c (n,p) e -> b (n,p) f -> m ()) ->
        a (n,p) e -> b (n,p) f -> m (c (n,p) e)
unsafeGetBinaryMatrixOp f a b = do
    c <- newCopyMatrix a
    f c b
    return c
