{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances,
        TypeFamilies, Rank2Types, ScopedTypeVariables, FunctionalDependencies #-}
{-# OPTIONS_HADDOCK hide #-}
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
import Control.Monad.Interleave
import Control.Monad.ST
import Data.AEq
import Foreign
import System.IO.Unsafe
import Text.Printf

import BLAS.Internal( checkBinaryOp, checkedSubmatrix, checkedDiag,
    checkedRow, checkedCol, inlinePerformIO )

import Data.Elem.BLAS( Elem, BLAS1, BLAS2, BLAS3, conjugate )
import qualified Data.Elem.BLAS.Level1 as BLAS
import qualified Data.Elem.BLAS.Level2 as BLAS
import qualified Data.Elem.BLAS.Level3 as BLAS

import Data.Tensor.Class
import Data.Tensor.Class.ITensor
import Data.Tensor.Class.MTensor

import Data.Matrix.Class
import Data.Matrix.Herm
import Data.Matrix.TriBase

import Data.Vector.Dense.IOBase
import Data.Vector.Dense.Base

import Data.Matrix.Dense.IOBase

-- | Immutable dense matrices.  The type arguments are as follows:
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
newtype Matrix e = Matrix (IOMatrix e)

freezeIOMatrix :: IOMatrix e -> IO (Matrix e)
freezeIOMatrix x@(IOMatrix _ _ _ _ _ _) = do
    y <- newCopyIOMatrix x
    return (Matrix y)
{-# INLINE freezeIOMatrix #-}

thawIOMatrix :: Matrix e -> IO (IOMatrix e)
thawIOMatrix (Matrix x@(IOMatrix _ _ _ _ _ _)) =
    newCopyIOMatrix x
{-# INLINE thawIOMatrix #-}

unsafeFreezeIOMatrix :: IOMatrix e -> IO (Matrix e)
unsafeFreezeIOMatrix = return . Matrix
{-# INLINE unsafeFreezeIOMatrix #-}

unsafeThawIOMatrix :: Matrix e -> IO (IOMatrix e)
unsafeThawIOMatrix (Matrix x) = return x
{-# INLINE unsafeThawIOMatrix #-}

-- | Common functionality for all dense matrix types.
class ( HasVectorView a, 
        MatrixShaped a, 
        HasHerm a,
        BaseVector (VectorView a) ) => BaseMatrix a where
      
    -- | Get the leading dimension of the storage of the matrix.    
    ldaMatrix :: a e -> Int
    
    -- | Get the storage type of the matrix.
    transEnumMatrix :: a e -> TransEnum
    
    -- | Indicate whether or not the underlying matrix storage is
    -- transposed and conjugated.
    isHermMatrix :: a e -> Bool
    isHermMatrix = (ConjTrans ==) . transEnumMatrix
    {-# INLINE isHermMatrix #-}

    unsafeSubmatrixView :: a e -> (Int,Int) -> (Int,Int) -> a e
    unsafeDiagView :: a e -> Int -> VectorView a e
    unsafeRowView :: a e -> Int -> VectorView a e
    unsafeColView :: a e -> Int -> VectorView a e

    -- | Possibly create a vector view of a matrix.  This will fail if the
    -- matrix is hermed or if the lda of the matrix is not equal to the
    -- number of rows in the matrix.
    maybeViewMatrixAsVector :: a e -> Maybe (VectorView a e)
    
    -- | Possibly create a matrix view of a row vector.  This will fail if
    -- the stride of the vector is not @1@ or the vector is conjugated.
    maybeViewVectorAsRow  :: VectorView a e -> Maybe (a e)
    
    -- | Possibly create a matrix view of a column vector.  This will fail
    -- if the stride of the vector is not @1@ or the vector is not conjugated.
    maybeViewVectorAsCol  :: VectorView a e -> Maybe (a e)
    
    -- | Possible create a matrix view of the vector.  This will fail if
    -- the stride of the vector is not @1@ or the vector is conjugated.
    -- An error will be called if the vector does not have the same number
    -- of elements as the desired matrix.
    maybeViewVectorAsMatrix :: (Int,Int) 
                            -> VectorView a e 
                            -> Maybe (a e)
    
    -- | Unsafe cast from a matrix to an 'IOMatrix'.
    unsafeMatrixToIOMatrix :: a e -> IOMatrix e
    unsafeIOMatrixToMatrix :: IOMatrix e -> a e


-- | Dense matrices that can be read in a monad.
class ( BaseMatrix a,
        MonadInterleave m,
        ReadTensor a (Int,Int) m,
        MMatrix a m, 
        MMatrix (Herm a) m, 
        MMatrix (Tri a) m,
        MSolve (Tri a) m,
        ReadVector (VectorView a) m ) => ReadMatrix a m where

    -- | Convert a mutable matrix to an immutable one by taking a complete
    -- copy of it.
    freezeMatrix :: a e -> m (Matrix e)
    unsafeFreezeMatrix :: a e -> m (Matrix e)

    -- | Cast the matrix to an 'IOMatrix', perform an @IO@ action, and
    -- convert the @IO@ action to an action in the monad @m@.
    unsafePerformIOWithMatrix :: a e -> (IOMatrix e -> IO r) -> m r


-- | Dense matrices that can be created or modified in a monad.
class ( ReadMatrix a m, 
        WriteTensor a (Int,Int) m, 
        WriteVector (VectorView a) m ) => WriteMatrix a m | m -> a where

    -- | Creates a new matrix of the given shape.  The elements will be 
    -- uninitialized.
    newMatrix_ :: (Elem e) => (Int,Int) -> m (a e)

    -- | Convert an immutable matrix to a mutable one by taking a complete
    -- copy of it.
    thawMatrix :: Matrix e -> m (a e)
    unsafeThawMatrix :: Matrix e -> m (a e)

    -- | Unsafely convert an 'IO' action that creates an 'IOMatrix' into
    -- an action in @m@ that creates a matrix.
    unsafeConvertIOMatrix :: IO (IOMatrix e) -> m (a e)


-- | Creates a new matrix with the given association list.  Unspecified
-- indices will get initialized to zero.
newMatrix :: (WriteMatrix a m, Elem e)
          => (Int,Int) -> [((Int,Int), e)] -> m (a e)
newMatrix mn ies = unsafeConvertIOMatrix $ 
    newIOMatrix mn ies
{-# INLINE newMatrix #-}
    
unsafeNewMatrix :: (WriteMatrix a m, Elem e)
                => (Int,Int) -> [((Int,Int), e)] -> m (a e)
unsafeNewMatrix mn ies = unsafeConvertIOMatrix $
    unsafeNewIOMatrix mn ies
{-# INLINE unsafeNewMatrix #-}

-- | Create a new matrix with the given elements in column-major order.
newListMatrix :: (WriteMatrix a m, Elem e) 
              => (Int,Int) -> [e] -> m (a e)
newListMatrix mn es = unsafeConvertIOMatrix $
    newListIOMatrix mn es
{-# INLINE newListMatrix #-}

-- | Form a matrix from a list of column vectors.
newColsMatrix :: (ReadVector x m, WriteMatrix a m, Elem e)
              => (Int,Int) -> [x e] -> m (a e)
newColsMatrix mn cs = unsafeConvertIOMatrix $
    newColsIOMatrix mn (map unsafeVectorToIOVector cs)
{-# INLINE newColsMatrix #-}

-- | Form a matrix from a list of row vectors.
newRowsMatrix :: (ReadVector x m, WriteMatrix a m, Elem e)
              => (Int,Int) -> [x e] -> m (a e)
newRowsMatrix mn rs = unsafeConvertIOMatrix $
    newRowsIOMatrix mn (map unsafeVectorToIOVector rs)
{-# INLINE newRowsMatrix #-}

-- | Create a new matrix from a column vector.
newColMatrix :: (ReadVector x m, WriteMatrix a m, Elem e) 
             => x e -> m (a e)
newColMatrix x = unsafeConvertIOMatrix $ 
    newColIOMatrix (unsafeVectorToIOVector x)
{-# INLINE newColMatrix #-}

-- | Create a new matrix from a row vector.
newRowMatrix :: (ReadVector x m, WriteMatrix a m, Elem e) 
             => x e -> m (a e)
newRowMatrix x = unsafeConvertIOMatrix $ 
    newRowIOMatrix (unsafeVectorToIOVector x)
{-# INLINE newRowMatrix #-}

-- | Create a zero matrix of the specified shape.
newZeroMatrix :: (WriteMatrix a m, Elem e) 
              => (Int,Int) -> m (a e)
newZeroMatrix mn = unsafeConvertIOMatrix $ newZeroIOMatrix mn
{-# INLINE newZeroMatrix #-}

-- | Set every element in the matrix to zero.
setZeroMatrix :: (WriteMatrix a m) => a e -> m ()
setZeroMatrix a = unsafePerformIOWithMatrix a $ setZeroIOMatrix
{-# INLINE setZeroMatrix #-}

-- | Create a constant matrix of the specified shape.
newConstantMatrix :: (WriteMatrix a m, Elem e) 
                  => (Int,Int) -> e -> m (a e)
newConstantMatrix mn e = unsafeConvertIOMatrix $ newConstantIOMatrix mn e
{-# INLINE newConstantMatrix #-}

-- | Set every element in the matrix to the given constant.
setConstantMatrix :: (WriteMatrix a m) => e -> a e -> m ()
setConstantMatrix e a = unsafePerformIOWithMatrix a $ setConstantIOMatrix e
{-# INLINE setConstantMatrix #-}

-- | Create a new matrix of the given shape with ones along the diagonal, 
-- and zeros everywhere else.
newIdentityMatrix :: (WriteMatrix a m, Elem e) => (Int,Int) -> m (a e)
newIdentityMatrix = unsafeConvertIOMatrix . newIdentityIOMatrix
{-# INLINE newIdentityMatrix #-}

-- | Set diagonal elements to one and all other elements to zero.
setIdentityMatrix :: (WriteMatrix a m) => a e -> m ()
setIdentityMatrix a =
    unsafePerformIOWithMatrix a $ setIdentityIOMatrix
{-# INLINE setIdentityMatrix #-}

-- | Get a copy of a matrix.
newCopyMatrix :: (ReadMatrix a m, WriteMatrix b m) =>
    a e -> m (b e)
newCopyMatrix a = unsafeConvertIOMatrix $ 
    newCopyIOMatrix (unsafeMatrixToIOMatrix a)
{-# INLINE newCopyMatrix #-}

-- | Get a copy of a matrix and make sure the returned matrix is not
-- a view.  Specififially, the returned matrix will have @isHermMatrix@
-- equal to @False@.
newCopyMatrix' :: (ReadMatrix a m, WriteMatrix b m) =>
    a e -> m (b e)
newCopyMatrix' a = unsafeConvertIOMatrix $
    newCopyIOMatrix' (unsafeMatrixToIOMatrix a)
{-# INLINE newCopyMatrix' #-}

-- | @copyMatrix dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyMatrix :: (WriteMatrix b m, ReadMatrix a m) => 
    b e -> a e -> m ()
copyMatrix b a = checkBinaryOp (shape b) (shape a) $ unsafeCopyMatrix b a
{-# INLINE copyMatrix #-}

unsafeCopyMatrix :: (WriteMatrix b m, ReadMatrix a m) =>
    b e -> a e -> m ()
unsafeCopyMatrix = liftMatrix2 unsafeCopyVector
{-# INLINE unsafeCopyMatrix #-}

-- | @swapMatrix x y@ swaps the values stored in two matrices.
swapMatrix :: (WriteMatrix a m) =>
    a e -> a e -> m ()
swapMatrix a b = checkBinaryOp (shape b) (shape a) $ unsafeSwapMatrix a b
{-# INLINE swapMatrix #-}

unsafeSwapMatrix :: (WriteMatrix a m) =>
    a e -> a e -> m ()
unsafeSwapMatrix = liftMatrix2 unsafeSwapVector
{-# INLINE unsafeSwapMatrix #-}

-- | Swap the elements in two rows of a matrix.
swapRows :: (WriteMatrix a m) => a e -> Int -> Int -> m ()
swapRows a i j = 
    when (i /= j) $ unsafeSwapVector (rowView a i) (rowView a j)
{-# INLINE swapRows #-}

-- | Swap the elements in two columns of a matrix.
swapCols :: (WriteMatrix a m) => a e -> Int -> Int -> m ()
swapCols a i j = 
    when (i /= j) $ unsafeSwapVector (colView a i) (colView a j)
{-# INLINE swapCols #-}

unsafeSwapRows :: (WriteMatrix a m) => a e -> Int -> Int -> m ()
unsafeSwapRows a i j = 
    when (i /= j) $ unsafeSwapVector (unsafeRowView a i) (unsafeRowView a j)
{-# INLINE unsafeSwapRows #-}

unsafeSwapCols :: (WriteMatrix a m) => a e -> Int -> Int -> m ()
unsafeSwapCols a i j = 
    when (i /= j) $ unsafeSwapVector (unsafeColView a i) (unsafeColView a j)
{-# INLINE unsafeSwapCols #-}

-- | @submatrixView a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrixView :: (BaseMatrix a) => a e -> (Int,Int) -> (Int,Int) -> a e
submatrixView a = checkedSubmatrix (shape a) (unsafeSubmatrixView a)
{-# INLINE submatrixView #-}

-- | Divide the rows of a matrix into two blocks and return views into the
-- blocks.  The integer argument indicates how many rows should be in the
-- first block.
splitRowsAt :: (BaseMatrix a) =>
    Int -> a e -> (a e, a e)
splitRowsAt m1 a = ( submatrixView a (0,0)  (m1,n)
                   , submatrixView a (m1,0) (m2,n)
                   )
  where 
    (m,n) = shape a
    m2    = m - m1
{-# INLINE splitRowsAt #-}

unsafeSplitRowsAt :: (BaseMatrix a) =>
    Int -> a e -> (a e, a e)
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
splitColsAt :: (BaseMatrix a) =>
    Int -> a e -> (a e, a e)
splitColsAt n1 a = ( submatrixView a (0,0)  (m,n1)
                   , submatrixView a (0,n1) (m,n2)
                   )
  where
    (m,n) = shape a
    n2    = n - n1
{-# INLINE splitColsAt #-}

unsafeSplitColsAt :: (BaseMatrix a) =>
    Int -> a e -> (a e, a e)
unsafeSplitColsAt n1 a = ( unsafeSubmatrixView a (0,0)  (m,n1)
                         , unsafeSubmatrixView a (0,n1) (m,n2)
                         )
  where
    (m,n) = shape a
    n2    = n - n1
{-# INLINE unsafeSplitColsAt #-}

-- | Get a list of vector views of the rows of the matrix.
rowViews :: (BaseMatrix a) => a e -> [VectorView a e]
rowViews a = [ unsafeRowView a i | i <- [0..numRows a - 1] ]
{-# INLINE rowViews #-}

-- | Get a list of vector views of the columns of the matrix.
colViews :: (BaseMatrix a) => a e -> [VectorView a e]
colViews a = [ unsafeColView a j | j <- [0..numCols a - 1] ]
{-# INLINE colViews #-}

-- | Get a vector view of the given row in a matrix.
rowView :: (BaseMatrix a) => a e -> Int -> VectorView a e
rowView a = checkedRow (shape a) (unsafeRowView a)
{-# INLINE rowView #-}

unsafeGetRowMatrix :: (ReadMatrix a m, WriteVector y m) =>
    a e -> Int -> m (y e)
unsafeGetRowMatrix a i = newCopyVector (unsafeRowView a i)
{-# INLINE unsafeGetRowMatrix #-}

-- | Get a vector view of the given column in a matrix.
colView :: (BaseMatrix a) => a e -> Int -> VectorView a e
colView a = checkedCol (shape a) (unsafeColView a)
{-# INLINE colView #-}

unsafeGetColMatrix :: (ReadMatrix a m, WriteVector y m) =>
    a e -> Int -> m (y e)
unsafeGetColMatrix a j = newCopyVector (unsafeColView a j)
{-# INLINE unsafeGetColMatrix #-}

-- | Get a vector view of the given diagonal in a matrix.
diagView :: (BaseMatrix a) => a e -> Int -> VectorView a e
diagView a = checkedDiag (shape a) (unsafeDiagView a)
{-# INLINE diagView #-}

-- | Get the given diagonal in a matrix.  Negative indices correspond
-- to sub-diagonals.
getDiag :: (ReadMatrix a m, WriteVector y m) => 
    a e -> Int -> m (y e)
getDiag a = checkedDiag (shape a) (unsafeGetDiag a)
{-# INLINE getDiag #-}

-- Same as 'getDiag' but not range-checked.
unsafeGetDiag :: (ReadMatrix a m, WriteVector y m) => 
    a e -> Int -> m (y e)
unsafeGetDiag a i = newCopyVector (unsafeDiagView a i)
{-# INLINE unsafeGetDiag #-}

-- | Conjugate every element of a matrix.
doConjMatrix :: (WriteMatrix a m, BLAS1 e) => a e -> m ()
doConjMatrix = liftMatrix doConjVector
{-# INLINE doConjMatrix #-}

-- | Get a new matrix equal to the elementwise conjugate of another
-- matrix.
getConjMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    a e -> m (b e)
getConjMatrix = getUnaryMatrixOp doConjMatrix
{-# INLINE getConjMatrix #-}

-- | Scale every element of a matrix by the given value.
scaleByMatrix :: (WriteMatrix a m, BLAS1 e) => e -> a e -> m ()
scaleByMatrix k = liftMatrix (scaleByVector k)
{-# INLINE scaleByMatrix #-}

-- | Get a new matrix by scaling the elements of another matrix
-- by a given value.
getScaledMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    e -> a e -> m (b e)
getScaledMatrix e = getUnaryMatrixOp (scaleByMatrix e)
{-# INLINE getScaledMatrix #-}

-- | Add a constant to every element in a matrix.
shiftByMatrix :: (WriteMatrix a m, BLAS1 e) => e -> a e -> m ()
shiftByMatrix k = liftMatrix (shiftByVector k)
{-# INLINE shiftByMatrix #-}

-- | Get a new matrix by shifting the elements of another matrix
-- by a given value.
getShiftedMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    e -> a e -> m (b e)
getShiftedMatrix e = getUnaryMatrixOp (shiftByMatrix e)
{-# INLINE getShiftedMatrix #-}

-- | Replace the first argument with the elementwise sum.
addMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b e -> a e -> m ()
addMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeAddMatrix b a
{-# INLINE addMatrix #-}

unsafeAddMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b e -> a e -> m ()
unsafeAddMatrix b a = unsafeAxpyMatrix 1 a b
{-# INLINE unsafeAddMatrix #-}

-- | @getAddMatrix a b@ creates a new matrix equal to the sum @a+b@.  The 
-- operands must have the same shape.
getAddMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a e -> b e -> m (c e)
getAddMatrix = checkMatrixOp2 unsafeGetAddMatrix
{-# INLINE getAddMatrix #-}

unsafeGetAddMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a e -> b e -> m (c e)
unsafeGetAddMatrix = unsafeGetBinaryMatrixOp unsafeAddMatrix
{-# INLINE unsafeGetAddMatrix #-}

-- | Replace the first argument with the elementwise difference.
subMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b e -> a e -> m ()    
subMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeSubMatrix b a
{-# INLINE subMatrix #-}

unsafeSubMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b e -> a e -> m ()
unsafeSubMatrix b a = unsafeAxpyMatrix (-1) a b
{-# INLINE unsafeSubMatrix #-}

-- | @getSubMatrix a b@ creates a new matrix equal to the difference @a-b@.  The 
-- operands must have the same shape.
getSubMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a e -> b e -> m (c e)
getSubMatrix = checkMatrixOp2 unsafeGetSubMatrix
{-# INLINE getSubMatrix #-}

unsafeGetSubMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a e -> b e -> m (c e)
unsafeGetSubMatrix = unsafeGetBinaryMatrixOp unsafeSubMatrix
{-# INLINE unsafeGetSubMatrix #-}

-- | @axpyMatrix e a b@ replaces @b := e a + b@.
axpyMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    e -> a e -> b e -> m ()
axpyMatrix alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyMatrix alpha x y
{-# INLINE axpyMatrix #-}

unsafeAxpyMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    e -> a e -> b e -> m ()
unsafeAxpyMatrix alpha = liftMatrix2 (unsafeAxpyVector alpha)
{-# INLINE unsafeAxpyMatrix #-}

-- | Replace the first argument with the elementwise product.
mulMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b e -> a e -> m ()    
mulMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeMulMatrix b a
{-# INLINE mulMatrix #-}

unsafeMulMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b e -> a e -> m ()
unsafeMulMatrix = liftMatrix2 unsafeMulVector
{-# INLINE unsafeMulMatrix #-}

-- | @getMulMatrix a b@ creates a new matrix equal to the elementwise product 
-- @a*b@.  The operands must have the same shape.
getMulMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a e -> b e -> m (c e)
getMulMatrix = checkMatrixOp2 unsafeGetMulMatrix
{-# INLINE getMulMatrix #-}

unsafeGetMulMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a e -> b e -> m (c e)
unsafeGetMulMatrix = unsafeGetBinaryMatrixOp unsafeMulMatrix
{-# INLINE unsafeGetMulMatrix #-}

-- | Replace the first argument with the elementwise quotient.
divMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b e -> a e -> m ()    
divMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeDivMatrix b a
{-# INLINE divMatrix #-}

unsafeDivMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b e -> a e -> m ()
unsafeDivMatrix = liftMatrix2 unsafeDivVector
{-# INLINE unsafeDivMatrix #-}

-- | @getDivMatrix a b@ creates a new matrix equal to the elementwise ratio
-- @a/b@.  The operands must have the same shape.
getDivMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a e -> b e -> m (c e)
getDivMatrix = checkMatrixOp2 unsafeGetDivMatrix
{-# INLINE getDivMatrix #-}

unsafeGetDivMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a e -> b e -> m (c e)
unsafeGetDivMatrix = unsafeGetBinaryMatrixOp unsafeDivMatrix
{-# INLINE unsafeGetDivMatrix #-}

-------------------------------- MMatrix ----------------------------------

-- | A type class for mutable matrices associated with a monad.  The member
-- functions of the type class do not perform any checks on the validity of
-- shapes or indices, so in general their safe counterparts should be
-- preferred.
--
-- Minimal complete definition: 
--    @unsafeGetSApply{Vector,Matrix}@ or @unsafeDoSApplyAdd{Vector,Matrix}@
--
-- Optional:
--    @unsafeDoSApply{Vector,Matrix}_@
--
class (HasHerm a, MonadInterleave m) => MMatrix a m where
    unsafeGetSApplyVector :: (ReadVector x m, WriteVector y m, BLAS2 e) =>
        e -> a e -> x e -> m (y e)
    unsafeGetSApplyVector alpha a x = do
        y  <- newVector_ (numRows a)
        unsafeDoSApplyAddVector alpha a x 0 y
        return y
    {-# INLINE unsafeGetSApplyVector #-}

    unsafeGetSApplyMatrix :: (ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
        e -> a e -> b e -> m (c e)
    unsafeGetSApplyMatrix alpha a b = do
        c  <- newMatrix_ (numRows a, numCols b)
        unsafeDoSApplyAddMatrix alpha a b 0 c
        return c
    {-# INLINE unsafeGetSApplyMatrix #-}

    unsafeDoSApplyAddVector :: (ReadVector x m, WriteVector y m, BLAS2 e) =>
        e -> a e -> x e -> e -> y e -> m ()
    unsafeDoSApplyAddVector alpha a x beta (y :: y e) = do
        (y' :: y e) <- unsafeGetSApplyVector alpha a x
        scaleByVector beta y
        unsafeAxpyVector 1 y' y
    {-# INLINE unsafeDoSApplyAddVector #-}

    unsafeDoSApplyAddMatrix :: (ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
        e -> a e -> b e -> e -> c e -> m ()
    unsafeDoSApplyAddMatrix alpha a b beta c = do
        c' <- unsafeGetSApplyMatrix alpha a b
        scaleByMatrix beta c
        unsafeAxpyMatrix 1 c' c
    {-# INLINE unsafeDoSApplyAddMatrix #-}

    unsafeDoSApplyVector_ :: (WriteVector y m, BLAS2 e) =>
        e -> a e -> y e -> m ()
    unsafeDoSApplyVector_ alpha a (x :: y e) = do
        y <- newVector_ (dim x) :: m (y e)
        unsafeDoSApplyAddVector alpha a x 0 y
        unsafeCopyVector x y
    {-# INLINE unsafeDoSApplyVector_  #-}

    unsafeDoSApplyMatrix_ :: (WriteMatrix b m, BLAS3 e) =>
        e -> a e -> b e -> m ()
    unsafeDoSApplyMatrix_ alpha a (b :: b e) = do
        c <- newMatrix_ (shape b)
        unsafeDoSApplyAddMatrix alpha a b 0 c
        unsafeCopyMatrix b c
    {-# INLINE unsafeDoSApplyMatrix_ #-}

    unsafeGetCol :: (WriteVector x m, Elem e) => a e -> Int -> m (x e)
    
    unsafeGetRow :: (WriteVector x m, Elem e) => a e -> Int -> m (x e)
    unsafeGetRow a i = liftM conj $ unsafeGetCol (herm a) i
    {-# INLINE unsafeGetRow #-}

-- | Get a lazy list of the column vectors in the matrix.
getCols :: (MMatrix a m, WriteVector x m, Elem e)
        => a e -> m [x e]
getCols a =
    let n    = numCols a
        go j | j == n    = return []
             | otherwise = unsafeInterleave $ do
                               c  <- unsafeGetCol a j
                               cs <- go (j+1)
                               return (c:cs)
    in go 0
{-# INLINE getCols #-}

-- | Get a lazy list the row vectors in the matrix.
getRows :: (MMatrix a m, WriteVector x m, Elem e)
        => a e -> m [x e]
getRows a =
    let m    = numRows a
        go i | i == m    = return []
             | otherwise = unsafeInterleave $ do
                                r  <- unsafeGetRow a i
                                rs <- go (i+1)
                                return (r:rs)
    in go 0
{-# INLINE getRows #-}

-- | Rank 1 update to a matrix.  Sets @a := a + alpha x y^H@.
rank1UpdateMatrix :: 
    (WriteMatrix a m, ReadVector x m, ReadVector y m, BLAS2 e) 
    => a e -> e -> x e -> y e -> m ()
rank1UpdateMatrix a alpha x y
    | shape a /= (dim x, dim y) =
        error $ printf ("rank1UpdateMatrix <matrix of shape %s> _"
                       ++ " <vector of dim %d> <vector of dim %d>:"
                       ++ " dimension mismatch.")
                       (show $ shape a) (dim x) (dim y)
    | otherwise =
        unsafeRank1UpdateMatrix a alpha x y
{-# INLINE rank1UpdateMatrix #-}

unsafeRank1UpdateMatrix :: 
    (WriteMatrix a m, ReadVector x m, ReadVector y m, BLAS2 e) 
    => a e -> e -> x e -> y e -> m ()
unsafeRank1UpdateMatrix a alpha x y =
    unsafePerformIOWithMatrix a $ \a' ->
        unsafeRank1UpdateIOMatrix a' alpha (unsafeVectorToIOVector x)
                                           (unsafeVectorToIOVector y)
{-# INLINE unsafeRank1UpdateMatrix #-}

-- | @gemv alpha a x beta y@ replaces @y := alpha a * x + beta y@.
gemv :: (ReadMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> a e -> x e -> e -> y e -> m ()
gemv alpha a x beta y
    | n == 0 =
        scaleByVector beta y
    | otherwise =
        let (m',n') = if isHermMatrix a then (n,m) else (m,n)
            ldA     = ldaMatrix a
            incX    = stride x
            incY    = stride y
        in unsafePerformIOWithVector y $ \y' ->
           withIOMatrix (unsafeMatrixToIOMatrix a) $ \pA ->
           withIOVector (unsafeVectorToIOVector x) $ \pX ->
           withIOVector y' $ \pY ->
               BLAS.gemv (transEnumMatrix a) (conjEnum x) (conjEnum y) m' n' 
                         alpha pA ldA pX incX beta pY incY
  where
    (m,n) = shape a
{-# INLINE gemv #-}

-- | @gemm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gemm :: (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
    e -> a e -> b e -> e -> c e -> m ()
gemm alpha a b beta c
    | numRows a == 0 || numCols a == 0 || numCols b == 0 = 
        scaleByMatrix beta c
    | isHermMatrix c = gemm (conjugate alpha) (herm b) (herm a) (conjugate beta) (herm c)
    | otherwise =
        let transA = transEnumMatrix a
            transB = transEnumMatrix b
            (m,n)  = shape c
            k      = numCols a
            ldA    = ldaMatrix a
            ldB    = ldaMatrix b
            ldC    = ldaMatrix c
        in 
            unsafePerformIOWithMatrix c $ \c' ->
            withIOMatrix c' $ \pC ->
            withIOMatrix (unsafeMatrixToIOMatrix a) $ \pA ->
            withIOMatrix (unsafeMatrixToIOMatrix b) $ \pB ->
                BLAS.gemm transA transB m n k alpha pA ldA pB ldB beta pC ldC
{-# INLINE gemm #-}

unsafeGetColHermMatrix :: (ReadMatrix a m, WriteVector x m, Elem e)
                       => Herm a e -> Int -> m (x e)
unsafeGetColHermMatrix (Herm l a) i =
    let n  = numRows a
        r  = unsafeRowView a i
        c  = unsafeColView a i 
        (r',c') = case l of { Lower -> (c,r) ; Upper -> (conj r, conj c) }
    in do
        x <- newVector_ n
        unsafeCopyVector (subvectorView x 0 i)     (conj $ subvectorView c' 0 i)
        unsafeCopyVector (subvectorView x i (n-i)) (subvectorView r' i (n-i))
        return x
{-# INLINE unsafeGetColHermMatrix #-}

hemv :: (ReadMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> Herm a e -> x e -> e -> y e -> m ()
hemv alpha h x beta y =
    let (u,a) = hermToBase h
        n     = numCols a
        uplo  = case isHermMatrix a of
                    True  -> flipUpLo u
                    False -> u
        ldA   = ldaMatrix a
        incX  = stride x
        incY  = stride y
    in unsafePerformIOWithVector y $ \y' ->
       withIOMatrix (unsafeMatrixToIOMatrix a) $ \pA ->
       withIOVector (unsafeVectorToIOVector x) $ \pX ->
       withIOVector y' $ \pY ->
           BLAS.hemv uplo (conjEnum x) (conjEnum y) n
                     alpha pA ldA pX incX beta pY incY
{-# INLINE hemv #-}

hemm :: (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
    e -> Herm a e -> b e -> e -> c e -> m ()
hemm alpha h b beta c
    | numRows b == 0 || numCols b == 0 || numCols c == 0 = return ()
    | (isHermMatrix a) /= (isHermMatrix c) || (isHermMatrix a) /= (isHermMatrix b) =
        zipWithM_ (\x y -> hemv alpha h x beta y) (colViews b) (colViews c)
    | otherwise =
        let (m,n)   = shape c
            (side,u',m',n')
                    = if isHermMatrix a
                          then (RightSide, flipUpLo u, n, m)
                          else (LeftSide,  u,          m, n)
            uploA   = u'
            (alpha',beta') = if isHermMatrix c then (conjugate alpha, conjugate beta) else (alpha,beta)
            ldA     = ldaMatrix a
            ldB     = ldaMatrix b
            ldC     = ldaMatrix c
        in 
            withMatrixPtr   a $ \pA ->
            withIOMatrix (unsafeMatrixToIOMatrix b) $ \pB ->
            withIOMatrix (unsafeMatrixToIOMatrix c) $ \pC ->
                BLAS.hemm side uploA m' n' alpha' pA ldA pB ldB beta' pC ldC
    where
      withMatrixPtr d f = unsafePerformIOWithMatrix d $ flip withIOMatrix f
      (u,a) = hermToBase h
{-# INLINE hemm #-}

unsafeDoSApplyAddVectorTriMatrix :: (ReadMatrix a m,
    ReadVector x m, WriteVector y m, BLAS2 e) =>
        e -> Tri a e -> x e -> e -> y e -> m ()
unsafeDoSApplyAddVectorTriMatrix alpha t x beta (y :: y e) =
    if beta == 0
        then unsafeDoSApplyTriMatrix alpha t x y
        else do
            (y' :: y e) <- newCopyVector y
            unsafeDoSApplyTriMatrix alpha t x y'
            scaleByVector beta y
            unsafeAxpyVector 1 y' y
{-# INLINE unsafeDoSApplyAddVectorTriMatrix #-}

unsafeDoSApplyAddMatrixTriMatrix :: (ReadMatrix a m,
    ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
        e -> Tri a e -> b e -> e -> c e -> m ()
unsafeDoSApplyAddMatrixTriMatrix alpha t b beta (c :: c e) =
    if beta == 0
        then unsafeDoSApplyMatrixTriMatrix alpha t b c
        else do
            (c' :: c e) <- newCopyMatrix c
            unsafeDoSApplyMatrixTriMatrix alpha t b c'
            scaleByMatrix beta c
            unsafeAxpyMatrix 1 c' c
{-# INLINE unsafeDoSApplyAddMatrixTriMatrix #-}

unsafeDoSApplyTriMatrix :: (ReadMatrix a m,
    ReadVector x m, WriteVector y m, BLAS2 e) =>
        e -> Tri a e -> x e -> y e -> m ()
unsafeDoSApplyTriMatrix alpha t x y =
    let mn      = min (numRows t) (numCols t)
        (x1,x2) = splitElemsAt mn x
        (y1,y2) = splitElemsAt mn y
        (u,d,a) = triToBase t
    in do
        unsafeCopyVector y1 x1
        case (u, toLower d a, toUpper d a) of
            (Lower,Left t',_) -> do
                trmv alpha t' y1
            (Lower,Right (t',r),_) -> do
                trmv alpha t'    y1
                gemv alpha r x 0 y2
            (Upper,_,Left t') -> do
                trmv alpha t' y1
                setZeroVector y2
            (Upper,_,Right (t',r)) -> do
                trmv alpha t'     y1
                gemv alpha r x2 1 y1
{-# INLINE unsafeDoSApplyTriMatrix #-}

unsafeDoSApplyMatrixTriMatrix :: (ReadMatrix a m,
    ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
        e -> Tri a e -> b e -> c e -> m ()
unsafeDoSApplyMatrixTriMatrix alpha t b c =
    let mn      = min (numRows t) (numCols t)
        (b1,b2) = splitRowsAt mn b
        (c1,c2) = splitRowsAt mn c
        (u,d,a) = triToBase t
    in do
        unsafeCopyMatrix c1 b1
        case (u, toLower d a, toUpper d a) of
            (Lower,Left t',_) -> do
                trmm alpha t' c1
            (Lower,Right (t',r),_) -> do
                trmm alpha t'    c1
                gemm alpha r b 0 c2
            (Upper,_,Left t') -> do
                trmm alpha t' c1
                setZeroMatrix c2
            (Upper,_,Right (t',r)) -> do
                trmm alpha t'     c1
                gemm alpha r b2 1 c1
{-# INLINE unsafeDoSApplyMatrixTriMatrix #-}

toLower :: (BaseMatrix a) => DiagEnum -> a e 
        -> Either (Tri a e) 
                  (Tri a e, a e)
toLower d a =
    if m <= n
        then Left $  triFromBase Lower d (unsafeSubmatrixView a (0,0) (m,m))
        else let t = triFromBase Lower d (unsafeSubmatrixView a (0,0) (n,n))
                 r = unsafeSubmatrixView a (n,0) (k,n)
             in Right (t,r)
  where
    (m,n) = shape a
    k     = m - n
{-# INLINE toLower #-}
    
toUpper :: (BaseMatrix a) => DiagEnum -> a e
        -> Either (Tri a e)
                  (Tri a e, a e)
toUpper d a =
    if n <= m
        then Left $  triFromBase Upper d (unsafeSubmatrixView a (0,0) (n,n))
        else let t = triFromBase Upper d (unsafeSubmatrixView a (0,0) (m,m))
                 r = unsafeSubmatrixView a (0,m) (m,k)
             in Right (t,r)
  where
    (m,n) = shape a
    k     = n - m
{-# INLINE toUpper #-}

unsafeGetColTriMatrix :: (ReadMatrix a m, WriteVector x m, Elem e)
                       => Tri a e -> Int -> m (x e)
unsafeGetColTriMatrix a@(Tri Upper _ _) j = 
    liftM conj $ unsafeGetRowTriMatrix (herm a) j

unsafeGetColTriMatrix (Tri Lower _ a) j | j >= numRows a =
    newZeroVector (numRows a)

unsafeGetColTriMatrix (Tri Lower Unit a) j = do
    x <- newVector_ m
    setZeroVector (unsafeSubvectorView x 0 j)
    unsafeWriteElem x j 1
    unsafeCopyVector (unsafeSubvectorView x (j+1) (m-j-1))
                     (unsafeSubvectorView c (j+1) (m-j-1))
    return x
  where
    m = numRows a
    c = unsafeColView a j
    
unsafeGetColTriMatrix (Tri Lower NonUnit a) j = do
    x <- newVector_ m
    setZeroVector (unsafeSubvectorView x 0 j)
    unsafeCopyVector (unsafeSubvectorView x j (m-j))
                     (unsafeSubvectorView c j (m-j))
    return x
  where
    m = numRows a
    c = unsafeColView a j
{-# INLINE unsafeGetColTriMatrix #-}

unsafeGetRowTriMatrix :: (ReadMatrix a m, WriteVector x m, Elem e)
                      => Tri a e -> Int -> m (x e)
unsafeGetRowTriMatrix a@(Tri Upper _ _) i = 
    liftM conj $ unsafeGetColTriMatrix (herm a) i
    
unsafeGetRowTriMatrix (Tri Lower Unit a) i = do
    x <- newVector_ n
    unsafeCopyVector (unsafeSubvectorView x 0 $ i)
                     (unsafeSubvectorView r 0 $ i)
    when (i < n) $ do
        unsafeWriteElem x i 1
        setZeroVector (unsafeSubvectorView x (i+1) (n-i-1))
    return x
  where
    n = numCols a
    r = unsafeRowView a i

unsafeGetRowTriMatrix (Tri Lower NonUnit a) i = do
    x <- newVector_ n
    unsafeCopyVector (unsafeSubvectorView x 0 $ i+1)
                     (unsafeSubvectorView r 0 $ i+1)
    when (i < n) $ do
        setZeroVector (unsafeSubvectorView x (i+1) (n-i-1))
    return x
  where
    n = numCols a
    r = unsafeRowView a i
{-# INLINE unsafeGetRowTriMatrix #-}

trmv :: (ReadMatrix a m, WriteVector y m, BLAS2 e) =>
    e -> Tri a e -> y e -> m ()
trmv alpha t x =
    let (u,d,a)      = triToBase t
        (trans,uplo) = if isHermMatrix a then (ConjTrans, flipUpLo u) 
                                         else (NoTrans  , u)
        n            = dim x
        ldA          = ldaMatrix a
        incX         = stride x
        in unsafePerformIOWithVector x $ \x' ->
           withIOMatrix (unsafeMatrixToIOMatrix a) $ \pA ->
           withIOVector x' $ \pX -> do
               when (alpha /= 1) $ scaleByVector alpha x'
               BLAS.trmv uplo trans d (conjEnum x) n pA ldA pX incX
{-# INLINE trmv #-}

trmm :: (ReadMatrix a m, WriteMatrix b m, BLAS3 e) =>
    e -> Tri a e -> b e -> m ()
trmm _ _ b
    | numRows b == 0 || numCols b == 0 = return ()
trmm alpha t b =
    let (u,d,a)   = triToBase t
        (h,u')    = if isHermMatrix a then (ConjTrans, flipUpLo u) else (NoTrans, u)
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if isHermMatrix b
                        then (RightSide, flipTrans h, n, m, conjugate alpha)
                        else (LeftSide , h          , m, n, alpha       )
        uploA     = u'
        transA    = h'
        diagA     = d
        ldA       = ldaMatrix a
        ldB       = ldaMatrix b
    in  
        withMatrixPtr   a $ \pA ->
        withIOMatrix (unsafeMatrixToIOMatrix b) $ \pB ->
            BLAS.trmm side uploA transA diagA m' n' alpha' pA ldA pB ldB
  where 
    withMatrixPtr d f = unsafePerformIOWithMatrix d $ flip withIOMatrix f
{-# INLINE trmm #-}
  
unsafeDoSSolveTriMatrix :: (ReadMatrix a m,
    ReadVector y m, WriteVector x m, BLAS2 e) =>
        e -> Tri a e -> y e -> x e -> m ()
unsafeDoSSolveTriMatrix alpha t y x =
    let mn      = min (numRows t) (numCols t)
        (x1,x2) = splitElemsAt mn x
        (y1, _) = splitElemsAt mn y
        (u,d,a) = triToBase t
    in do
        unsafeCopyVector x1 y1
        case (u, toLower d a, toUpper d a) of
            (Lower,Left t',_) -> do
                trsv alpha t' x1
                setZeroVector x2
            (Lower,Right (t',_),_) -> do
                trsv alpha t' x1
            (Upper,_,Left t') -> do
                trsv alpha t' x1
            (Upper,_,Right (t',_)) -> do
                trsv alpha t' x1
                setZeroVector x2
{-# INLINE unsafeDoSSolveTriMatrix #-}

unsafeDoSSolveMatrixTriMatrix :: (ReadMatrix a m,
    ReadMatrix c m, WriteMatrix b m, BLAS3 e) =>
        e -> Tri a e -> c e -> b e -> m ()
unsafeDoSSolveMatrixTriMatrix alpha t c b =
    let mn      = min (numRows t) (numCols t)
        (b1,b2) = splitRowsAt mn b
        (c1, _) = splitRowsAt mn c
        (u,d,a) = triToBase t
    in do
        unsafeCopyMatrix b1 c1
        case (u, toLower d a, toUpper d a) of
            (Lower,Left t',_) -> do
                trsm alpha t' b1
                setZeroMatrix b2
            (Lower,Right (t',_),_) -> do
                trsm alpha t' b1
            (Upper,_,Left t') -> do
                trsm alpha t' b1
            (Upper,_,Right (t',_)) -> do
                trsm alpha t' b1
                setZeroMatrix b2
{-# INLINE unsafeDoSSolveMatrixTriMatrix #-}

trsv :: (ReadMatrix a m, WriteVector y m, BLAS2 e) =>
    e -> Tri a e -> y e -> m ()
trsv alpha t x =
    let (u,d,a)      = triToBase t
        (trans,uplo) = if isHermMatrix a then (ConjTrans, flipUpLo u) 
                                         else (NoTrans  , u)
        n            = dim x
        ldA          = ldaMatrix a
        incX         = stride x
        in unsafePerformIOWithVector x $ \x' ->
           withIOMatrix (unsafeMatrixToIOMatrix a) $ \pA ->
           withIOVector x' $ \pX -> do
               when (alpha /= 1) $ scaleByVector alpha x'
               BLAS.trsv uplo trans d (conjEnum x) n pA ldA pX incX
{-# INLINE trsv #-}

trsm :: (ReadMatrix a m, WriteMatrix b m, BLAS3 e) =>
    e -> Tri a e -> b e -> m ()
trsm _ _ b
    | numRows b == 0 || numCols b == 0 = return ()
trsm alpha t b =
    let (u,d,a)   = triToBase t
        (h,u')    = if isHermMatrix a then (ConjTrans, flipUpLo u) else (NoTrans, u)
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if isHermMatrix b
                        then (RightSide, flipTrans h, n, m, conjugate alpha)
                        else (LeftSide , h          , m, n, alpha     )
        uploA     = u'
        transA    = h'
        diagA     = d
        ldA       = ldaMatrix a
        ldB       = ldaMatrix b
    in 
        withMatrixPtr   a $ \pA ->
        withIOMatrix (unsafeMatrixToIOMatrix b) $ \pB -> do
            BLAS.trsm side uploA transA diagA m' n' alpha' pA ldA pB ldB
  where 
    withMatrixPtr d f = unsafePerformIOWithMatrix d $ flip withIOMatrix f
{-# INLINE trsm #-}

------------------------------------ MSolve ------------------------------

-- | A type class for mutable matrices with inverses.  The member
-- functions of the type class do not perform any checks on the validity
-- of shapes or indices, so in general their safe counterparts should be
-- preferred.
class (MatrixShaped a, Monad m) => MSolve a m where
    unsafeDoSolveVector :: (ReadVector y m, WriteVector x m, BLAS2 e) =>
        a e -> y e -> x e -> m ()
    unsafeDoSolveVector = unsafeDoSSolveVector 1
    {-# INLINE unsafeDoSolveVector #-}
    
    unsafeDoSolveMatrix :: (ReadMatrix c m, WriteMatrix b m, BLAS3 e) =>
        a e -> c e -> b e -> m ()
    unsafeDoSolveMatrix = unsafeDoSSolveMatrix 1
    {-# INLINE unsafeDoSolveMatrix #-}    
    
    unsafeDoSSolveVector :: (ReadVector y m, WriteVector x m, BLAS2 e) =>
        e -> a e -> y e -> x e -> m ()
    unsafeDoSSolveVector alpha a y x = do
        unsafeDoSolveVector a y x
        scaleByVector alpha x
    {-# INLINE unsafeDoSSolveVector #-}        
    
    unsafeDoSSolveMatrix :: (ReadMatrix c m, WriteMatrix b m, BLAS3 e) =>
        e -> a e -> c e -> b e -> m ()
    unsafeDoSSolveMatrix alpha a c b = do
        unsafeDoSolveMatrix a c b
        scaleByMatrix alpha b
    {-# INLINE unsafeDoSSolveMatrix #-}

    unsafeDoSolveVector_ :: (WriteVector x m, BLAS2 e) => a e -> x e -> m ()
    unsafeDoSolveVector_ = unsafeDoSSolveVector_ 1
    {-# INLINE unsafeDoSolveVector_ #-}

    unsafeDoSSolveVector_ :: (WriteVector x m, BLAS2 e) => e -> a e -> x e -> m ()
    unsafeDoSSolveVector_ alpha a x = do
        scaleByVector alpha x
        unsafeDoSolveVector_ a x
    {-# INLINE unsafeDoSSolveVector_ #-}        
        
    unsafeDoSolveMatrix_ :: (WriteMatrix b m, BLAS3 e) => a e -> b e -> m ()
    unsafeDoSolveMatrix_ = unsafeDoSSolveMatrix_ 1
    {-# INLINE unsafeDoSolveMatrix_ #-}
        
    unsafeDoSSolveMatrix_ :: (WriteMatrix b m, BLAS3 e) => e -> a e -> b e -> m ()
    unsafeDoSSolveMatrix_ alpha a b = do
        scaleByMatrix alpha b
        unsafeDoSolveMatrix_ a b
    {-# INLINE unsafeDoSSolveMatrix_ #-}

------------------------------------ Instances ------------------------------


instance BaseMatrix IOMatrix where
    ldaMatrix = ldaIOMatrix
    {-# INLINE ldaMatrix #-}
    transEnumMatrix = transEnumIOMatrix
    {-# INLINE transEnumMatrix #-}
    unsafeSubmatrixView = unsafeSubmatrixViewIOMatrix
    {-# INLINE unsafeSubmatrixView #-}
    unsafeDiagView = unsafeDiagViewIOMatrix
    {-# INLINE unsafeDiagView #-}
    unsafeRowView = unsafeRowViewIOMatrix
    {-# INLINE unsafeRowView #-}
    unsafeColView = unsafeColViewIOMatrix
    {-# INLINE unsafeColView #-}
    maybeViewMatrixAsVector = maybeViewIOMatrixAsVector
    {-# INLINE maybeViewMatrixAsVector #-}
    maybeViewVectorAsMatrix = maybeViewVectorAsIOMatrix
    {-# INLINE maybeViewVectorAsMatrix #-}
    maybeViewVectorAsRow = maybeViewVectorAsRowIOMatrix
    {-# INLINE maybeViewVectorAsRow #-}    
    maybeViewVectorAsCol = maybeViewVectorAsColIOMatrix
    {-# INLINE maybeViewVectorAsCol #-}
    unsafeIOMatrixToMatrix = id
    {-# INLINE unsafeIOMatrixToMatrix #-}
    unsafeMatrixToIOMatrix = id
    {-# INLINE unsafeMatrixToIOMatrix #-}

instance MMatrix IOMatrix IO where
    unsafeDoSApplyAddVector = gemv
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = gemm
    {-# INLINE unsafeDoSApplyAddMatrix #-}
    unsafeGetRow = unsafeGetRowMatrix
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColMatrix
    {-# INLINE unsafeGetCol #-}

instance MMatrix (Herm IOMatrix) IO where
    unsafeDoSApplyAddVector = hemv
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = hemm
    {-# INLINE unsafeDoSApplyAddMatrix #-}    
    unsafeGetCol = unsafeGetColHermMatrix
    {-# INLINE unsafeGetCol #-}
    
instance MMatrix (Tri IOMatrix) IO where
    unsafeDoSApplyAddVector = unsafeDoSApplyAddVectorTriMatrix
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = unsafeDoSApplyAddMatrixTriMatrix
    {-# INLINE unsafeDoSApplyAddMatrix #-}
    unsafeDoSApplyVector_  = trmv
    {-# INLINE unsafeDoSApplyVector_  #-}
    unsafeDoSApplyMatrix_ = trmm
    {-# INLINE unsafeDoSApplyMatrix_ #-}
    unsafeGetCol = unsafeGetColTriMatrix
    {-# INLINE unsafeGetCol #-}
    unsafeGetRow = unsafeGetRowTriMatrix
    {-# INLINE unsafeGetRow #-}

instance MSolve  (Tri IOMatrix) IO where
    unsafeDoSSolveVector = unsafeDoSSolveTriMatrix
    {-# INLINE unsafeDoSSolveVector #-}
    unsafeDoSSolveMatrix = unsafeDoSSolveMatrixTriMatrix
    {-# INLINE unsafeDoSSolveMatrix #-}    
    unsafeDoSSolveVector_ = trsv
    {-# INLINE unsafeDoSSolveVector_ #-}
    unsafeDoSSolveMatrix_ = trsm
    {-# INLINE unsafeDoSSolveMatrix_ #-}

instance ReadMatrix IOMatrix IO where
    freezeMatrix = freezeIOMatrix
    {-# INLINE freezeMatrix #-}
    unsafeFreezeMatrix = unsafeFreezeIOMatrix
    {-# INLINE unsafeFreezeMatrix #-}
    unsafePerformIOWithMatrix a f = f a
    {-# INLINE unsafePerformIOWithMatrix #-}

instance WriteMatrix IOMatrix IO where
    newMatrix_ = newIOMatrix_
    {-# INLINE newMatrix_ #-}
    thawMatrix = thawIOMatrix
    {-# INLINE thawMatrix #-}
    unsafeThawMatrix = unsafeThawIOMatrix
    {-# INLINE unsafeThawMatrix #-}
    unsafeConvertIOMatrix = id
    {-# INLINE unsafeConvertIOMatrix #-}

-- | A type class for immutable matrices.  The member functions of the
-- type class do not perform any checks on the validity of shapes or
-- indices, so in general their safe counterparts should be preferred.
class (HasHerm a) => IMatrix a where
    unsafeSApplyVector :: (BLAS2 e) => e -> a e -> Vector e -> Vector e
    unsafeSApplyMatrix :: (BLAS3 e) => e -> a e -> Matrix e -> Matrix e

    unsafeRow :: (Elem e) => a e -> Int -> Vector e
    unsafeRow a = conj . unsafeCol (herm a)
    {-# INLINE unsafeRow #-}

    unsafeCol :: (Elem e) => a e -> Int -> Vector e

-- | Get a list the row vectors in the matrix.
rows :: (IMatrix a, Elem e) => a e -> [Vector e]
rows a = [ unsafeRow a i | i <- [0..numRows a - 1] ]
{-# INLINE rows #-}

-- | Get a list the column vectors in the matrix.
cols :: (IMatrix a, Elem e) => a e -> [Vector e]
cols a = [ unsafeCol a j | j <- [0..numCols a - 1] ]
{-# INLINE cols #-}

-- | Create a new matrix of the given size and initialize the given elements to
-- the given values.  All other elements get set to zero.
matrix :: (Elem e) => (Int,Int) -> [((Int,Int), e)] -> Matrix e
matrix mn ies = unsafePerformIO $
    unsafeFreezeIOMatrix =<< newIOMatrix mn ies
{-# INLINE matrix #-}

-- Same as 'matrix' but does not do any bounds checking.
unsafeMatrix :: (Elem e) => (Int,Int) -> [((Int,Int), e)] -> Matrix e
unsafeMatrix mn ies =  unsafePerformIO $ 
    unsafeFreezeIOMatrix =<< unsafeNewIOMatrix mn ies
{-# INLINE unsafeMatrix #-}

-- | Create a new matrix with the given elements in row-major order.
listMatrix :: (Elem e) => (Int,Int) -> [e] -> Matrix e
listMatrix mn es = unsafePerformIO $ 
    unsafeFreezeIOMatrix =<< newListIOMatrix mn es
{-# INLINE listMatrix #-}

replaceMatrix :: Matrix e -> [((Int,Int),e)] -> Matrix e
replaceMatrix (Matrix a@(IOMatrix _ m n _ _ _)) ies =
    let go b (((i,j),e):ies') = do
            when (i < 0 || i >= m || j < 0 || j >= n) $ error $ printf
                ("(//) <matrix of shape (%d,%d)> [ ..., ((%d,%d),_), ... ]:"
                ++ "invalid index") m n i j 
            unsafeWriteElemIOMatrix b (i,j) e
            go b ies'
        go _ [] = return ()
    in
    unsafePerformIO $ do
        b  <- newCopyIOMatrix a
        go b ies
        return (Matrix b)
{-# INLINE replaceMatrix #-}

unsafeReplaceMatrix :: Matrix e -> [((Int,Int),e)] -> Matrix e
unsafeReplaceMatrix (Matrix a@(IOMatrix _ _ _ _ _ _)) ies =
    let go b (((i,j),e):ies') = do
            unsafeWriteElemIOMatrix b (i,j) e
            go b ies'
        go _ [] = return ()
    in
    unsafePerformIO $ do
        b  <- newCopyIOMatrix a
        go b ies
        return (Matrix b)
{-# INLINE unsafeReplaceMatrix #-}

-- | Create a matrix of the given shape from a list of rows
rowsMatrix :: (Elem e) => (Int,Int) -> [Vector e] -> Matrix e
rowsMatrix mn rs = unsafePerformIO $ 
    unsafeFreezeIOMatrix =<< newRowsIOMatrix mn rs
{-# INLINE rowsMatrix #-}

-- | Create a matrix of the given shape from a list of columns
colsMatrix :: (Elem e) => (Int,Int) -> [Vector e] -> Matrix e
colsMatrix mn cs = unsafePerformIO $ 
    unsafeFreezeIOMatrix =<< newColsIOMatrix mn cs
{-# INLINE colsMatrix #-}

-- | Get a matrix from a row vector.
matrixFromRow :: Vector e -> Matrix e
matrixFromRow (Vector x@(IOVector _ _ _ _ _)) = 
    case maybeViewVectorAsRow x of
        Just x' -> Matrix x'
        Nothing -> unsafePerformIO $ unsafeFreezeIOMatrix =<< newRowIOMatrix x
{-# INLINE matrixFromRow #-}

-- | Get a matrix from a column vector.
matrixFromCol :: Vector e -> Matrix e
matrixFromCol (Vector x@(IOVector _ _ _ _ _)) = 
    case maybeViewVectorAsCol x of
        Just x' -> Matrix x'
        Nothing -> unsafePerformIO $ unsafeFreezeIOMatrix =<< newColIOMatrix x
{-# INLINE matrixFromCol #-}

-- | Get a matrix from the elements stored in columnwise order in the vector.
matrixFromVector :: (Int,Int) -> Vector e -> Matrix e
matrixFromVector (m,n) x@(Vector (IOVector _ _ _ _ _))
    | dim x /= m*n =
        error $ "matrixFromVector " ++ show (m,n) ++ "<vector of dim "
              ++ show (dim x) ++ ">: vector dimension must be equal to "
              ++ "the number of elements in the desired matrix"
    | otherwise =
        case maybeViewVectorAsMatrix (m,n) x of
            Just a  -> a
            Nothing -> listMatrix (m,n) (elems x)
{-# INLINE matrixFromVector #-}

-- | Get a vector by concatenating the columns of the matrix.
vectorFromMatrix :: Matrix e -> Vector e
vectorFromMatrix a@(Matrix (IOMatrix _ _ _ _ _ _)) =
    case maybeViewMatrixAsVector a of
        Just x  -> x
        Nothing -> listVector (size a) (concatMap elems (colViews a))
{-# INLINE vectorFromMatrix #-}

-- | Get a new zero of the given shape.
zeroMatrix :: (Elem e) => (Int,Int) -> Matrix e
zeroMatrix mn = unsafePerformIO $
    unsafeFreezeIOMatrix =<< newZeroIOMatrix mn
{-# INLINE zeroMatrix #-}

-- | Get a new constant of the given shape.
constantMatrix :: (Elem e) => (Int,Int) -> e -> Matrix e
constantMatrix mn e = unsafePerformIO $
    unsafeFreezeIOMatrix =<< newConstantIOMatrix mn e
{-# INLINE constantMatrix #-}

-- | Get a new matrix of the given shape with ones along the diagonal and
-- zeroes everywhere else.
identityMatrix :: (Elem e) => (Int,Int) -> Matrix e
identityMatrix mn = unsafePerformIO $
    unsafeFreezeIOMatrix =<< newIdentityIOMatrix mn
{-# INLINE identityMatrix #-}

-- | @submatrix a ij mn@ returns the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrix :: Matrix e -> (Int,Int) -> (Int,Int) -> Matrix e
submatrix (Matrix a) ij mn = 
    Matrix $ submatrixView a ij mn
{-# INLINE submatrix #-}

unsafeSubmatrix :: Matrix e -> (Int,Int) -> (Int,Int) -> Matrix e
unsafeSubmatrix (Matrix a) ij mn = 
    Matrix $ unsafeSubmatrixView a ij mn
{-# INLINE unsafeSubmatrix #-}

-- | Get a the given diagonal in a matrix.  Negative indices correspond to
-- sub-diagonals.
diag :: Matrix e -> Int -> Vector e
diag (Matrix a) i = Vector (diagView a i)
{-# INLINE diag #-}

-- Same as 'diag' but index is not range-checked.
unsafeDiag :: Matrix e -> Int -> Vector e
unsafeDiag (Matrix a) i = Vector (diagView a i)
{-# INLINE unsafeDiag #-}

unsafeAtMatrix :: Matrix e -> (Int,Int) -> e
unsafeAtMatrix (Matrix a) = inlinePerformIO . unsafeReadElemIOMatrix a
{-# INLINE unsafeAtMatrix #-}

indicesMatrix :: Matrix e -> [(Int,Int)]
indicesMatrix (Matrix a) = indicesIOMatrix a
{-# INLINE indicesMatrix #-}

elemsMatrix :: Matrix e -> [e]
elemsMatrix (Matrix a) = 
    case maybeViewIOMatrixAsVector a of
        (Just x) -> elemsVector (Vector x)
        Nothing  -> concatMap (elemsVector . Vector) (vecViews a)
  where
    vecViews = if isHermIOMatrix a
                   then rowViews
                   else colViews
{-# INLINE elemsMatrix #-}

assocsMatrix :: Matrix e -> [((Int,Int),e)]
assocsMatrix a = zip (indicesMatrix a) (elemsMatrix a)
{-# INLINE assocsMatrix #-}

tmapMatrix :: (e -> e) -> Matrix e -> Matrix e
tmapMatrix f a@(Matrix ma@(IOMatrix _ _ _ _ _ _))
    | isHermIOMatrix ma = herm $ 
                              listMatrix (n,m) $ map (conjugate . f) (elems a)
    | otherwise         = listMatrix (m,n) $ map f (elems a)
                              
  where
    (m,n) = shape a
{-# INLINE tmapMatrix #-}

tzipWithMatrix :: (e -> e -> e) -> Matrix e -> Matrix e -> Matrix e
tzipWithMatrix f a@(Matrix (IOMatrix _ _ _ _ _ _)) b
    | shape b /= mn =
        error ("tzipWith: matrix shapes differ; first has shape `" ++
                show mn ++ "' and second has shape `" ++
                show (shape b) ++ "'")
    | otherwise =
        listMatrix mn $ zipWith f (colElems a) (colElems b)
  where
    mn = shape a
    colElems = (concatMap elems) . colViews
{-# INLINE tzipWithMatrix #-}

instance Shaped Matrix (Int,Int) where
    shape (Matrix a) = shapeIOMatrix a
    {-# INLINE shape #-}
    bounds (Matrix a) = boundsIOMatrix a
    {-# INLINE bounds #-}

instance MatrixShaped Matrix where

instance HasHerm Matrix where
    herm (Matrix a) = Matrix (herm a)
    {-# INLINE herm #-}
    
instance HasVectorView Matrix where
    type VectorView Matrix = Vector
    
instance BaseMatrix Matrix where
    ldaMatrix (Matrix a) = ldaIOMatrix a
    {-# INLINE ldaMatrix #-}
    transEnumMatrix (Matrix a) = transEnumIOMatrix a
    {-# INLINE transEnumMatrix #-}
    unsafeSubmatrixView (Matrix a) ij mn =
        Matrix (unsafeSubmatrixViewIOMatrix a ij mn)
    {-# INLINE unsafeSubmatrixView #-}
    unsafeDiagView (Matrix a) i = Vector (unsafeDiagViewIOMatrix a i)
    {-# INLINE unsafeDiagView #-}
    unsafeRowView (Matrix a) i = Vector (unsafeRowViewIOMatrix a i)
    {-# INLINE unsafeRowView #-}
    unsafeColView (Matrix a) i = Vector (unsafeColViewIOMatrix a i)
    {-# INLINE unsafeColView #-}
    maybeViewMatrixAsVector (Matrix a) = liftM Vector (maybeViewMatrixAsVector a)
    {-# INLINE maybeViewMatrixAsVector #-}
    maybeViewVectorAsMatrix mn (Vector x) = 
        liftM Matrix $ maybeViewVectorAsIOMatrix mn x
    {-# INLINE maybeViewVectorAsMatrix #-}
    maybeViewVectorAsRow (Vector x) = liftM Matrix (maybeViewVectorAsRow x)
    {-# INLINE maybeViewVectorAsRow #-}
    maybeViewVectorAsCol (Vector x) = liftM Matrix (maybeViewVectorAsCol x)
    {-# INLINE maybeViewVectorAsCol #-}
    unsafeIOMatrixToMatrix = Matrix
    {-# INLINE unsafeIOMatrixToMatrix #-}
    unsafeMatrixToIOMatrix (Matrix a) = a
    {-# INLINE unsafeMatrixToIOMatrix #-}

-- The NOINLINE pragmas and the strictness annotations here are *really*
-- important.  Otherwise, the compiler might think that certain actions
-- don't need to be run.
instance (MonadInterleave m) => ReadMatrix Matrix m where
    freezeMatrix (Matrix a) = return $! unsafePerformIO $ freezeIOMatrix a
    {-# NOINLINE freezeMatrix #-}
    unsafeFreezeMatrix = return . id
    {-# INLINE unsafeFreezeMatrix #-}
    unsafePerformIOWithMatrix (Matrix a) f = return $! unsafePerformIO $ f a
    {-# NOINLINE unsafePerformIOWithMatrix #-}

instance ITensor Matrix (Int,Int) where
    size (Matrix a) = sizeIOMatrix a
    {-# INLINE size #-}
    (//)          = replaceMatrix
    {-# INLINE (//) #-}
    unsafeReplace = unsafeReplaceMatrix
    {-# INLINE unsafeReplace #-}
    unsafeAt      = unsafeAtMatrix
    {-# INLINE unsafeAt #-}
    indices       = indicesMatrix
    {-# INLINE indices #-}
    elems         = elemsMatrix
    {-# INLINE elems #-}
    assocs        = assocsMatrix
    {-# INLINE assocs #-}
    tmap          = tmapMatrix
    {-# INLINE tmap #-}
    (*>) k (Matrix a) = unsafePerformIO $
        unsafeFreezeIOMatrix =<< getScaledIOMatrix k a
    {-# INLINE (*>) #-}
    shift k (Matrix a) = unsafePerformIO $
        unsafeFreezeIOMatrix =<< getShiftedIOMatrix k a
    {-# INLINE shift #-}

instance (Monad m) => ReadTensor Matrix (Int,Int) m where
    getSize = return . size
    {-# INLINE getSize #-}
    getAssocs = return . assocs
    {-# INLINE getAssocs #-}
    getIndices = return . indices
    {-# INLINE getIndices #-}
    getElems = return . elems
    {-# INLINE getElems #-}
    getAssocs' = return . assocs
    {-# INLINE getAssocs' #-}
    getIndices' = return . indices
    {-# INLINE getIndices' #-}
    getElems' = return . elems
    {-# INLINE getElems' #-}
    unsafeReadElem x i = return $ unsafeAt x i
    {-# INLINE unsafeReadElem #-}

instance IMatrix Matrix where
    unsafeSApplyVector alpha a x = unsafePerformIO $
        unsafeFreezeIOVector =<< unsafeGetSApplyVector alpha a x
    {-# INLINE unsafeSApplyVector #-}    
    unsafeSApplyMatrix alpha a b = unsafePerformIO $
        unsafeFreezeIOMatrix =<< unsafeGetSApplyMatrix alpha a b
    {-# INLINE unsafeSApplyMatrix #-}   
    unsafeRow = unsafeRowView
    {-# INLINE unsafeRow #-}   
    unsafeCol = unsafeColView
    {-# INLINE unsafeCol #-}

instance IMatrix (Herm Matrix) where
    unsafeSApplyVector alpha a x = unsafePerformIO $ 
        unsafeFreezeIOVector =<< unsafeGetSApplyVector alpha a x
    {-# INLINE unsafeSApplyVector #-}    
    unsafeSApplyMatrix alpha a b = unsafePerformIO $
        unsafeFreezeIOMatrix =<< unsafeGetSApplyMatrix alpha a b    
    {-# INLINE unsafeSApplyMatrix #-}   
    unsafeRow a i = unsafePerformIO $
        unsafeFreezeIOVector =<< unsafeGetRow a i
    {-# INLINE unsafeRow #-}
    unsafeCol a j = unsafePerformIO $
        unsafeFreezeIOVector =<< unsafeGetCol a j
    {-# INLINE unsafeCol #-}

instance IMatrix (Tri Matrix) where
    unsafeSApplyVector alpha a x = unsafePerformIO $ 
        unsafeFreezeIOVector =<< unsafeGetSApplyVector alpha a x
    {-# INLINE unsafeSApplyVector #-}
    unsafeSApplyMatrix alpha a b = unsafePerformIO $
        unsafeFreezeIOMatrix =<< unsafeGetSApplyMatrix alpha a b
    {-# INLINE unsafeSApplyMatrix #-}   
    unsafeRow a i = unsafePerformIO $
        unsafeFreezeIOVector =<< unsafeGetRow a i
    {-# INLINE unsafeRow #-}
    unsafeCol a j = unsafePerformIO $
        unsafeFreezeIOVector =<< unsafeGetCol a j
    {-# INLINE unsafeCol #-}

instance (MonadInterleave m) => MMatrix Matrix m where
    unsafeDoSApplyAddVector = gemv
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = gemm
    {-# INLINE unsafeDoSApplyAddMatrix #-}
    unsafeGetRow = unsafeGetRowMatrix
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColMatrix
    {-# INLINE unsafeGetCol #-}

instance (MonadInterleave m) => MMatrix (Herm Matrix) m where
    unsafeDoSApplyAddVector = hemv
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = hemm
    {-# INLINE unsafeDoSApplyAddMatrix #-}    
    unsafeGetCol = unsafeGetColHermMatrix
    {-# INLINE unsafeGetCol #-}

instance (MonadInterleave m) => MMatrix (Tri Matrix) m where
    unsafeDoSApplyAddVector = unsafeDoSApplyAddVectorTriMatrix
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = unsafeDoSApplyAddMatrixTriMatrix
    {-# INLINE unsafeDoSApplyAddMatrix #-}
    unsafeDoSApplyVector_  = trmv
    {-# INLINE unsafeDoSApplyVector_  #-}
    unsafeDoSApplyMatrix_ = trmm
    {-# INLINE unsafeDoSApplyMatrix_ #-}
    unsafeGetCol = unsafeGetColTriMatrix
    {-# INLINE unsafeGetCol #-}
    unsafeGetRow = unsafeGetRowTriMatrix
    {-# INLINE unsafeGetRow #-}

instance (MonadInterleave m) => MSolve (Tri Matrix) m where
    unsafeDoSSolveVector = unsafeDoSSolveTriMatrix
    {-# INLINE unsafeDoSSolveVector #-}
    unsafeDoSSolveMatrix = unsafeDoSSolveMatrixTriMatrix
    {-# INLINE unsafeDoSSolveMatrix #-}    
    unsafeDoSSolveVector_ = trsv
    {-# INLINE unsafeDoSSolveVector_ #-}
    unsafeDoSSolveMatrix_ = trsm
    {-# INLINE unsafeDoSSolveMatrix_ #-}

compareMatrixWith :: (e -> e -> Bool) 
                  -> Matrix e -> Matrix e -> Bool
compareMatrixWith cmp a b
    | shape a /= shape b =
        False
    | isHermMatrix a == isHermMatrix b =
        let elems' = if isHermMatrix a then elems . herm
                                       else elems
        in
            and $ zipWith cmp (elems' a) (elems' b)
    | otherwise =
        and $ zipWith cmp (colElems a) (colElems b)
  where
    colElems c = concatMap elems (colViews c)

instance (Eq e) => Eq (Matrix e) where
    (==) = compareMatrixWith (==)
    {-# INLINE (==) #-}

instance (AEq e) => AEq (Matrix e) where
    (===) = compareMatrixWith (===)
    {-# INLINE (===) #-}
    (~==) = compareMatrixWith (~==)
    {-# INLINE (~==) #-}
instance (Show e) => Show (Matrix e) where
    show a | isHermMatrix a = 
                "herm (" ++ show (herm a) ++ ")"
           | otherwise =
                "listMatrix " ++ show (shape a) ++ " " ++ show (elems a)
        
instance (BLAS1 e) => Num (Matrix e) where
    (+) x y     = unsafePerformIO $ unsafeFreezeIOMatrix =<< getAddMatrix x y
    {-# INLINE (+) #-}
    (-) x y     = unsafePerformIO $ unsafeFreezeIOMatrix =<< getSubMatrix x y
    {-# INLINE (-) #-}
    (*) x y     = unsafePerformIO $ unsafeFreezeIOMatrix =<< getMulMatrix x y
    {-# INLINE (*) #-}
    negate      = ((-1) *>)
    {-# INLINE negate #-}
    abs         = tmap abs
    {-# INLINE abs #-}
    signum      = tmap signum
    {-# INLINE signum #-}
    fromInteger = (constantMatrix (1,1)) . fromInteger
    
instance (BLAS3 e) => Fractional (Matrix e) where
    (/) x y      = unsafePerformIO $ unsafeFreezeIOMatrix =<< getDivMatrix x y
    {-# INLINE (/) #-}
    recip        = tmap recip
    {-# INLINE recip #-}
    fromRational = (constantMatrix (1,1)) . fromRational 

instance (BLAS3 e, Floating e) => Floating (Matrix e) where
    pi    = constantMatrix (1,1) pi
    exp   = tmap exp
    {-# INLINE exp #-}
    sqrt  = tmap sqrt
    {-# INLINE sqrt #-}
    log   = tmap log
    {-# INLINE log #-}
    (**)  = tzipWithMatrix (**)
    {-# INLINE (**) #-}
    sin   = tmap sin
    {-# INLINE sin #-}
    cos   = tmap cos
    {-# INLINE cos #-}
    tan   = tmap tan
    {-# INLINE tan #-}
    asin  = tmap asin
    {-# INLINE asin #-}
    acos  = tmap acos
    {-# INLINE acos #-}
    atan  = tmap atan
    {-# INLINE atan #-}
    sinh  = tmap sinh
    {-# INLINE sinh #-}
    cosh  = tmap cosh
    {-# INLINE cosh #-}
    tanh  = tmap tanh
    {-# INLINE tanh #-}
    asinh = tmap asinh
    {-# INLINE asinh #-}
    acosh = tmap acosh
    {-# INLINE acosh #-}
    atanh = tmap atanh
    {-# INLINE atanh #-}


-- | Take a unary elementwise vector operation and apply it to the elements
-- of a matrix.
liftMatrix :: (ReadMatrix a m) =>
    (VectorView a e -> m ()) -> a e -> m ()
liftMatrix f a =
    case maybeViewMatrixAsVector a of
        Just x -> f x
        _ ->
            let xs = case isHermMatrix a of
                          True ->  rowViews a
                          False -> colViews a
                go (y:ys) = do f y
                               go ys
                go []     = return ()
            in go xs
{-# INLINE liftMatrix #-}

-- | Take a binary elementwise vector operation and apply it to the elements
-- of a pair of matrices.
liftMatrix2 :: (ReadMatrix a m, ReadMatrix b m) =>
    (VectorView a e -> VectorView b f -> m ()) ->
        a e -> b f -> m ()
liftMatrix2 f a b =
    if isHermMatrix a == isHermMatrix b
        then case (maybeViewMatrixAsVector a, maybeViewMatrixAsVector b) of
                 ((Just x), (Just y)) -> f x y
                 _                    -> elementwise
        else elementwise
  where
    elementwise =
        let vecsA = if isHermMatrix a then rowViews else colViews
            vecsB = if isHermMatrix a then rowViews else colViews
            go (x:xs) (y:ys) = do f x y
                                  go xs ys
            go []     []     = return ()
            go _      _      = error $ printf 
                ("liftMatrix2 <matrix of shape %s> <matrix of shape %s>:"
                ++ " shape mismatch") (show $ shape a) (show $ shape b)
        in go (vecsA a) (vecsB b)
{-# INLINE liftMatrix2 #-}

checkMatrixOp2 :: (BaseMatrix x, BaseMatrix y) => 
    (x e -> y f -> a) ->
        x e -> y f -> a
checkMatrixOp2 f x y = 
    checkBinaryOp (shape x) (shape y) $ f x y
{-# INLINE checkMatrixOp2 #-}

getUnaryMatrixOp :: (ReadMatrix a m, WriteMatrix b m) =>
    (b e -> m ()) -> a e -> m (b e)
getUnaryMatrixOp f a = do
    b  <- newCopyMatrix a
    f b
    return b
{-# INLINE getUnaryMatrixOp #-}

unsafeGetBinaryMatrixOp :: 
    (WriteMatrix c m, ReadMatrix a m, ReadMatrix b m) =>
    (c e -> b f -> m ()) ->
        a e -> b f -> m (c e)
unsafeGetBinaryMatrixOp f a b = do
    c  <- newCopyMatrix a
    f c b
    return c

indexOfMatrix :: (BaseMatrix a) => a e -> (Int,Int) -> Int
indexOfMatrix a (i,j) = 
    let (i',j') = case isHermMatrix a of
                        True  -> (j,i)
                        False -> (i,j)
        l = ldaMatrix a
    in i' + j'*l
{-# INLINE indexOfMatrix #-}

-- | Dense matrix in the 'ST' monad.  The type arguments are as follows:
--
--     * @s@: the state variable argument for the 'ST' type
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
newtype STMatrix s e = STMatrix (IOMatrix e)

-- | A safe way to create and work with a mutable matrix before returning 
-- an immutable matrix for later perusal. This function avoids copying
-- the matrix before returning it - it uses unsafeFreezeMatrix internally,
-- but this wrapper is a safe interface to that function. 
runSTMatrix :: (forall s . ST s (STMatrix s e)) -> Matrix e
runSTMatrix mx = 
    runST $ mx >>= \(STMatrix x) -> return (Matrix x)

instance HasVectorView (STMatrix s) where
    type VectorView (STMatrix s) = STVector s

instance Shaped (STMatrix s) (Int,Int) where
    shape (STMatrix a) = shapeIOMatrix a
    {-# INLINE shape #-}
    bounds (STMatrix a) = boundsIOMatrix a
    {-# INLINE bounds #-}

instance ReadTensor (STMatrix s) (Int,Int) (ST s) where
    getSize (STMatrix a) = unsafeIOToST $ getSizeIOMatrix a
    {-# INLINE getSize #-}
    unsafeReadElem (STMatrix a) i = unsafeIOToST $ unsafeReadElemIOMatrix a i
    {-# INLINE unsafeReadElem #-}
    getIndices (STMatrix a) = unsafeIOToST $ getIndicesIOMatrix a
    {-# INLINE getIndices #-}
    getIndices' (STMatrix a) = unsafeIOToST $ getIndicesIOMatrix' a
    {-# INLINE getIndices' #-}
    getElems (STMatrix a) = unsafeIOToST $ getElemsIOMatrix a
    {-# INLINE getElems #-}
    getElems' (STMatrix a) = unsafeIOToST $ getElemsIOMatrix' a
    {-# INLINE getElems' #-}
    getAssocs (STMatrix a) = unsafeIOToST $ getAssocsIOMatrix a
    {-# INLINE getAssocs #-}
    getAssocs' (STMatrix a) = unsafeIOToST $ getAssocsIOMatrix' a
    {-# INLINE getAssocs' #-}

instance WriteTensor (STMatrix s) (Int,Int) (ST s) where
    getMaxSize (STMatrix a) =  unsafeIOToST $ getMaxSizeIOMatrix a
    {-# INLINE getMaxSize #-}
    setZero (STMatrix a) = unsafeIOToST $ setZeroIOMatrix a
    {-# INLINE setZero #-}
    setConstant e (STMatrix a) = unsafeIOToST $ setConstantIOMatrix e a
    {-# INLINE setConstant #-}
    canModifyElem (STMatrix a) i = unsafeIOToST $ canModifyElemIOMatrix a i
    {-# INLINE canModifyElem #-}
    unsafeWriteElem (STMatrix a) i e = unsafeIOToST $ unsafeWriteElemIOMatrix a i e
    {-# INLINE unsafeWriteElem #-}
    unsafeModifyElem (STMatrix a) i f = unsafeIOToST $ unsafeModifyElemIOMatrix a i f
    {-# INLINE unsafeModifyElem #-}
    modifyWith f (STMatrix a) = unsafeIOToST $ modifyWithIOMatrix f a
    {-# INLINE modifyWith #-}
    doConj (STMatrix a) = unsafeIOToST $ doConjIOMatrix a
    {-# INLINE doConj #-}
    scaleBy k (STMatrix a) = unsafeIOToST $ scaleByIOMatrix k a
    {-# INLINE scaleBy #-}
    shiftBy k (STMatrix a) = unsafeIOToST $ shiftByIOMatrix k a
    {-# INLINE shiftBy #-}

instance MatrixShaped (STMatrix s) where

instance HasHerm (STMatrix s) where
    herm (STMatrix a) = STMatrix (herm a)
    {-# INLINE herm #-}
    
instance MMatrix (STMatrix s) (ST s) where
    unsafeDoSApplyAddVector = gemv
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = gemm
    {-# INLINE unsafeDoSApplyAddMatrix #-}
    unsafeGetRow = unsafeGetRowMatrix
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColMatrix
    {-# INLINE unsafeGetCol #-}

instance MMatrix (Herm (STMatrix s)) (ST s) where
    unsafeDoSApplyAddVector = hemv
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = hemm
    {-# INLINE unsafeDoSApplyAddMatrix #-}    
    unsafeGetCol = unsafeGetColHermMatrix
    {-# INLINE unsafeGetCol #-}

instance MMatrix (Tri (STMatrix s)) (ST s) where
    unsafeDoSApplyAddVector = unsafeDoSApplyAddVectorTriMatrix
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = unsafeDoSApplyAddMatrixTriMatrix
    {-# INLINE unsafeDoSApplyAddMatrix #-}
    unsafeDoSApplyVector_  = trmv
    {-# INLINE unsafeDoSApplyVector_  #-}
    unsafeDoSApplyMatrix_ = trmm
    {-# INLINE unsafeDoSApplyMatrix_ #-}
    unsafeGetRow = unsafeGetRowTriMatrix
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColTriMatrix
    {-# INLINE unsafeGetCol #-}

instance MSolve (Tri (STMatrix s)) (ST s) where
    unsafeDoSSolveVector = unsafeDoSSolveTriMatrix
    {-# INLINE unsafeDoSSolveVector #-}
    unsafeDoSSolveMatrix = unsafeDoSSolveMatrixTriMatrix
    {-# INLINE unsafeDoSSolveMatrix #-}    
    unsafeDoSSolveVector_ = trsv
    {-# INLINE unsafeDoSSolveVector_ #-}
    unsafeDoSSolveMatrix_ = trsm
    {-# INLINE unsafeDoSSolveMatrix_ #-}

instance BaseMatrix (STMatrix s) where
    ldaMatrix (STMatrix a) = ldaIOMatrix a
    {-# INLINE ldaMatrix #-}
    transEnumMatrix (STMatrix a) = transEnumIOMatrix a
    {-# INLINE transEnumMatrix #-}
    unsafeSubmatrixView (STMatrix a) ij mn =
        STMatrix (unsafeSubmatrixViewIOMatrix a ij mn)
    {-# INLINE unsafeSubmatrixView #-}
    unsafeDiagView (STMatrix a) i = STVector (unsafeDiagViewIOMatrix a i)
    {-# INLINE unsafeDiagView #-}
    unsafeRowView (STMatrix a) i = STVector (unsafeRowViewIOMatrix a i)
    {-# INLINE unsafeRowView #-}
    unsafeColView (STMatrix a) i = STVector (unsafeColViewIOMatrix a i)
    {-# INLINE unsafeColView #-}
    maybeViewMatrixAsVector (STMatrix a) = liftM STVector (maybeViewMatrixAsVector a)
    {-# INLINE maybeViewMatrixAsVector #-}
    maybeViewVectorAsMatrix mn (STVector x) = 
        liftM STMatrix $ maybeViewVectorAsIOMatrix mn x
    {-# INLINE maybeViewVectorAsMatrix #-}    
    maybeViewVectorAsRow (STVector x) = liftM STMatrix (maybeViewVectorAsRow x)
    {-# INLINE maybeViewVectorAsRow #-}
    maybeViewVectorAsCol (STVector x) = liftM STMatrix (maybeViewVectorAsCol x)
    {-# INLINE maybeViewVectorAsCol #-}
    unsafeIOMatrixToMatrix = STMatrix
    {-# INLINE unsafeIOMatrixToMatrix #-}
    unsafeMatrixToIOMatrix (STMatrix a) = a
    {-# INLINE unsafeMatrixToIOMatrix #-}

instance ReadMatrix (STMatrix s) (ST s) where
    freezeMatrix (STMatrix a) = unsafeIOToST $ freezeIOMatrix a
    {-# INLINE freezeMatrix #-}
    unsafeFreezeMatrix (STMatrix a) = unsafeIOToST $ unsafeFreezeIOMatrix a
    {-# INLINE unsafeFreezeMatrix #-}
    unsafePerformIOWithMatrix (STMatrix a) f = unsafeIOToST $ f a
    {-# INLINE unsafePerformIOWithMatrix #-}

instance WriteMatrix (STMatrix s) (ST s) where
    newMatrix_ = unsafeIOToST . liftM STMatrix . newIOMatrix_
    {-# INLINE newMatrix_ #-}
    thawMatrix = unsafeIOToST . liftM STMatrix . thawIOMatrix
    {-# INLINE thawMatrix #-}
    unsafeThawMatrix = unsafeIOToST . liftM STMatrix . unsafeThawIOMatrix
    {-# INLINE unsafeThawMatrix #-}
    unsafeConvertIOMatrix = unsafeIOToST . liftM STMatrix
    {-# INLINE unsafeConvertIOMatrix #-}
