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
import Unsafe.Coerce

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
--     * @np@: a phantom type for the shape of the matrix.  Most functions
--       will demand that this be specified as a pair.  When writing a function
--       signature, you should always prefer @Matrix (n,p) e@ to
--       @Matrix np e@.
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
newtype Matrix np e = Matrix (IOMatrix np e)

freezeIOMatrix :: IOMatrix np e -> IO (Matrix np e)
freezeIOMatrix x@(IOMatrix _ _ _ _ _ _) = do
    y <- newCopyIOMatrix x
    return (Matrix y)
{-# INLINE freezeIOMatrix #-}

thawIOMatrix :: Matrix np e -> IO (IOMatrix np e)
thawIOMatrix (Matrix x@(IOMatrix _ _ _ _ _ _)) =
    newCopyIOMatrix x
{-# INLINE thawIOMatrix #-}

unsafeFreezeIOMatrix :: IOMatrix np e -> IO (Matrix np e)
unsafeFreezeIOMatrix = return . Matrix
{-# INLINE unsafeFreezeIOMatrix #-}

unsafeThawIOMatrix :: Matrix np e -> IO (IOMatrix np e)
unsafeThawIOMatrix (Matrix x) = return x
{-# INLINE unsafeThawIOMatrix #-}

-- | Common functionality for all dense matrix types.
class ( HasVectorView a, 
        MatrixShaped a, 
        BaseVector (VectorView a) ) => BaseMatrix a where
      
    -- | Get the leading dimension of the storage of the matrix.    
    ldaMatrix :: a (n,p) e -> Int
    
    -- | Get the storage type of the matrix.
    transEnumMatrix :: a (n,p) e -> TransEnum
    
    -- | Indicate whether or not the underlying matrix storage is
    -- transposed and conjugated.
    isHermMatrix :: a (n,p) e -> Bool
    isHermMatrix = (ConjTrans ==) . transEnumMatrix
    {-# INLINE isHermMatrix #-}

    -- | Cast the shape type of the matrix.
    coerceMatrix :: a np e -> a np' e
    coerceMatrix = unsafeCoerce
    {-# INLINE coerceMatrix #-}

    unsafeSubmatrixView :: a (n,p) e -> (Int,Int) -> (Int,Int) -> a (n',p') e
    unsafeDiagView :: a (n,p) e -> Int -> VectorView a k e
    unsafeRowView :: a (n,p) e -> Int -> VectorView a p e
    unsafeColView :: a (n,p) e -> Int -> VectorView a n e

    -- | Possibly create a vector view of a matrix.  This will fail if the
    -- matrix is hermed or if the lda of the matrix is not equal to the
    -- number of rows in the matrix.
    maybeViewMatrixAsVector :: a (n,p) e -> Maybe (VectorView a np e)
    
    -- | Possibly create a matrix view of a row vector.  This will fail if
    -- the stride of the vector is not @1@ or the vector is conjugated.
    maybeViewVectorAsRow  :: VectorView a p e -> Maybe (a (one,p) e)
    
    -- | Possibly create a matrix view of a column vector.  This will fail
    -- if the stride of the vector is not @1@ or the vector is not conjugated.
    maybeViewVectorAsCol  :: VectorView a n e -> Maybe (a (n,one) e)
    
    -- | Possible create a matrix view of the vector.  This will fail if
    -- the stride of the vector is not @1@ or the vector is conjugated.
    -- An error will be called if the vector does not have the same number
    -- of elements as the desired matrix.
    maybeViewVectorAsMatrix :: (Int,Int) 
                            -> VectorView a np e 
                            -> Maybe (a (n,p) e)
    
    -- | Unsafe cast from a matrix to an 'IOMatrix'.
    unsafeMatrixToIOMatrix :: a (n,p) e -> IOMatrix (n,p) e
    unsafeIOMatrixToMatrix :: IOMatrix (n,p) e -> a (n,p) e


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
    freezeMatrix :: a (n,p) e -> m (Matrix (n,p) e)
    unsafeFreezeMatrix :: a (n,p) e -> m (Matrix (n,p) e)

    -- | Cast the matrix to an 'IOMatrix', perform an @IO@ action, and
    -- convert the @IO@ action to an action in the monad @m@.
    unsafePerformIOWithMatrix :: a (n,p) e -> (IOMatrix (n,p) e -> IO r) -> m r


-- | Dense matrices that can be created or modified in a monad.
class ( ReadMatrix a m, 
        WriteTensor a (Int,Int) m, 
        WriteVector (VectorView a) m ) => WriteMatrix a m | m -> a where

    -- | Creates a new matrix of the given shape.  The elements will be 
    -- uninitialized.
    newMatrix_ :: (Elem e) => (Int,Int) -> m (a (n,p) e)

    -- | Convert an immutable matrix to a mutable one by taking a complete
    -- copy of it.
    thawMatrix :: Matrix (n,p) e -> m (a (n,p) e)
    unsafeThawMatrix :: Matrix (n,p) e -> m (a (n,p) e)

    -- | Unsafely convert an 'IO' action that creates an 'IOMatrix' into
    -- an action in @m@ that creates a matrix.
    unsafeConvertIOMatrix :: IO (IOMatrix (n,p) e) -> m (a (n,p) e)


-- | Creates a new matrix with the given association list.  Unspecified
-- indices will get initialized to zero.
newMatrix :: (WriteMatrix a m, Elem e)
          => (Int,Int) -> [((Int,Int), e)] -> m (a (n,p) e)
newMatrix mn ies = unsafeConvertIOMatrix $ 
    newIOMatrix mn ies
{-# INLINE newMatrix #-}
    
unsafeNewMatrix :: (WriteMatrix a m, Elem e)
                => (Int,Int) -> [((Int,Int), e)] -> m (a (n,p) e)
unsafeNewMatrix mn ies = unsafeConvertIOMatrix $
    unsafeNewIOMatrix mn ies
{-# INLINE unsafeNewMatrix #-}

-- | Create a new matrix with the given elements in column-major order.
newListMatrix :: (WriteMatrix a m, Elem e) 
              => (Int,Int) -> [e] -> m (a (n,p) e)
newListMatrix mn es = unsafeConvertIOMatrix $
    newListIOMatrix mn es
{-# INLINE newListMatrix #-}

-- | Form a matrix from a list of column vectors.
newColsMatrix :: (ReadVector x m, WriteMatrix a m, Elem e)
              => (Int,Int) -> [x n e] -> m (a (n,p) e)
newColsMatrix mn cs = unsafeConvertIOMatrix $
    newColsIOMatrix mn (map unsafeVectorToIOVector cs)
{-# INLINE newColsMatrix #-}

-- | Form a matrix from a list of row vectors.
newRowsMatrix :: (ReadVector x m, WriteMatrix a m, Elem e)
              => (Int,Int) -> [x p e] -> m (a (n,p) e)
newRowsMatrix mn rs = unsafeConvertIOMatrix $
    newRowsIOMatrix mn (map unsafeVectorToIOVector rs)
{-# INLINE newRowsMatrix #-}

-- | Create a new matrix from a column vector.
newColMatrix :: (ReadVector x m, WriteMatrix a m, Elem e) 
             => x n e -> m (a (n,one) e)
newColMatrix x = unsafeConvertIOMatrix $ 
    newColIOMatrix (unsafeVectorToIOVector x)
{-# INLINE newColMatrix #-}

-- | Create a new matrix from a row vector.
newRowMatrix :: (ReadVector x m, WriteMatrix a m, Elem e) 
             => x p e -> m (a (one,p) e)
newRowMatrix x = unsafeConvertIOMatrix $ 
    newRowIOMatrix (unsafeVectorToIOVector x)
{-# INLINE newRowMatrix #-}

-- | Create a zero matrix of the specified shape.
newZeroMatrix :: (WriteMatrix a m, Elem e) 
              => (Int,Int) -> m (a (n,p) e)
newZeroMatrix mn = unsafeConvertIOMatrix $ newZeroIOMatrix mn
{-# INLINE newZeroMatrix #-}

-- | Set every element in the matrix to zero.
setZeroMatrix :: (WriteMatrix a m) => a (n,p) e -> m ()
setZeroMatrix a = unsafePerformIOWithMatrix a $ setZeroIOMatrix
{-# INLINE setZeroMatrix #-}

-- | Create a constant matrix of the specified shape.
newConstantMatrix :: (WriteMatrix a m, Elem e) 
                  => (Int,Int) -> e -> m (a (n,p) e)
newConstantMatrix mn e = unsafeConvertIOMatrix $ newConstantIOMatrix mn e
{-# INLINE newConstantMatrix #-}

-- | Set every element in the matrix to the given constant.
setConstantMatrix :: (WriteMatrix a m) => e -> a (n,p) e -> m ()
setConstantMatrix e a = unsafePerformIOWithMatrix a $ setConstantIOMatrix e
{-# INLINE setConstantMatrix #-}

-- | Create a new matrix of the given shape with ones along the diagonal, 
-- and zeros everywhere else.
newIdentityMatrix :: (WriteMatrix a m, Elem e) => (Int,Int) -> m (a (n,p) e)
newIdentityMatrix = unsafeConvertIOMatrix . newIdentityIOMatrix
{-# INLINE newIdentityMatrix #-}

-- | Set diagonal elements to one and all other elements to zero.
setIdentityMatrix :: (WriteMatrix a m) => a (n,p) e -> m ()
setIdentityMatrix a =
    unsafePerformIOWithMatrix a $ setIdentityIOMatrix
{-# INLINE setIdentityMatrix #-}

-- | Get a copy of a matrix.
newCopyMatrix :: (ReadMatrix a m, WriteMatrix b m) =>
    a (n,p) e -> m (b (n,p) e)
newCopyMatrix a = unsafeConvertIOMatrix $ 
    newCopyIOMatrix (unsafeMatrixToIOMatrix a)
{-# INLINE newCopyMatrix #-}

-- | Get a copy of a matrix and make sure the returned matrix is not
-- a view.  Specififially, the returned matrix will have @isHermMatrix@
-- equal to @False@.
newCopyMatrix' :: (ReadMatrix a m, WriteMatrix b m) =>
    a (n,p) e -> m (b (n,p) e)
newCopyMatrix' a = unsafeConvertIOMatrix $
    newCopyIOMatrix' (unsafeMatrixToIOMatrix a)
{-# INLINE newCopyMatrix' #-}

-- | @copyMatrix dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyMatrix :: (WriteMatrix b m, ReadMatrix a m) => 
    b (n,p) e -> a (n,p) e -> m ()
copyMatrix b a = checkBinaryOp (shape b) (shape a) $ unsafeCopyMatrix b a
{-# INLINE copyMatrix #-}

unsafeCopyMatrix :: (WriteMatrix b m, ReadMatrix a m) =>
    b (n,p) e -> a (n,p) e -> m ()
unsafeCopyMatrix = liftMatrix2 unsafeCopyVector
{-# INLINE unsafeCopyMatrix #-}

-- | @swapMatrix x y@ swaps the values stored in two matrices.
swapMatrix :: (WriteMatrix a m) =>
    a (n,p) e -> a (n,p) e -> m ()
swapMatrix a b = checkBinaryOp (shape b) (shape a) $ unsafeSwapMatrix a b
{-# INLINE swapMatrix #-}

unsafeSwapMatrix :: (WriteMatrix a m) =>
    a (n,p) e -> a (n,p) e -> m ()
unsafeSwapMatrix = liftMatrix2 unsafeSwapVector
{-# INLINE unsafeSwapMatrix #-}

-- | Swap the elements in two rows of a matrix.
swapRows :: (WriteMatrix a m) => a (n,p) e -> Int -> Int -> m ()
swapRows a i j = 
    when (i /= j) $ unsafeSwapVector (rowView a i) (rowView a j)
{-# INLINE swapRows #-}

-- | Swap the elements in two columns of a matrix.
swapCols :: (WriteMatrix a m) => a (n,p) e -> Int -> Int -> m ()
swapCols a i j = 
    when (i /= j) $ unsafeSwapVector (colView a i) (colView a j)
{-# INLINE swapCols #-}

unsafeSwapRows :: (WriteMatrix a m) => a (n,p) e -> Int -> Int -> m ()
unsafeSwapRows a i j = 
    when (i /= j) $ unsafeSwapVector (unsafeRowView a i) (unsafeRowView a j)
{-# INLINE unsafeSwapRows #-}

unsafeSwapCols :: (WriteMatrix a m) => a (n,p) e -> Int -> Int -> m ()
unsafeSwapCols a i j = 
    when (i /= j) $ unsafeSwapVector (unsafeColView a i) (unsafeColView a j)
{-# INLINE unsafeSwapCols #-}

-- | @submatrixView a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrixView :: (BaseMatrix a) => a (n,p) e -> (Int,Int) -> (Int,Int) -> a (n',p') e
submatrixView a = checkedSubmatrix (shape a) (unsafeSubmatrixView a)
{-# INLINE submatrixView #-}

-- | Divide the rows of a matrix into two blocks and return views into the
-- blocks.  The integer argument indicates how many rows should be in the
-- first block.
splitRowsAt :: (BaseMatrix a) =>
    Int -> a (n,p) e -> (a (n1,p) e, a (n2,p) e)
splitRowsAt m1 a = ( submatrixView a (0,0)  (m1,n)
                   , submatrixView a (m1,0) (m2,n)
                   )
  where 
    (m,n) = shape a
    m2    = m - m1
{-# INLINE splitRowsAt #-}

unsafeSplitRowsAt :: (BaseMatrix a) =>
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
splitColsAt :: (BaseMatrix a) =>
    Int -> a (n,p) e -> (a (n,p1) e, a (n,p2) e)
splitColsAt n1 a = ( submatrixView a (0,0)  (m,n1)
                   , submatrixView a (0,n1) (m,n2)
                   )
  where
    (m,n) = shape a
    n2    = n - n1
{-# INLINE splitColsAt #-}

unsafeSplitColsAt :: (BaseMatrix a) =>
    Int -> a (n,p) e -> (a (n,p1) e, a (n,p2) e)
unsafeSplitColsAt n1 a = ( unsafeSubmatrixView a (0,0)  (m,n1)
                         , unsafeSubmatrixView a (0,n1) (m,n2)
                         )
  where
    (m,n) = shape a
    n2    = n - n1
{-# INLINE unsafeSplitColsAt #-}

-- | Get a list of vector views of the rows of the matrix.
rowViews :: (BaseMatrix a) => a (n,p) e -> [VectorView a p e]
rowViews a = [ unsafeRowView a i | i <- [0..numRows a - 1] ]
{-# INLINE rowViews #-}

-- | Get a list of vector views of the columns of the matrix.
colViews :: (BaseMatrix a) => a (n,p) e -> [VectorView a n e]
colViews a = [ unsafeColView a j | j <- [0..numCols a - 1] ]
{-# INLINE colViews #-}

-- | Get a vector view of the given row in a matrix.
rowView :: (BaseMatrix a) => a (n,p) e -> Int -> VectorView a p e
rowView a = checkedRow (shape a) (unsafeRowView a)
{-# INLINE rowView #-}

unsafeGetRowMatrix :: (ReadMatrix a m, WriteVector y m) =>
    a (n,p) e -> Int -> m (y p e)
unsafeGetRowMatrix a i = newCopyVector (unsafeRowView a i)
{-# INLINE unsafeGetRowMatrix #-}

-- | Get a vector view of the given column in a matrix.
colView :: (BaseMatrix a) => a (n,p) e -> Int -> VectorView a n e
colView a = checkedCol (shape a) (unsafeColView a)
{-# INLINE colView #-}

unsafeGetColMatrix :: (ReadMatrix a m, WriteVector y m) =>
    a (n,p) e -> Int -> m (y n e)
unsafeGetColMatrix a j = newCopyVector (unsafeColView a j)
{-# INLINE unsafeGetColMatrix #-}

-- | Get a vector view of the given diagonal in a matrix.
diagView :: (BaseMatrix a) => a (n,p) e -> Int -> VectorView a k e
diagView a = checkedDiag (shape a) (unsafeDiagView a)
{-# INLINE diagView #-}

-- | Get the given diagonal in a matrix.  Negative indices correspond
-- to sub-diagonals.
getDiag :: (ReadMatrix a m, WriteVector y m) => 
    a (n,p) e -> Int -> m (y k e)
getDiag a = checkedDiag (shape a) (unsafeGetDiag a)
{-# INLINE getDiag #-}

-- Same as 'getDiag' but not range-checked.
unsafeGetDiag :: (ReadMatrix a m, WriteVector y m) => 
    a (n,p) e -> Int -> m (y k e)
unsafeGetDiag a i = newCopyVector (unsafeDiagView a i)
{-# INLINE unsafeGetDiag #-}

-- | Conjugate every element of a matrix.
doConjMatrix :: (WriteMatrix a m, BLAS1 e) => a (n,p) e -> m ()
doConjMatrix = liftMatrix doConjVector
{-# INLINE doConjMatrix #-}

-- | Get a new matrix with elements with the conjugates of the elements
-- of the given matrix.
getConjMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    a (n,p) e -> m (b (n,p) e)
getConjMatrix = getUnaryMatrixOp doConjMatrix
{-# INLINE getConjMatrix #-}

-- | Scale every element of a matrix by the given value.
scaleByMatrix :: (WriteMatrix a m, BLAS1 e) => e -> a (n,p) e -> m ()
scaleByMatrix k = liftMatrix (scaleByVector k)
{-# INLINE scaleByMatrix #-}

-- | Get a new matrix by scaling the elements of another matrix
-- by a given value.
getScaledMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    e -> a (n,p) e -> m (b (n,p) e)
getScaledMatrix e = getUnaryMatrixOp (scaleByMatrix e)
{-# INLINE getScaledMatrix #-}

-- | Add a constant to every element in a matrix.
shiftByMatrix :: (WriteMatrix a m, BLAS1 e) => e -> a (n,p) e -> m ()
shiftByMatrix k = liftMatrix (shiftByVector k)
{-# INLINE shiftByMatrix #-}

-- | Get a new matrix by shifting the elements of another matrix
-- by a given value.
getShiftedMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    e -> a (n,p) e -> m (b (n,p) e)
getShiftedMatrix e = getUnaryMatrixOp (shiftByMatrix e)
{-# INLINE getShiftedMatrix #-}

-- | Replace the first argument with the elementwise sum.
addMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b (n,p) e -> a (n,p) e -> m ()
addMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeAddMatrix b a
{-# INLINE addMatrix #-}

unsafeAddMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b (n,p) e -> a (n,p) e -> m ()
unsafeAddMatrix b a = unsafeAxpyMatrix 1 a b
{-# INLINE unsafeAddMatrix #-}

-- | @getAddMatrix a b@ creates a new matrix equal to the sum @a+b@.  The 
-- operands must have the same shape.
getAddMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
getAddMatrix = checkMatrixOp2 unsafeGetAddMatrix
{-# INLINE getAddMatrix #-}

unsafeGetAddMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
unsafeGetAddMatrix = unsafeGetBinaryMatrixOp unsafeAddMatrix
{-# INLINE unsafeGetAddMatrix #-}

-- | Replace the first argument with the elementwise sum.
subMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b (n,p) e -> a (n,p) e -> m ()    
subMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeSubMatrix b a
{-# INLINE subMatrix #-}

unsafeSubMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b (n,p) e -> a (n,p) e -> m ()
unsafeSubMatrix b a = unsafeAxpyMatrix (-1) a b
{-# INLINE unsafeSubMatrix #-}

-- | @getSubMatrix a b@ creates a new matrix equal to the difference @a-b@.  The 
-- operands must have the same shape.
getSubMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
getSubMatrix = checkMatrixOp2 unsafeGetSubMatrix
{-# INLINE getSubMatrix #-}

unsafeGetSubMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
unsafeGetSubMatrix = unsafeGetBinaryMatrixOp unsafeSubMatrix
{-# INLINE unsafeGetSubMatrix #-}

-- | @axpyMatrix a x y@ replaces @y := a x + y@.
axpyMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    e -> a (n,p) e -> b (n,p) e -> m ()
axpyMatrix alpha x y = 
    checkBinaryOp (shape x) (shape y) $ unsafeAxpyMatrix alpha x y
{-# INLINE axpyMatrix #-}

unsafeAxpyMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS1 e) =>
    e -> a (n,p) e -> b (n,p) e -> m ()
unsafeAxpyMatrix alpha = liftMatrix2 (unsafeAxpyVector alpha)
{-# INLINE unsafeAxpyMatrix #-}

-- | Replace the first argument with the elementwise product.
mulMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b (n,p) e -> a (n,p) e -> m ()    
mulMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeMulMatrix b a
{-# INLINE mulMatrix #-}

unsafeMulMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b (n,p) e -> a (n,p) e -> m ()
unsafeMulMatrix = liftMatrix2 unsafeMulVector
{-# INLINE unsafeMulMatrix #-}

-- | @getMulMatrix a b@ creates a new matrix equal to the elementwise product 
-- @a*b@.  The operands must have the same shape.
getMulMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
getMulMatrix = checkMatrixOp2 unsafeGetMulMatrix
{-# INLINE getMulMatrix #-}

unsafeGetMulMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
unsafeGetMulMatrix = unsafeGetBinaryMatrixOp unsafeMulMatrix
{-# INLINE unsafeGetMulMatrix #-}

-- | Replace the first argument with the elementwise quotient.
divMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b (n,p) e -> a (n,p) e -> m ()    
divMatrix b a = 
    checkBinaryOp (shape b) (shape a) $ unsafeDivMatrix b a
{-# INLINE divMatrix #-}

unsafeDivMatrix :: (WriteMatrix b m, ReadMatrix a m, BLAS1 e) =>
    b (n,p) e -> a (n,p) e -> m ()
unsafeDivMatrix = liftMatrix2 unsafeDivVector
{-# INLINE unsafeDivMatrix #-}

-- | @getDivMatrix a b@ creates a new matrix equal to the elementwise ratio
-- @a/b@.  The operands must have the same shape.
getDivMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
getDivMatrix = checkMatrixOp2 unsafeGetDivMatrix
{-# INLINE getDivMatrix #-}

unsafeGetDivMatrix :: 
    (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    a (n,p) e -> b (n,p) e -> m (c (n,p) e)
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
class (MatrixShaped a, MonadInterleave m) => MMatrix a m where
    unsafeGetSApplyVector :: (ReadVector x m, WriteVector y m, BLAS3 e) =>
        e -> a (n,p) e -> x p e -> m (y n e)
    unsafeGetSApplyVector alpha a x = do
        y  <- newVector_ (numRows a)
        unsafeDoSApplyAddVector alpha a x 0 y
        return y
    {-# INLINE unsafeGetSApplyVector #-}

    unsafeGetSApplyMatrix :: (ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
        e -> a (n,p) e -> b (p,q) e -> m (c (n,q) e)
    unsafeGetSApplyMatrix alpha a b = do
        c  <- newMatrix_ (numRows a, numCols b)
        unsafeDoSApplyAddMatrix alpha a b 0 c
        return c
    {-# INLINE unsafeGetSApplyMatrix #-}

    unsafeDoSApplyAddVector :: (ReadVector x m, WriteVector y m, BLAS3 e) =>
        e -> a (n,p) e -> x p e -> e -> y n e -> m ()
    unsafeDoSApplyAddVector alpha a x beta (y :: y k e) = do
        (y' :: y k e) <- unsafeGetSApplyVector alpha a x
        scaleByVector beta y
        unsafeAxpyVector 1 y' y
    {-# INLINE unsafeDoSApplyAddVector #-}

    unsafeDoSApplyAddMatrix :: (ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
        e -> a (n,p) e -> b (p,q) e -> e -> c (n,q) e -> m ()
    unsafeDoSApplyAddMatrix alpha a b beta c = do
        c' <- unsafeGetSApplyMatrix alpha a b
        scaleByMatrix beta c
        unsafeAxpyMatrix 1 c' c
    {-# INLINE unsafeDoSApplyAddMatrix #-}

    unsafeDoSApplyVector_ :: (WriteVector y m, BLAS3 e) =>
        e -> a (n,n) e -> y n e -> m ()
    unsafeDoSApplyVector_ alpha a (x :: y n e) = do
        y <- newVector_ (dim x)
        unsafeDoSApplyAddVector alpha a x 0 y
        unsafeCopyVector x y
    {-# INLINE unsafeDoSApplyVector_  #-}

    unsafeDoSApplyMatrix_ :: (WriteMatrix b m, BLAS3 e) =>
        e -> a (n,n) e -> b (n,p) e -> m ()
    unsafeDoSApplyMatrix_ alpha a (b :: b (k,l) e) = do
        c <- newMatrix_ (shape b)
        unsafeDoSApplyAddMatrix alpha a b 0 c
        unsafeCopyMatrix b c
    {-# INLINE unsafeDoSApplyMatrix_ #-}

    unsafeGetCol :: (WriteVector x m, Elem e) => a (n,p) e -> Int -> m (x n e)
    
    unsafeGetRow :: (WriteVector x m, Elem e) => a (n,p) e -> Int -> m (x p e)
    unsafeGetRow a i = liftM conj $ unsafeGetCol (herm a) i
    {-# INLINE unsafeGetRow #-}

-- | Get a lazy list of the column vectors in the matrix.
getCols :: (MMatrix a m, WriteVector x m, Elem e)
        => a (n,p) e -> m [x n e]
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
        => a (n,p) e -> m [x p e]
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
    => a (n,p) e -> e -> x n e -> y p e -> m ()
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
    => a (n,p) e -> e -> x n e -> y p e -> m ()
unsafeRank1UpdateMatrix a alpha x y =
    unsafePerformIOWithMatrix a $ \a' ->
        unsafeRank1UpdateIOMatrix a' alpha (unsafeVectorToIOVector x)
                                           (unsafeVectorToIOVector y)
{-# INLINE unsafeRank1UpdateMatrix #-}

-- | @gemv alpha a x beta y@ replaces @y := alpha a * x + beta y@.
gemv :: (ReadMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
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
    e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
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
                       => Herm a (n,p) e -> Int -> m (x n e)
unsafeGetColHermMatrix (Herm l a) i =
    let a' = coerceMatrix a
        n  = numRows a'
        r  = unsafeRowView a' i
        c  = unsafeColView a' i 
        (r',c') = case l of { Lower -> (c,r) ; Upper -> (conj r, conj c) }
    in do
        x <- newVector_ n
        unsafeCopyVector (subvectorView x 0 i)     (conj $ subvectorView c' 0 i)
        unsafeCopyVector (subvectorView x i (n-i)) (subvectorView r' i (n-i))
        return x
{-# INLINE unsafeGetColHermMatrix #-}

hemv :: (ReadMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> Herm a (k,k) e -> x k e -> e -> y k e -> m ()
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
    e -> Herm a (k,k) e -> b (k,l) e -> e -> c (k,l) e -> m ()
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

hemv' :: (ReadMatrix a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> Herm a (r,s) e -> x s e -> e -> y r e -> m ()
hemv' alpha a x beta y = 
    hemv alpha (coerceHerm a) x beta (coerceVector y)
{-# INLINE hemv' #-}

hemm' :: (ReadMatrix a m, ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
    e -> Herm a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
hemm' alpha a b beta c = 
    hemm alpha (coerceHerm a) b beta (coerceMatrix c)
{-# INLINE hemm' #-}

unsafeDoSApplyAddVectorTriMatrix :: (ReadMatrix a m,
    ReadVector x m, WriteVector y m, BLAS2 e) =>
        e -> Tri a (k,l) e -> x l e -> e -> y k e -> m ()
unsafeDoSApplyAddVectorTriMatrix alpha t x beta (y :: y k e) =
    if beta == 0
        then unsafeDoSApplyTriMatrix alpha t x y
        else do
            (y' :: y k e) <- newCopyVector y
            unsafeDoSApplyTriMatrix alpha t x y'
            scaleByVector beta y
            unsafeAxpyVector 1 y' y
{-# INLINE unsafeDoSApplyAddVectorTriMatrix #-}

unsafeDoSApplyAddMatrixTriMatrix :: (ReadMatrix a m,
    ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
        e -> Tri a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
unsafeDoSApplyAddMatrixTriMatrix alpha t b beta (c :: c (r,t) e) =
    if beta == 0
        then unsafeDoSApplyMatrixTriMatrix alpha t b c
        else do
            (c' :: c (r,t) e) <- newCopyMatrix c
            unsafeDoSApplyMatrixTriMatrix alpha t b c'
            scaleByMatrix beta c
            unsafeAxpyMatrix 1 c' c
{-# INLINE unsafeDoSApplyAddMatrixTriMatrix #-}

unsafeDoSApplyTriMatrix :: (ReadMatrix a m,
    ReadVector x m, WriteVector y m, BLAS2 e) =>
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
            gemv alpha r x 0 y2
            
        (Upper,_,Left t') -> do
            unsafeCopyVector (coerceVector y) x
            trmv alpha t' (coerceVector y)

        (Upper,_,Right (t',r)) ->
            let x1 = unsafeSubvectorView x 0            (numCols t')
                x2 = unsafeSubvectorView x (numCols t') (numCols r)
            in do
                unsafeCopyVector y x1
                trmv alpha t' y
                gemv alpha r x2 1 y
  where
    (u,d,a) = triToBase t
{-# INLINE unsafeDoSApplyTriMatrix #-}

unsafeDoSApplyMatrixTriMatrix :: (ReadMatrix a m,
    ReadMatrix b m, WriteMatrix c m, BLAS3 e) =>
        e -> Tri a (r,s) e -> b (s,t) e -> c (r,t) e -> m ()
unsafeDoSApplyMatrixTriMatrix alpha t b c =
    case (u, toLower d a, toUpper d a) of
        (Lower,Left t',_) -> do
            unsafeCopyMatrix c (coerceMatrix b)
            trmm alpha t' c
            
        (Lower,Right (t',r),_) -> do
            let c1 = unsafeSubmatrixView c (0,0)          (numRows t',numCols c)
                c2 = unsafeSubmatrixView c (numRows t',0) (numRows r ,numCols c)
            unsafeCopyMatrix c1 b
            trmm alpha t' c1
            gemm alpha r b 0 c2
            
        (Upper,_,Left t') -> do
            unsafeCopyMatrix (coerceMatrix c) b
            trmm alpha t' (coerceMatrix c)

        (Upper,_,Right (t',r)) ->
            let b1 = unsafeSubmatrixView b (0,0)          (numCols t',numCols b)
                b2 = unsafeSubmatrixView b (numCols t',0) (numCols r ,numCols b)
            in do
                unsafeCopyMatrix c b1
                trmm alpha t' c
                gemm alpha r b2 1 c
  where
    (u,d,a) = triToBase t
{-# INLINE unsafeDoSApplyMatrixTriMatrix #-}

toLower :: (BaseMatrix a) => DiagEnum -> a (m,n) e 
        -> Either (Tri a (m,m) e) 
                  (Tri a (n,n) e, a (k,n) e)
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
    
toUpper :: (BaseMatrix a) => DiagEnum -> a (m,n) e
        -> Either (Tri a (n,n) e)
                  (Tri a (m,m) e, a (m,k) e)
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
                       => Tri a (n,p) e -> Int -> m (x n e)
unsafeGetColTriMatrix a@(Tri Upper _ _) j = 
    liftM conj $ unsafeGetRowTriMatrix (herm a) j

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
                      => Tri a (n,p) e -> Int -> m (x p e)
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
    e -> Tri a (k,k) e -> y n e -> m ()
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
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
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
                setZeroVector x2
  where
    (u,d,a) = triToBase t
{-# INLINE unsafeDoSSolveTriMatrix #-}

unsafeDoSSolveMatrixTriMatrix :: (ReadMatrix a m,
    ReadMatrix c m, WriteMatrix b m, BLAS3 e) =>
        e -> Tri a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
unsafeDoSSolveMatrixTriMatrix alpha t c b =
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
                setZeroMatrix b2
  where
    (u,d,a) = triToBase t
{-# INLINE unsafeDoSSolveMatrixTriMatrix #-}

trsv :: (ReadMatrix a m, WriteVector y m, BLAS2 e) =>
    e -> Tri a (k,k) e -> y n e -> m ()
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
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
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
    unsafeDoSolveVector :: (ReadVector y m, WriteVector x m, BLAS3 e) =>
        a (n,p) e -> y n e -> x p e -> m ()
    unsafeDoSolveVector = unsafeDoSSolveVector 1
    {-# INLINE unsafeDoSolveVector #-}
    
    unsafeDoSolveMatrix :: (ReadMatrix c m, WriteMatrix b m, BLAS3 e) =>
        a (n,p) e -> c (n,q) e -> b (p,q) e -> m ()
    unsafeDoSolveMatrix = unsafeDoSSolveMatrix 1
    {-# INLINE unsafeDoSolveMatrix #-}    
    
    unsafeDoSSolveVector :: (ReadVector y m, WriteVector x m, BLAS3 e) =>
        e -> a (n,p) e -> y n e -> x p e -> m ()
    unsafeDoSSolveVector alpha a y x = do
        unsafeDoSolveVector a y x
        scaleByVector alpha x
    {-# INLINE unsafeDoSSolveVector #-}        
    
    unsafeDoSSolveMatrix :: (ReadMatrix c m, WriteMatrix b m, BLAS3 e) =>
        e -> a (n,p) e -> c (n,q) e -> b (p,q) e -> m ()
    unsafeDoSSolveMatrix alpha a c b = do
        unsafeDoSolveMatrix a c b
        scaleByMatrix alpha b
    {-# INLINE unsafeDoSSolveMatrix #-}

    unsafeDoSolveVector_ :: (WriteVector x m, BLAS3 e) => a (n,n) e -> x n e -> m ()
    unsafeDoSolveVector_ = unsafeDoSSolveVector_ 1
    {-# INLINE unsafeDoSolveVector_ #-}

    unsafeDoSSolveVector_ :: (WriteVector x m, BLAS3 e) => e -> a (n,n) e -> x n e -> m ()
    unsafeDoSSolveVector_ alpha a x = do
        scaleByVector alpha x
        unsafeDoSolveVector_ a x
    {-# INLINE unsafeDoSSolveVector_ #-}        
        
    unsafeDoSolveMatrix_ :: (WriteMatrix b m, BLAS3 e) => a (n,n) e -> b (n,p) e -> m ()
    unsafeDoSolveMatrix_ = unsafeDoSSolveMatrix_ 1
    {-# INLINE unsafeDoSolveMatrix_ #-}
        
    unsafeDoSSolveMatrix_ :: (WriteMatrix b m, BLAS3 e) => e -> a (n,n) e -> b (n,p) e -> m ()
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
    unsafeDoSApplyAddVector = hemv'
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = hemm'
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
class (MatrixShaped a) => IMatrix a where
    unsafeSApplyVector :: (BLAS3 e) => e -> a (m,n) e -> Vector n e -> Vector m e
    unsafeSApplyMatrix :: (BLAS3 e) => e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e

    unsafeRow :: (Elem e) => a (m,n) e -> Int -> Vector n e
    unsafeRow a = conj . unsafeCol (herm a)
    {-# INLINE unsafeRow #-}

    unsafeCol :: (Elem e) => a (m,n) e -> Int -> Vector m e

-- | Get a list the row vectors in the matrix.
rows :: (IMatrix a, Elem e) => a (m,n) e -> [Vector n e]
rows a = [ unsafeRow a i | i <- [0..numRows a - 1] ]
{-# INLINE rows #-}

-- | Get a list the column vectors in the matrix.
cols :: (IMatrix a, Elem e) => a (m,n) e -> [Vector m e]
cols a = [ unsafeCol a j | j <- [0..numCols a - 1] ]
{-# INLINE cols #-}

-- | Create a new matrix of the given size and initialize the given elements to
-- the given values.  All other elements get set to zero.
matrix :: (Elem e) => (Int,Int) -> [((Int,Int), e)] -> Matrix (n,p) e
matrix mn ies = unsafePerformIO $
    unsafeFreezeIOMatrix =<< newIOMatrix mn ies
{-# INLINE matrix #-}

-- Same as 'matrix' but does not do any bounds checking.
unsafeMatrix :: (Elem e) => (Int,Int) -> [((Int,Int), e)] -> Matrix (n,p) e
unsafeMatrix mn ies =  unsafePerformIO $ 
    unsafeFreezeIOMatrix =<< unsafeNewIOMatrix mn ies
{-# INLINE unsafeMatrix #-}

-- | Create a new matrix with the given elements in row-major order.
listMatrix :: (Elem e) => (Int,Int) -> [e] -> Matrix (n,p) e
listMatrix mn es = unsafePerformIO $ 
    unsafeFreezeIOMatrix =<< newListIOMatrix mn es
{-# INLINE listMatrix #-}

replaceMatrix :: Matrix np e -> [((Int,Int),e)] -> Matrix np e
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

unsafeReplaceMatrix :: Matrix np e -> [((Int,Int),e)] -> Matrix np e
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
rowsMatrix :: (Elem e) => (Int,Int) -> [Vector p e] -> Matrix (n,p) e
rowsMatrix mn rs = unsafePerformIO $ 
    unsafeFreezeIOMatrix =<< newRowsIOMatrix mn rs
{-# INLINE rowsMatrix #-}

-- | Create a matrix of the given shape from a list of columns
colsMatrix :: (Elem e) => (Int,Int) -> [Vector n e] -> Matrix (n,p) e
colsMatrix mn cs = unsafePerformIO $ 
    unsafeFreezeIOMatrix =<< newColsIOMatrix mn cs
{-# INLINE colsMatrix #-}

-- | Get a matrix from a row vector.
matrixFromRow :: Vector p e -> Matrix (one,p) e
matrixFromRow (Vector x@(IOVector _ _ _ _ _)) = 
    case maybeViewVectorAsRow x of
        Just x' -> Matrix x'
        Nothing -> unsafePerformIO $ unsafeFreezeIOMatrix =<< newRowIOMatrix x
{-# INLINE matrixFromRow #-}

-- | Get a matrix from a column vector.
matrixFromCol :: Vector n e -> Matrix (n,one) e
matrixFromCol (Vector x@(IOVector _ _ _ _ _)) = 
    case maybeViewVectorAsCol x of
        Just x' -> Matrix x'
        Nothing -> unsafePerformIO $ unsafeFreezeIOMatrix =<< newColIOMatrix x
{-# INLINE matrixFromCol #-}

-- | Get a matrix from the elements stored in columnwise order in the vector.
matrixFromVector :: (Int,Int) -> Vector np e -> Matrix (n,p) e
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
vectorFromMatrix :: Matrix (n,p) e -> Vector np e
vectorFromMatrix a@(Matrix (IOMatrix _ _ _ _ _ _)) =
    case maybeViewMatrixAsVector a of
        Just x  -> x
        Nothing -> listVector (size a) (concatMap elems (colViews a))
{-# INLINE vectorFromMatrix #-}

-- | Get a new zero of the given shape.
zeroMatrix :: (Elem e) => (Int,Int) -> Matrix (n,p) e
zeroMatrix mn = unsafePerformIO $
    unsafeFreezeIOMatrix =<< newZeroIOMatrix mn
{-# INLINE zeroMatrix #-}

-- | Get a new constant of the given shape.
constantMatrix :: (Elem e) => (Int,Int) -> e -> Matrix (n,p) e
constantMatrix mn e = unsafePerformIO $
    unsafeFreezeIOMatrix =<< newConstantIOMatrix mn e
{-# INLINE constantMatrix #-}

-- | Get a new matrix of the given shape with ones along the diagonal and
-- zeroes everywhere else.
identityMatrix :: (Elem e) => (Int,Int) -> Matrix (n,p) e
identityMatrix mn = unsafePerformIO $
    unsafeFreezeIOMatrix =<< newIdentityIOMatrix mn
{-# INLINE identityMatrix #-}

-- | @submatrix a ij mn@ returns the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrix :: Matrix (n,p) e -> (Int,Int) -> (Int,Int) -> Matrix (n',p') e
submatrix (Matrix a) ij mn = 
    Matrix $ submatrixView a ij mn
{-# INLINE submatrix #-}

unsafeSubmatrix :: Matrix (n,p) e -> (Int,Int) -> (Int,Int) -> Matrix (n',p') e
unsafeSubmatrix (Matrix a) ij mn = 
    Matrix $ unsafeSubmatrixView a ij mn
{-# INLINE unsafeSubmatrix #-}

-- | Get a the given diagonal in a matrix.  Negative indices correspond to
-- sub-diagonals.
diag :: Matrix (n,p) e -> Int -> Vector k e
diag (Matrix a) i = Vector (diagView a i)
{-# INLINE diag #-}

-- Same as 'diag' but index is not range-checked.
unsafeDiag :: Matrix (n,p) e -> Int -> Vector k e
unsafeDiag (Matrix a) i = Vector (diagView a i)
{-# INLINE unsafeDiag #-}

unsafeAtMatrix :: Matrix np e -> (Int,Int) -> e
unsafeAtMatrix (Matrix a) = inlinePerformIO . unsafeReadElemIOMatrix a
{-# INLINE unsafeAtMatrix #-}

indicesMatrix :: Matrix np e -> [(Int,Int)]
indicesMatrix (Matrix a) = indicesIOMatrix a
{-# INLINE indicesMatrix #-}

elemsMatrix :: Matrix np e -> [e]
elemsMatrix (Matrix a) = 
    case maybeViewIOMatrixAsVector a of
        (Just x) -> elemsVector (Vector x)
        Nothing  -> concatMap (elemsVector . Vector) (vecViews a)
  where
    vecViews = if isHermIOMatrix a
                   then rowViews . coerceMatrix
                   else colViews . coerceMatrix
{-# INLINE elemsMatrix #-}

assocsMatrix :: Matrix np e -> [((Int,Int),e)]
assocsMatrix a = zip (indicesMatrix a) (elemsMatrix a)
{-# INLINE assocsMatrix #-}

tmapMatrix :: (e -> e) -> Matrix np e -> Matrix np e
tmapMatrix f a@(Matrix ma@(IOMatrix _ _ _ _ _ _))
    | isHermIOMatrix ma = coerceMatrix $ herm $ 
                              listMatrix (n,m) $ map (conjugate . f) (elems a)
    | otherwise         = coerceMatrix $
                              listMatrix (m,n) $ map f (elems a)
  where
    (m,n) = shape a
{-# INLINE tmapMatrix #-}

tzipWithMatrix :: (e -> e -> e) -> Matrix np e -> Matrix np e -> Matrix np e
tzipWithMatrix f a@(Matrix (IOMatrix _ _ _ _ _ _)) b
    | shape b /= mn =
        error ("tzipWith: matrix shapes differ; first has shape `" ++
                show mn ++ "' and second has shape `" ++
                show (shape b) ++ "'")
    | otherwise =
        coerceMatrix $
            listMatrix mn $ zipWith f (colElems a) (colElems b)
  where
    mn = shape a
    colElems = (concatMap elems) . colViews . coerceMatrix
{-# INLINE tzipWithMatrix #-}

instance Shaped Matrix (Int,Int) where
    shape (Matrix a) = shapeIOMatrix a
    {-# INLINE shape #-}
    bounds (Matrix a) = boundsIOMatrix a
    {-# INLINE bounds #-}

instance MatrixShaped Matrix where
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
    unsafeDoSApplyAddVector = hemv'
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = hemm'
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
                  -> Matrix (n,p) e -> Matrix (n,p) e -> Bool
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
    colElems c = concatMap elems (colViews $ coerceMatrix c)

instance (Eq e) => Eq (Matrix (n,p) e) where
    (==) = compareMatrixWith (==)
    {-# INLINE (==) #-}

instance (AEq e) => AEq (Matrix (n,p) e) where
    (===) = compareMatrixWith (===)
    {-# INLINE (===) #-}
    (~==) = compareMatrixWith (~==)
    {-# INLINE (~==) #-}
instance (Show e) => Show (Matrix (n,p) e) where
    show a | isHermMatrix a = 
                "herm (" ++ show (herm a) ++ ")"
           | otherwise =
                "listMatrix " ++ show (shape a) ++ " " ++ show (elems a)
        
instance (BLAS1 e) => Num (Matrix (n,p) e) where
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
    fromInteger = coerceMatrix . (constantMatrix (1,1)) . fromInteger
    
instance (BLAS3 e) => Fractional (Matrix (n,p) e) where
    (/) x y      = unsafePerformIO $ unsafeFreezeIOMatrix =<< getDivMatrix x y
    {-# INLINE (/) #-}
    recip        = tmap recip
    {-# INLINE recip #-}
    fromRational = coerceMatrix . (constantMatrix (1,1)) . fromRational 

instance (BLAS3 e, Floating e) => Floating (Matrix (m,n) e) where
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
    (forall k. VectorView a k e -> m ()) -> a (n,p) e -> m ()
liftMatrix f a =
    case maybeViewMatrixAsVector a of
        Just x -> f x
        _ ->
            let xs = case isHermMatrix a of
                          True ->  rowViews (coerceMatrix a)
                          False -> colViews (coerceMatrix a)
                go (y:ys) = do f y
                               go ys
                go []     = return ()
            in go xs
{-# INLINE liftMatrix #-}

-- | Take a binary elementwise vector operation and apply it to the elements
-- of a pair of matrices.
liftMatrix2 :: (ReadMatrix a m, ReadMatrix b m) =>
    (forall k. VectorView a k e -> VectorView b k f -> m ()) ->
        a (n,p) e -> b (n,p) f -> m ()
liftMatrix2 f a b =
    if isHermMatrix a == isHermMatrix b
        then case (maybeViewMatrixAsVector a, maybeViewMatrixAsVector b) of
                 ((Just x), (Just y)) -> f x y
                 _                    -> elementwise
        else elementwise
  where
    elementwise =
        let vecsA = if isHermMatrix a then rowViews . coerceMatrix
                                      else colViews . coerceMatrix
            vecsB = if isHermMatrix a then rowViews . coerceMatrix
                                      else colViews . coerceMatrix
            go (x:xs) (y:ys) = do f x y
                                  go xs ys
            go []     []     = return ()
            go _      _      = error $ printf 
                ("liftMatrix2 <matrix of shape %s> <matrix of shape %s>:"
                ++ " shape mismatch") (show $ shape a) (show $ shape b)
        in go (vecsA a) (vecsB b)
{-# INLINE liftMatrix2 #-}

checkMatrixOp2 :: (BaseMatrix x, BaseMatrix y) => 
    (x n e -> y n f -> a) ->
        x n e -> y n f -> a
checkMatrixOp2 f x y = 
    checkBinaryOp (shape x) (shape y) $ f x y
{-# INLINE checkMatrixOp2 #-}

getUnaryMatrixOp :: (ReadMatrix a m, WriteMatrix b m) =>
    (b (n,p) e -> m ()) -> a (n,p) e -> m (b (n,p) e)
getUnaryMatrixOp f a = do
    b  <- newCopyMatrix a
    f b
    return b
{-# INLINE getUnaryMatrixOp #-}

unsafeGetBinaryMatrixOp :: 
    (WriteMatrix c m, ReadMatrix a m, ReadMatrix b m) =>
    (c (n,p) e -> b (n,p) f -> m ()) ->
        a (n,p) e -> b (n,p) f -> m (c (n,p) e)
unsafeGetBinaryMatrixOp f a b = do
    c  <- newCopyMatrix a
    f c b
    return c

indexOfMatrix :: (BaseMatrix a) => a (n,p) e -> (Int,Int) -> Int
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
--     * @np@: a phantom type for the shape of the matrix.  Most functions
--       will demand that this be specified as a pair.  When writing a function
--       signature, you should always prefer @STMatrix s (n,p) e@ to
--       @STMatrix s np e@.
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
newtype STMatrix s np e = STMatrix (IOMatrix np e)

-- | A safe way to create and work with a mutable matrix before returning 
-- an immutable matrix for later perusal. This function avoids copying
-- the matrix before returning it - it uses unsafeFreezeMatrix internally,
-- but this wrapper is a safe interface to that function. 
runSTMatrix :: (forall s . ST s (STMatrix s (n,p) e)) -> Matrix (n,p) e
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
    unsafeDoSApplyAddVector = hemv'
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = hemm'
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
