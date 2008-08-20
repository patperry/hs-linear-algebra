{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense (
    -- * Dense matrix type
    Matrix,

    -- * Creating matrices
    matrix, 
    listMatrix,
    rowsMatrix,
    colsMatrix,
    rowMatrix,
    colMatrix,
    
    -- * Special matrices
    identityMatrix,

    -- * Converting between mutable and immutable matrices
    UnsafeFreezeMatrix(..),
    UnsafeThawMatrix(..),
    freezeMatrix,
    thawMatrix,

    module Data.Matrix.Dense.Class.Read,
    module BLAS.Tensor.Immutable,
    module BLAS.Matrix.RowCol.Immutable,
    module BLAS.Matrix.Diag.Immutable,
    module BLAS.Numeric.Immutable,
    
    -- * Unsafe operations
    unsafeMatrix,
    
    ) where

import Data.AEq
import System.IO.Unsafe

import BLAS.Elem ( Elem, BLAS1 )
import BLAS.Internal ( UnsafeIOToM(..), inlinePerformIO )

import BLAS.Tensor.Immutable
import BLAS.Matrix.RowCol.Immutable
import BLAS.Matrix.Diag.Immutable
import BLAS.Numeric.Immutable

import Data.Vector.Dense
import Data.Vector.Dense.IO( IOVector )

import Data.Matrix.Dense.IO hiding ( IOMatrix, liftMatrix, liftMatrix2 )
import qualified Data.Matrix.Dense.IO as IO
import Data.Matrix.Dense.Class.Write( newCopyMatrix )
import Data.Matrix.Dense.Class.Read hiding ( liftMatrix, liftMatrix2 )

newtype Matrix mn e = M (IO.IOMatrix mn e)

unsafeFreezeIOMatrix :: IO.IOMatrix mn e -> Matrix mn e
unsafeFreezeIOMatrix = M

unsafeThawIOMatrix :: Matrix mn e -> IO.IOMatrix mn e
unsafeThawIOMatrix (M a) = a

class UnsafeFreezeMatrix a where
    unsafeFreezeMatrix :: a mn e -> Matrix mn e
    
class UnsafeThawMatrix a where
    unsafeThawMatrix :: Matrix mn e -> a mn e
    
instance UnsafeFreezeMatrix IO.IOMatrix where
    unsafeFreezeMatrix = unsafeFreezeIOMatrix
instance UnsafeThawMatrix IO.IOMatrix where
    unsafeThawMatrix = unsafeThawIOMatrix
    
freezeMatrix :: (ReadMatrix a x e m, WriteMatrix b y e m, UnsafeFreezeMatrix b, BLAS1 e) =>
    a mn e -> m (Matrix mn e)
freezeMatrix x = do
    x' <- newCopyMatrix x
    return (unsafeFreezeMatrix x')

thawMatrix :: (WriteMatrix a y e m, BLAS1 e) =>
    Matrix mn e -> m (a mn e)
thawMatrix = undefined -- newCopyMatrix

liftMatrix :: (IO.IOMatrix n e -> a) -> Matrix n e -> a
liftMatrix f (M x) = f x
{-# INLINE liftMatrix #-}

liftMatrix2 :: 
    (IO.IOMatrix n e -> IO.IOMatrix n e -> a) -> 
        Matrix n e -> Matrix n e -> a
liftMatrix2 f x = liftMatrix (liftMatrix f x)
{-# INLINE liftMatrix2 #-}

unsafeLiftMatrix :: (IO.IOMatrix n e -> IO a) -> Matrix n e -> a
unsafeLiftMatrix f = unsafePerformIO . liftMatrix f
{-# NOINLINE unsafeLiftMatrix #-}

unsafeLiftMatrix2 :: 
    (IO.IOMatrix n e -> IO.IOMatrix n e -> IO a) -> 
        Matrix n e -> Matrix n e -> a
unsafeLiftMatrix2 f x y = unsafePerformIO $ liftMatrix2 f x y
{-# NOINLINE unsafeLiftMatrix2 #-}

inlineLiftMatrix :: (IO.IOMatrix n e -> IO a) -> Matrix n e -> a
inlineLiftMatrix f = inlinePerformIO . liftMatrix f
{-# INLINE inlineLiftMatrix #-}


-- | Create a new matrix of the given size and initialize the given elements to
-- the given values.  All other elements get set to zero.
matrix :: (BLAS1 e) => (Int,Int) -> [((Int,Int), e)] -> Matrix (m,n) e
matrix mn ies = unsafeFreezeIOMatrix $ unsafePerformIO $ newMatrix mn ies
{-# NOINLINE matrix #-}

-- | Create a new matrix with the given elements in row-major order.
listMatrix :: (Elem e) => (Int,Int) -> [e] -> Matrix (m,n) e
listMatrix mn es = unsafeFreezeIOMatrix $ unsafePerformIO $ newListMatrix mn es
{-# NOINLINE listMatrix #-}

-- | Same as 'matrix' but does not do any bounds checking.
unsafeMatrix :: (BLAS1 e) => (Int,Int) -> [((Int,Int), e)] -> Matrix (m,n) e
unsafeMatrix mn ies = unsafeFreezeIOMatrix $ unsafePerformIO $ unsafeNewMatrix mn ies
{-# NOINLINE unsafeMatrix #-}

-- | Create a matrix of the given shape from a list of rows
rowsMatrix :: (BLAS1 e) => (Int,Int) -> [Vector n e] -> Matrix (m,n) e
rowsMatrix mn rs = unsafeFreezeIOMatrix $ unsafePerformIO $ newRowsMatrix mn rs
{-# NOINLINE rowsMatrix #-}

-- | Create a matrix of the given shape from a list of columns
colsMatrix :: (BLAS1 e) => (Int,Int) -> [Vector m e] -> Matrix (m,n) e
colsMatrix mn cs = unsafeFreezeIOMatrix $ unsafePerformIO $ newColsMatrix mn cs
{-# NOINLINE colsMatrix #-}

-- | Get a new matrix of the given shape with ones along the diagonal and
-- zeroes everywhere else.
identityMatrix :: (BLAS1 e) => (Int,Int) -> Matrix (m,n) e
identityMatrix mn = unsafeFreezeIOMatrix $ unsafePerformIO $ newIdentityMatrix mn
{-# NOINLINE identityMatrix #-}

-- | Get a matrix from a row vector.
rowMatrix :: (BLAS1 e) => Vector n e -> Matrix (one,n) e
rowMatrix x = 
    case maybeFromRow $ unsafeThawIOVector x of
        Just x' -> unsafeFreezeIOMatrix x'
        Nothing -> unsafeFreezeIOMatrix $ unsafePerformIO $ newRowMatrix x
  where
    unsafeThawIOVector :: Vector n e -> IOVector n e
    unsafeThawIOVector = unsafeThawVector
{-# NOINLINE rowMatrix #-}

-- | Get a matrix from a column vector.
colMatrix :: (BLAS1 e) => Vector m e -> Matrix (m,one) e
colMatrix x = 
    case maybeFromCol $ unsafeThawIOVector x of
        Just x' -> unsafeFreezeIOMatrix x'
        Nothing -> unsafeFreezeIOMatrix $ unsafePerformIO $ newColMatrix x
  where
    unsafeThawIOVector :: Vector n e -> IOVector n e
    unsafeThawIOVector = unsafeThawVector
{-# NOINLINE colMatrix #-}

{-

instance (BLAS1 e, Show e) => Show (DMatrix Imm (m,n) e) where
    show a | isHerm a = 
                "herm (" ++ show (herm a) ++ ")"
           | otherwise =
                "listMatrix " ++ show (shape a) ++ " " ++ show (elems a)
        
compareHelp :: (BLAS1 e) => 
    (e -> e -> Bool) -> Matrix (m,n) e -> Matrix (m,n) e -> Bool
compareHelp cmp x y
    | isHerm x && isHerm y =
        compareHelp cmp (herm x) (herm y)
compareHelp cmp x y =
    (shape x == shape y) && (and $ zipWith cmp (elems x) (elems y))

instance (BLAS1 e, Eq e) => Eq (DMatrix Imm (m,n) e) where
    (==) = compareHelp (==)

instance (BLAS1 e, AEq e) => AEq (DMatrix Imm (m,n) e) where
    (===) = compareHelp (===)
    (~==) = compareHelp (~==)
-}
