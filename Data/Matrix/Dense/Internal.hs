{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Internal (
    -- * Dense matrix type
    Matrix(..),

    -- * Matrix shape
    module BLAS.Tensor.Base,
    module BLAS.Matrix.Base,
    coerceMatrix,

    -- * Creating matrices
    matrix, 
    listMatrix,
    rowsMatrix,
    colsMatrix,
    rowMatrix,
    colMatrix,
    unsafeMatrix,
    
    -- * Reading matrix elements
    module BLAS.Tensor.Immutable,
    
    -- * Special matrices
    zeroMatrix,
    constantMatrix,
    identityMatrix,

    -- * Matrix views
    submatrix,
    unsafeSubmatrix,
    
    -- * Vector views
    diag,
    unsafeDiag,

    -- * Matrix operations
    module BLAS.Numeric.Immutable,

    -- * Low-level properties
    lda,
    isHerm,
    ) where

import Data.AEq
import System.IO.Unsafe

import BLAS.Elem ( Elem, BLAS1 )
import BLAS.Internal ( inlinePerformIO )
import BLAS.UnsafeInterleaveM
import BLAS.UnsafeIOToM

import BLAS.Numeric.Immutable
import BLAS.Numeric

import BLAS.Tensor.Base
import BLAS.Tensor.Immutable
import BLAS.Tensor

import BLAS.Matrix.Base hiding ( BaseMatrix )
import qualified BLAS.Matrix.Base as BLAS
import BLAS.Numeric.Immutable

import Data.Matrix.Dense.Class.Creating
import Data.Matrix.Dense.Class.Special
import Data.Matrix.Dense.Class.Views( submatrix, unsafeSubmatrix,
    diagView, unsafeDiagView )
import Data.Matrix.Dense.Class.Internal( coerceMatrix, isHerm, lda, colViews,
    BaseMatrix(..), IOMatrix, maybeFromRow, maybeFromCol, newCopyMatrix,
    ReadMatrix )
import Data.Vector.Dense.Class.Internal
import Data.Vector.Dense

newtype Matrix mn e = M (IOMatrix mn e)

unsafeFreezeIOMatrix :: IOMatrix mn e -> Matrix mn e
unsafeFreezeIOMatrix = M

unsafeThawIOMatrix :: Matrix mn e -> IOMatrix mn e
unsafeThawIOMatrix (M a) = a


liftMatrix :: (IOMatrix n e -> a) -> Matrix n e -> a
liftMatrix f (M x) = f x
{-# INLINE liftMatrix #-}

liftMatrix2 :: 
    (IOMatrix n e -> IOMatrix n e -> a) -> 
        Matrix n e -> Matrix n e -> a
liftMatrix2 f x = liftMatrix (liftMatrix f x)
{-# INLINE liftMatrix2 #-}

unsafeLiftMatrix :: (IOMatrix n e -> IO a) -> Matrix n e -> a
unsafeLiftMatrix f = unsafePerformIO . liftMatrix f
{-# NOINLINE unsafeLiftMatrix #-}

unsafeLiftMatrix2 :: 
    (IOMatrix n e -> IOMatrix n e -> IO a) -> 
        Matrix n e -> Matrix n e -> a
unsafeLiftMatrix2 f x y = unsafePerformIO $ liftMatrix2 f x y
{-# NOINLINE unsafeLiftMatrix2 #-}

inlineLiftMatrix :: (IOMatrix n e -> IO a) -> Matrix n e -> a
inlineLiftMatrix f = inlinePerformIO . liftMatrix f
{-# INLINE inlineLiftMatrix #-}


-- | Create a new matrix of the given size and initialize the given elements to
-- the given values.  All other elements get set to zero.
matrix :: (BLAS1 e) => (Int,Int) -> [((Int,Int), e)] -> Matrix (m,n) e
matrix mn ies = unsafeFreezeIOMatrix $ unsafePerformIO $ newMatrix mn ies
{-# NOINLINE matrix #-}

-- | Create a new matrix with the given elements in row-major order.
listMatrix :: (BLAS1 e) => (Int,Int) -> [e] -> Matrix (m,n) e
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

-- | Get a new zero of the given shape.
zeroMatrix :: (BLAS1 e) => (Int,Int) -> Matrix (m,n) e
zeroMatrix mn = unsafeFreezeIOMatrix $ unsafePerformIO $ newZeroMatrix mn
{-# NOINLINE zeroMatrix #-}

-- | Get a new constant of the given shape.
constantMatrix :: (BLAS1 e) => (Int,Int) -> e -> Matrix (m,n) e
constantMatrix mn e = unsafeFreezeIOMatrix $ unsafePerformIO $ newConstantMatrix mn e
{-# NOINLINE constantMatrix #-}

-- | Get a new matrix of the given shape with ones along the diagonal and
-- zeroes everywhere else.
identityMatrix :: (BLAS1 e) => (Int,Int) -> Matrix (m,n) e
identityMatrix mn = unsafeFreezeIOMatrix $ unsafePerformIO $ newIdentityMatrix mn
{-# NOINLINE identityMatrix #-}


-- | Get a the given diagonal in a matrix.  Negative indices correspond to
-- sub-diagonals.
diag :: (Elem e) => Matrix mn e -> Int -> Vector k e
diag = diagView

-- | Same as 'diag' but index is not range-checked.
unsafeDiag :: (Elem e) => Matrix mn e -> Int -> Vector k e
unsafeDiag = unsafeDiagView


instance (Elem e) => BaseTensor Matrix (Int,Int) e where
    shape  = liftMatrix shape
    bounds = liftMatrix bounds

instance (BLAS1 e) => ITensor Matrix (Int,Int) e where
    constant mn e = coerceMatrix $ constantMatrix mn e
    zero mn       = coerceMatrix $ zeroMatrix mn

    (//)          = replaceHelp writeElem
    unsafeReplace = replaceHelp unsafeWriteElem
    
    unsafeAt a i  = inlineLiftMatrix (flip unsafeReadElem i) a
    {-# INLINE unsafeAt #-}
    
    size          = inlineLiftMatrix getSize
    elems         = inlineLiftMatrix getElems
    indices       = inlineLiftMatrix getIndices
    assocs        = inlineLiftMatrix getAssocs

    tmap f a 
        | isHerm a  = coerceMatrix $ herm $ 
                          listMatrix (n,m) $ map (conj . f) (elems a)
        | otherwise = coerceMatrix $
                          listMatrix (m,n) $ map f (elems a)
      where
        (m,n) = shape a


replaceHelp :: (BLAS1 e) => 
    (IOMatrix mn e -> (Int,Int) -> e -> IO ()) ->
        Matrix mn e -> [((Int,Int), e)] -> Matrix mn e
replaceHelp set x ies =
    unsafePerformIO $ do
        y  <- newCopyMatrix (unsafeThawIOMatrix x)
        mapM_ (uncurry $ set y) ies
        return (unsafeFreezeIOMatrix y)
{-# NOINLINE replaceHelp #-}

instance (BLAS1 e, Monad m) => ReadTensor Matrix (Int,Int) e m where
    getSize        = return . size
    getAssocs      = return . assocs
    getIndices     = return . indices
    getElems       = return . elems
    getAssocs'     = getAssocs
    getIndices'    = getIndices
    getElems'      = getElems
    unsafeReadElem x i = return (unsafeAt x i)

instance (BLAS1 e) => INumeric Matrix (Int,Int) e where
    (*>) k x = unsafeFreezeIOMatrix $ unsafeLiftMatrix (getScaled k) x
    {-# NOINLINE (*>) #-}

    shift k x = unsafeFreezeIOMatrix $ unsafeLiftMatrix (getShifted k) x
    {-# NOINLINE shift #-}

instance (BLAS1 e, Monad m) => ReadNumeric Matrix (Int,Int) e m where

instance (Elem e) => BLAS.BaseMatrix Matrix e where
    herm (M a) = M (herm a)
    
instance (Elem e) => BaseMatrix Matrix Vector e where
    matrixViewArray f o mn l h  = M $ matrixViewArray f o mn l h
    arrayFromMatrix (M a )      = arrayFromMatrix a

instance (BLAS1 e, UnsafeIOToM m, UnsafeInterleaveM m) => 
    ReadMatrix Matrix Vector e m where

instance (BLAS1 e) => Num (Matrix mn e) where
    (+) x y     = unsafeFreezeIOMatrix $ unsafeLiftMatrix2 getAdd x y
    (-) x y     = unsafeFreezeIOMatrix $ unsafeLiftMatrix2 getSub x y
    (*) x y     = unsafeFreezeIOMatrix $ unsafeLiftMatrix2 getMul x y
    negate      = ((-1) *>)
    abs         = tmap abs
    signum      = tmap signum
    fromInteger = (constant (1,1)) . fromInteger
    
instance (BLAS1 e) => Fractional (Matrix mn e) where
    (/) x y      = unsafeFreezeIOMatrix $ unsafeLiftMatrix2 getDiv x y
    recip        = tmap recip
    fromRational = (constant (1,1)) . fromRational 

instance (BLAS1 e, Floating e) => Floating (Matrix (m,n) e) where
    pi    = constant (1,1) pi
    exp   = tmap exp
    sqrt  = tmap sqrt
    log   = tmap log
    (**)  = tzipWith (**)
    sin   = tmap sin
    cos   = tmap cos
    tan   = tmap tan
    asin  = tmap asin
    acos  = tmap acos
    atan  = tmap atan
    sinh  = tmap sinh
    cosh  = tmap cosh
    tanh  = tmap tanh
    asinh = tmap asinh
    acosh = tmap acosh
    atanh = tmap atanh

tzipWith :: (BLAS1 e) =>
    (e -> e -> e) -> Matrix mn e -> Matrix mn e -> Matrix mn e
tzipWith f a b
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

instance (BLAS1 e, Show e) => Show (Matrix mn e) where
    show a | isHerm a = 
                "herm (" ++ show (herm $ coerceMatrix a) ++ ")"
           | otherwise =
                "listMatrix " ++ show (shape a) ++ " " ++ show (elems a)
        
compareHelp :: (BLAS1 e) => 
    (e -> e -> Bool) -> Matrix mn e -> Matrix mn e -> Bool
compareHelp cmp a b
    | shape a /= shape b =
        False
    | isHerm a == isHerm b =
        let elems' = if isHerm a then elems . herm .coerceMatrix
                                 else elems
        in
            and $ zipWith cmp (elems' a) (elems' b)
    | otherwise =
        and $ zipWith cmp (colElems a) (colElems b)
  where
    colElems c = concatMap elems (colViews $ coerceMatrix c)

instance (BLAS1 e, Eq e) => Eq (Matrix mn e) where
    (==) = compareHelp (==)

instance (BLAS1 e, AEq e) => AEq (Matrix mn e) where
    (===) = compareHelp (===)
    (~==) = compareHelp (~==)
