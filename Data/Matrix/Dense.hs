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
    
    module BLAS.Matrix.Base,
    module BLAS.Matrix.Immutable,
    module BLAS.Tensor.Base,
    module BLAS.Tensor.Dense.Immutable,
    module BLAS.Tensor.Immutable,
    module BLAS.Tensor.Scalable,

    -- * Creating matrices
    matrix, 
    listMatrix,
    fromCols,
    fromRows,
    
    -- * Special matrices
    identity,

    -- * Rows and columns
    row,
    col,
    rows,
    cols,

    -- * Diagonals
    diag,

    -- * Augmenting matrices
    submatrix,

    -- * Matrix arithmetic
    shift,
    scale,
    invScale,

    -- * Casting matrices
    coerceMatrix,
    
    -- * Converting between vectors and matrices
    fromRow,
    fromCol,
    
    -- * Unsafe operations
    unsafeMatrix,
    unsafeRow,
    unsafeCol,
    unsafeDiag,
    unsafeSubmatrix,
    
    ) where

import Data.Maybe                  ( fromJust )
import System.IO.Unsafe            ( unsafePerformIO )

import BLAS.Access
import BLAS.Elem ( BLAS1, BLAS2 )
import BLAS.Matrix.Base hiding ( Matrix )
import BLAS.Matrix.Immutable
import BLAS.Tensor.Base
import BLAS.Tensor.Dense.Immutable
import BLAS.Tensor.Immutable
import BLAS.Tensor.Scalable

import Data.Matrix.Dense.Internal
import qualified Data.Matrix.Dense.Internal as M
import Data.Matrix.Dense.Operations ( shift, scale, invScale, plus, minus, 
    times, divide )
import Data.Vector.Dense hiding ( scale, invScale, shift )


-- | Create a new matrix of the given size and initialize the given elements to
-- the given values.  All other elements get set to zero.
matrix :: (BLAS1 e) => (Int,Int) -> [((Int,Int), e)] -> Matrix (m,n) e
matrix mn ies = unsafePerformIO $ newMatrix mn ies
{-# NOINLINE matrix #-}

-- | Create a new matrix with the given elements in row-major order.
listMatrix :: (Elem e) => (Int,Int) -> [e] -> Matrix (m,n) e
listMatrix mn es = unsafePerformIO $ newListMatrix mn es
{-# NOINLINE listMatrix #-}



-- | Same as 'matrix' but does not do any bounds checking.
unsafeMatrix :: (BLAS1 e) => (Int,Int) -> [((Int,Int), e)] -> Matrix (m,n) e
unsafeMatrix mn ies = unsafePerformIO $ unsafeNewMatrix mn ies
{-# NOINLINE unsafeMatrix #-}

-- | Create a matrix of the given shape from a list of columns
fromCols :: (BLAS1 e) => (Int,Int) -> [Vector m e] -> Matrix (m,n) e
fromCols mn cs = unsafePerformIO $ newColsMatrix mn cs
{-# NOINLINE fromCols #-}

-- | Create a matrix of the given shape from a list of rows
fromRows :: (BLAS1 e) => (Int,Int) -> [Vector n e] -> Matrix (m,n) e
fromRows mn rs = unsafePerformIO $ newRowsMatrix mn rs
{-# NOINLINE fromRows #-}

-- | Get a new matrix of the given shape with ones along the diagonal and
-- zeroes everywhere else.
identity :: (BLAS1 e) => (Int,Int) -> Matrix (m,n) e
identity mn = unsafePerformIO $ newIdentity mn
{-# NOINLINE identity #-}

-- | Get a matrix from a row vector.
fromRow :: (BLAS1 e) => Vector n e -> Matrix (one,n) e
fromRow x = 
    case maybeFromRow x of
        Just x' -> x'
        Nothing -> fromJust $ maybeFromRow $ unsafePerformIO $ newCopy x
{-# NOINLINE fromRow #-}

-- | Get a matrix from a column vector.
fromCol :: (BLAS1 e) => Vector m e -> Matrix (m,one) e
fromCol x = 
    case maybeFromCol x of
        Just x' -> x'
        Nothing -> fromJust $ maybeFromCol $ unsafePerformIO $ newCopy x
{-# NOINLINE fromCol #-}


instance (BLAS1 e) => Scalable (DMatrix Imm (m,n)) e where
    (*>) = scale
        
instance (BLAS2 e) => Num (DMatrix Imm (m,n) e) where
    (+)           = plus
    (-)           = minus
    (*)           = times
    negate        = scale (-1)
    abs           = amap abs
    signum        = amap signum
    fromInteger n = constant (1,1) (fromInteger n)
    
instance (BLAS2 e) => Fractional (DMatrix Imm (m,n) e) where
    (/) a b        = divide a b
    recip          = amap recip
    fromRational q = constant (1,1) (fromRational q)
    
instance (BLAS2 e, Floating e) => Floating (DMatrix Imm (m,n) e) where
    pi       = constant (1,1) pi
    exp      = amap exp 
    sqrt     = amap sqrt
    log      = amap log
    (**)     = azipWith (**)
    sin      = amap sin
    cos      = amap cos
    tan      = amap tan
    asin     = amap asin
    acos     = amap acos
    atan     = amap atan
    sinh     = amap sinh
    cosh     = amap cosh
    tanh     = amap tanh
    asinh    = amap asinh
    acosh    = amap acosh
    atanh    = amap atanh
    
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
    