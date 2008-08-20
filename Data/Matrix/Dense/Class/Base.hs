-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Base (
    BaseMatrix(..),
    
    submatrix,
    coerceMatrix,
    
    -- * Converting to and from vectors
    maybeFromRow,
    maybeFromCol,
    maybeToVector,
    liftMatrix,
    liftMatrix2,
    ) where

import Control.Monad( zipWithM_ )
import Foreign ( ForeignPtr, Ptr )
import Unsafe.Coerce

import BLAS.Internal( checkedSubmatrix )
import BLAS.Tensor( shape )

import BLAS.Matrix hiding ( BaseMatrix )
import qualified BLAS.Matrix as BLAS

import Data.Vector.Dense.Class

class (BLAS.BaseMatrix a e) => BaseMatrix a e where
    lda :: a mn e -> Int
    isHerm :: a mn e -> Bool
    unsafeSubmatrix :: a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
    withMatrixPtr :: a mn e -> (Ptr e -> IO b) -> IO b

    matrixViewArray :: ForeignPtr e -> Int -> (Int,Int) -> Int -> Bool -> a mn e
    arrayFromMatrix :: a mn e -> (ForeignPtr e, Int, (Int,Int), Int, Bool)
    
-- | @submatrix a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrix :: (BaseMatrix a e) => a mn e -> (Int,Int) -> (Int,Int) -> a mn' e
submatrix a = checkedSubmatrix (shape a) (unsafeSubmatrix a)
{-# INLINE submatrix #-}

-- | Cast the shape type of the matrix.
coerceMatrix :: (BaseMatrix a e) => a mn e -> a mn' e
coerceMatrix = unsafeCoerce
{-# INLINE coerceMatrix #-}

-- | Create a matrix view of a row vector.  This will fail if the
-- vector is conjugated and the stride is not @1@.
maybeFromRow :: (BaseMatrix a e, BaseVector x e) => 
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
maybeFromCol :: (BaseMatrix a e, BaseVector x e) => 
    x n e -> Maybe (a (n,one) e)
maybeFromCol x
    | c = maybeFromRow (conj x) >>= return . herm
    | s == 1 =
        Just $ matrixViewArray f o (n,1) (max 1 n) False
    | otherwise =
        Nothing
  where
    (f,o,n,s,c) = arrayFromVector x

maybeToVector :: (BaseMatrix a e, RowColView a x e) => 
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
liftMatrix :: (Monad m, BaseMatrix a e, RowColView a x e) => 
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
liftMatrix2 :: (Monad m, BaseMatrix a e, BaseMatrix b e,
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

