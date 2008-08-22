-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Creating
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Creating (
    -- * Creating matrices
    newMatrix_,
    newMatrix,
    newListMatrix,
    newRowsMatrix,
    newColsMatrix,
    newRowMatrix,
    newColMatrix,
    unsafeNewMatrix,
    ) where

import Control.Monad( forM_ )
import Foreign( pokeArray )

import BLAS.Elem
import BLAS.UnsafeIOToM
import BLAS.Matrix

import Data.Vector.Dense.Class
import Data.Matrix.Dense.Class.Internal


-- | Creates a new matrix with the given association list.  Unspecified
-- indices will get initialized to zero.
newMatrix :: (WriteMatrix a x e m) => (Int,Int) -> [((Int,Int), e)] -> m (a mn e)
newMatrix = newMatrixHelp writeElem

-- | Same as 'newMatrix' but indices are not range-checked.
unsafeNewMatrix :: (WriteMatrix a x e m) => (Int,Int) -> [((Int,Int), e)] -> m (a mn e)
unsafeNewMatrix = newMatrixHelp unsafeWriteElem

newMatrixHelp :: (WriteMatrix a x e m) => 
    (a mn e -> (Int,Int) -> e -> m ()) -> (Int,Int) -> [((Int,Int),e)] -> m (a mn e)
newMatrixHelp set n ies = do
    a <- newZeroMatrix n
    mapM_ (uncurry $ set a) ies
    return a

-- | Create a new matrix with the given elements in column-major order.
newListMatrix :: (WriteMatrix a x e m) => (Int,Int) -> [e] -> m (a mn e)
newListMatrix (m,n) es = do
    a <- newZeroMatrix (m,n)
    unsafeIOToM $ withMatrixPtr a $ flip pokeArray (take (m*n) es)
    return a

-- | Form a matrix from a list of column vectors.
newColsMatrix :: (ReadVector x e m, WriteMatrix a y e m, BLAS1 e) => 
    (Int,Int) -> [x k e] -> m (a (k,l) e)
newColsMatrix (m,n) cs = do
    a <- newZero (m,n)
    forM_ (zip [0..(n-1)] cs) $ \(j,c) ->
        copyVector (unsafeColView a j) c
    return a

-- | Form a matrix from a list of row vectors.
newRowsMatrix :: (ReadVector x e m, WriteMatrix a y e m, BLAS1 e) => 
    (Int,Int) -> [x l e] -> m (a (k,l) e)
newRowsMatrix (m,n) rs = do
    a <- newZero (m,n)
    forM_ (zip [0..(m-1)] rs) $ \(i,r) ->
        copyVector (unsafeRowView a i) r
    return a

-- | Create a new matrix from a column vector.
newColMatrix :: (ReadVector x e m, WriteMatrix a y e m, BLAS1 e) => 
    x k e -> m (a (k,one) e)
newColMatrix x = newColsMatrix (dim x,1) [x]

-- | Create a new matrix from a row vector.
newRowMatrix :: (ReadVector x e m, WriteMatrix a y e m, BLAS1 e) => 
    x l e -> m (a (one,l) e)
newRowMatrix x = newRowsMatrix (1,dim x) [x]

