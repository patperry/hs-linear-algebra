-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Class.Write
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Class.Write (
    WriteMatrix(..),
    
    -- * Creating new matrices
    newMatrix,
    newListMatrix,
    newRowsMatrix,
    newColsMatrix,
    unsafeNewMatrix,
    
    -- * Special matrices
    newZeroMatrix,
    newConstantMatrix,
    newIdentityMatrix,
    setIdentity,
        
    module Data.Matrix.Dense.Class.Read
    ) where

import Foreign

import BLAS.Elem
import BLAS.Matrix
import BLAS.Tensor
import Data.Vector.Dense.Class
import Data.Matrix.Dense.Class.Read


class (WriteTensor a (Int,Int) e m, WriteVector x e m, ReadMatrix a x e m) => WriteMatrix a x e m | a -> m, m -> a where
    -- | Creates a new matrix of the given shape.  The elements will be 
    -- uninitialized.
    newMatrix_ :: (Int,Int) -> m (a mn e)

---------------------------  Creating Matrices --------------------------------

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
    a <- newMatrix_ (m,n)
    unsafeIOToM $ withMatrixPtr a $ flip pokeArray (take (m*n) es)
    return a

-- | Form a matrix from a list of column vectors.
newColsMatrix :: (ReadVector x e m, WriteMatrix a y e m, BLAS1 e) => 
    (Int,Int) -> [x k e] -> m (a (k,l) e)
newColsMatrix (m,n) _ = do
    a <- newZero (m,n)
    --forM_ (zip [0..(n-1)] cs) $ \(j,c) ->
    --    V.copyVector (unsafeCol (unsafeThaw a) j) c
    return a

-- | Form a matrix from a list of row vectors.
newRowsMatrix :: (ReadVector x e m, WriteMatrix a y e m, BLAS1 e) => 
    (Int,Int) -> [x l e] -> m (a (k,l) e)
newRowsMatrix (m,n) _ = do
    a <- newZero (m,n)
    --forM_ (zip [0..(m-1)] rs) $ \(i,r) ->
    --    copy (unsafeRow (unsafeThaw a) i) r
    return a


---------------------------  Special Matrices --------------------------------

-- | Create a new matrix of the given shape with ones along the diagonal, 
-- and zeros everywhere else.
newIdentityMatrix :: (WriteMatrix a x e m) => (Int,Int) -> m (a mn e)
newIdentityMatrix mn = do
    a <- newMatrix_ mn
    setIdentity a
    return a


setIdentity :: (WriteMatrix a x e m) => a mn e -> m ()
setIdentity a = do
    setZero a
    mapM_ (\i -> unsafeWriteElem a (i,i) 1) [0..(mn-1)]
  where
    mn = min (numRows a) (numCols a)

-- | Create a zero matrix of the specified shape.
newZeroMatrix :: (WriteMatrix a x e m) => (Int,Int) -> m (a mn e)
newZeroMatrix = newZero
    
-- | Create a constant matrix of the specified shape.
newConstantMatrix :: (WriteMatrix a x e m) => (Int,Int) -> e -> m (a mn e)
newConstantMatrix = newConstant
