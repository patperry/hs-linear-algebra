-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.RowCol.Read
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.RowCol.Read (
    getRow,
    getCol,
    
    getRows,
    getCols,
    getRows',
    getCols',
    
    RowColRead(..),
    ) where

import BLAS.Internal( checkedRow, checkedCol )
import BLAS.UnsafeInterleaveM
import BLAS.Tensor( shape )

import BLAS.Matrix.Base
import Data.Vector.Dense.Class( WriteVector )

class (BaseMatrix a e, Monad m) => RowColRead a e m where
    unsafeGetRow :: (WriteVector x e m) => a (k,l) e -> Int -> m (x l e)
    unsafeGetCol :: (WriteVector x e m) => a (k,l) e -> Int -> m (x k e)

-- | Get the given row in a matrix.
getRow :: (RowColRead a e m, WriteVector x e m) => a (k,l) e -> Int -> m (x l e)
getRow a = checkedRow (shape a) (unsafeGetRow a)

-- | Get the given column in a matrix.
getCol :: (RowColRead a e m, WriteVector x e m) => a (k,l) e -> Int -> m (x k e)
getCol a = checkedCol (shape a) (unsafeGetCol a)

-- | Get a lazy list the row vectors in the matrix.  See also "getRows'".
getRows :: (RowColRead a e m, WriteVector x e m, UnsafeInterleaveM m) => 
    a (k,l) e -> m [x l e]
getRows = unsafeInterleaveM . getRows'

-- | Get a lazy list of the column vectors in the matrix.  See also "getCols'".
getCols :: (RowColRead a e m, WriteVector x e m, UnsafeInterleaveM m) => 
    a (k,l) e -> m [x k e]
getCols = unsafeInterleaveM . getCols'

-- | Get a strict list the row vectors in the matrix.  See also "getRows".
getRows' :: (RowColRead a e m, WriteVector x e m) => a (k,l) e -> m [x l e]
getRows' a = mapM (unsafeGetRow a) [0..numRows a - 1]

-- | Get a strict list of the column vectors in the matrix.  See also "getCols".
getCols' :: (RowColRead a e m, WriteVector x e m) => a (k,l) e -> m [x k e]
getCols' a = mapM (unsafeGetCol a) [0..numCols a - 1]
