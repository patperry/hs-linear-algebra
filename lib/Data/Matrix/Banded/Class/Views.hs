-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Class.Views
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Class.Views (
    -- * Row and column views
    diagViewBanded,
    rowViewBanded,
    colViewBanded,
    
    unsafeDiagViewBanded,
    unsafeRowViewBanded,
    unsafeColViewBanded,
    
    ) where

import BLAS.Internal( checkedRow, checkedCol, checkedDiag, diagStart, diagLen )
import Data.Matrix.Class( HasVectorView(..) )
import Data.Tensor.Class( shape )

import Data.Matrix.Banded.Class.Internal
import Data.Vector.Dense.Base
import Data.Vector.Dense.IOBase
import Foreign

diagViewBanded :: (BaseBanded a e) => 
    a mn e -> Int -> VectorView a k e
diagViewBanded a = checkedDiag (shape a) (unsafeDiagViewBanded a) 

rowViewBanded :: (BaseBanded a e) => 
    a mn e -> Int -> (Int, VectorView a k e, Int)
rowViewBanded a = checkedRow (shape a) (unsafeRowViewBanded a) 

colViewBanded :: (BaseBanded a e) => 
    a mn e -> Int -> (Int, VectorView a k e, Int)
colViewBanded a = checkedCol (shape a) (unsafeColViewBanded a)

unsafeDiagViewBanded :: (BaseBanded a e) => 
    a mn e -> Int -> VectorView a k e
unsafeDiagViewBanded a d
    | isHermBanded a = conj $ unsafeDiagViewBanded a' (negate d)
    | otherwise =
        let (fp,p,m,n,_,_,ld,_) = arrayFromBanded a
            off = indexOfBanded a (diagStart d)
            p'  = p `advancePtr` off
            len = diagLen (m,n) d
            inc = ld
            c   = False
        in unsafeIOVectorToVector (IOVector fp p' len inc c)
  where
    a' = (hermBanded . coerceBanded) a
