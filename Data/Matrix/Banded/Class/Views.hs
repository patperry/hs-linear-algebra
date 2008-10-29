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

import Data.Matrix.Banded.Class.Internal
import Data.Vector.Dense.Class

diagViewBanded :: (BaseBanded a x e) => a mn e -> Int -> x k e
diagViewBanded a = checkedDiag (shape a) (unsafeDiagViewBanded a) 

rowViewBanded :: (BaseBanded a x e) => a mn e -> Int -> (Int, x k e, Int)
rowViewBanded a = checkedRow (shape a) (unsafeRowViewBanded a) 

colViewBanded :: (BaseBanded a x e) => a mn e -> Int -> (Int, x k e, Int)
colViewBanded a = checkedCol (shape a) (unsafeColViewBanded a)

unsafeDiagViewBanded :: (BaseBanded a x e) => a mn e -> Int -> x k e
unsafeDiagViewBanded a d
    | isHermBanded a = conj $ unsafeDiagViewBanded a' (negate d)
    | otherwise =
        let f   = fptrOfBanded a
            off = indexOfBanded a (diagStart d)
            len = diagLen (shape a) d
            inc = ldaOfBanded a
            c   = False
        in vectorViewArray f off len inc c
  where
    a' = (hermBanded . coerceBanded) a
