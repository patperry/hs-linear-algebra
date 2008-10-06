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
    bandDiagView,
    bandRowView,
    bandColView,
    unsafeBandDiagView,
    unsafeBandRowView,
    unsafeBandColView,
    
    ) where

import BLAS.Internal( checkedRow, checkedCol, checkedDiag, diagStart, diagLen )

import Data.Matrix.Banded.Class.Internal
import Data.Vector.Dense.Class

bandDiagView :: (BaseBanded a x e) => a mn e -> Int -> x k e
bandDiagView a = checkedDiag (shape a) (unsafeBandDiagView a) 

bandRowView :: (BaseBanded a x e) => a mn e -> Int -> (Int, x k e, Int)
bandRowView a = checkedRow (shape a) (unsafeBandRowView a) 

bandColView :: (BaseBanded a x e) => a mn e -> Int -> (Int, x k e, Int)
bandColView a = checkedCol (shape a) (unsafeBandColView a)

unsafeBandDiagView :: (BaseBanded a x e) => a mn e -> Int -> x k e
unsafeBandDiagView a d
    | isHermBanded a = conj $ unsafeBandDiagView a' (negate d)
    | otherwise =
        let f   = fptrOfBanded a
            off = indexOfBanded a (diagStart d)
            len = diagLen (shape a) d
            inc = ldaOfBanded a
            c   = False
        in vectorViewArray f off len inc c
  where
    a' = (hermBanded . coerceBanded) a
