-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Class.Creating
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Class.Creating (
    -- * Creating banded matrices
    newBanded,
    newListsBanded,
    unsafeNewBanded,
    ) where

import Control.Monad
import Foreign

import BLAS.Internal( clearArray )
import Data.Tensor.Class.MTensor( writeElem, unsafeWriteElem )
import BLAS.UnsafeIOToM

import Data.Vector.Dense.Class( dim )

import Data.Matrix.Banded.Class.Internal
import Data.Matrix.Banded.Class.Views


newBanded :: (WriteBanded a e m) => 
    (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> m (a mn e)
newBanded = newBandedHelp writeElem

unsafeNewBanded :: (WriteBanded a e m) => 
    (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> m (a mn e)
unsafeNewBanded = newBandedHelp unsafeWriteElem

newBandedHelp :: (WriteBanded a e m) => 
       (a mn e -> (Int,Int) -> e -> m ()) 
    -> (Int,Int) -> (Int,Int) -> [((Int,Int),e)] -> m (a mn e)
newBandedHelp set (m,n) (kl,ku) ijes = do
    x <- newBanded_ (m,n) (kl,ku)
    unsafeIOToM $ withBandedPtr x $ flip clearArray ((kl+1+ku)*n)
    mapM_ (uncurry $ set x) ijes
    return x

newListsBanded :: (WriteBanded a e m) => 
    (Int,Int) -> (Int,Int) -> [[e]] -> m (a mn e)
newListsBanded (m,n) (kl,ku) xs = do
    a <- newBanded_ (m,n) (kl,ku)
    zipWithM_ (writeDiagElems a) [(negate kl)..ku] xs
    return a
  where
    writeDiagElems :: (WriteBanded a e m) => a mn e -> Int -> [e] -> m ()
    writeDiagElems a i es =
        let d   = diagViewBanded a i
            nb  = max 0 (negate i)
            es' = drop nb es
        in zipWithM_ (unsafeWriteElem d) [0..(dim d - 1)] es'

