-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Vector.Statistics
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Basic multivariate statistics.
--

module Numeric.LinearAlgebra.Vector.Statistics (
    -- * Immutable interface
    sum,
    mean,
    weightedSum,
    weightedMean,

    -- * Mutable interface
    sumTo,
    meanTo,
    weightedSumTo,
    weightedMeanTo,

    ) where

import Prelude hiding ( sum )

import Control.Monad( forM_ )
import Control.Monad.ST( ST )

import Numeric.LinearAlgebra.Types

import Numeric.LinearAlgebra.Vector.Base( Vector )
import Numeric.LinearAlgebra.Vector.STBase( RVector, STVector )
import qualified Numeric.LinearAlgebra.Vector.Base as V
import qualified Numeric.LinearAlgebra.Vector.STBase as V


-- | Returns the sum of the vectors.  The first argument gives the dimension
-- of the vectors.
sum :: (VNum e) => Int -> [Vector e] -> Vector e
sum p xs = V.create $ do
    s <- V.new_ p
    sumTo xs s
    return s

-- | Returns the mean of the vectors.  The first argument gives the dimension
-- of the vectors.
mean :: (VNum e, Fractional e) => Int -> [Vector e] -> Vector e
mean p xs = V.create $ do
      m <- V.new_ p
      meanTo xs m
      return m

-- | Returns the weighted sum of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedSum :: (VNum e) => Int -> [(e, Vector e)] -> Vector e
weightedSum p wxs = V.create $ do
    s <- V.new_ p
    weightedSumTo wxs s
    return s

-- | Returns the weighted mean of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedMean :: (VNum e, Fractional e) => Int -> [(Double, Vector e)] -> Vector e
weightedMean p wxs = V.create $ do
       s <- V.new_ p
       weightedMeanTo wxs s
       return s

-- | Sets the target vector to the sum of the vectors.
sumTo :: (RVector v, VNum e) => [v e] -> STVector s e -> ST s ()
sumTo = weightedSumTo . zip (repeat 1)

-- | Sets the target vector to the mean of the vectors.
meanTo :: (RVector v, VNum e, Fractional e)
       => [v e] -> STVector s e -> ST s()
meanTo = weightedMeanTo . zip (repeat 1)

-- | Sets the target vector to the weigthed sum of the vectors.
weightedSumTo :: (RVector v, VNum e) => [(e, v e)] -> STVector s e -> ST s ()
weightedSumTo wxs s = do
    err <- V.new n 0
    old_s <- V.new_ n
    diff <- V.new_ n
    val <- V.new_ n
    
    V.setElems s (replicate n 0)
    forM_ wxs $ \(w,x) -> do
        V.unsafeCopyTo s old_s -- old_s := s
        V.scaleTo w x val      -- val := w * x
        V.addTo err val err    -- err := err + val
        V.addTo s err s        -- s := s + err
        
        V.subTo old_s s diff   -- diff := old_s - s
        V.addTo diff val err   -- err := diff + val
  where
    n = V.dim s

-- | Sets the target vector to the weighted mean of the vectors.
weightedMeanTo :: (RVector v, VNum e, Fractional e)
               => [(Double, v e)] -> STVector s e -> ST s ()
weightedMeanTo wxs m = let
    go _ _ [] = return ()
    go diff w_sum ((w,x):wxs') | w == 0    = go diff w_sum wxs'
                               | otherwise = let w_sum' = w_sum + w
                                             in do
                                    V.subTo x m diff
                                    V.addToWithScales
                                        (realToFrac $ w/w_sum') diff 1 m m
                                    go diff w_sum' wxs'
    in do
        diff <- V.new_ n
        V.setElems m (replicate n 0)
        go diff 0 wxs
  where
    n = V.dim m
