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
sum :: (BLAS1 e, VNum e) => Int -> [Vector e] -> Vector e
sum p xs = V.create $ do
    s <- V.new_ p
    sumTo xs s
    return s

-- | Returns the mean of the vectors.  The first argument gives the dimension
-- of the vectors.
mean :: (BLAS1 e, VNum e) => Int -> [Vector e] -> Vector e
mean p xs = V.create $ do
      m <- V.new_ p
      meanTo xs m
      return m

-- | Returns the weighted sum of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedSum :: (BLAS1 e, VNum e) => Int -> [(e, Vector e)] -> Vector e
weightedSum p wxs = V.create $ do
    s <- V.new_ p
    weightedSumTo wxs s
    return s

-- | Returns the weighted mean of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedMean :: (BLAS1 e, VNum e)
             => Int -> [(Double, Vector e)] -> Vector e
weightedMean p wxs = V.create $ do
       s <- V.new_ p
       weightedMeanTo wxs s
       return s

-- | Sets the target vector to the sum of the vectors.
sumTo :: (RVector v, BLAS1 e, VNum e) => [v e] -> STVector s e -> ST s ()
sumTo = weightedSumTo . zip (repeat 1)

-- | Sets the target vector to the mean of the vectors.
meanTo :: (RVector v, BLAS1 e, VNum e)
       => [v e] -> STVector s e -> ST s()
meanTo = weightedMeanTo . zip (repeat 1)

-- | Sets the target vector to the weigthed sum of the vectors.
weightedSumTo :: (RVector v, BLAS1 e, VNum e)
              => [(e, v e)] -> STVector s e -> ST s ()
weightedSumTo wxs s = do
    err <- V.new n 0
    old_s <- V.new_ n
    diff <- V.new_ n
    val <- V.new_ n
    
    V.setElems s (replicate n 0)
    forM_ wxs $ \(w,x) -> do
        V.unsafeCopyTo old_s s -- old_s := s
        
        V.unsafeCopyTo val x   -- val := w * x
        V.scaleM val w

        V.addTo err err val    -- err := err + val
        V.addTo s s err        -- s := s + err
        
        V.subTo diff old_s s   -- diff := old_s - s
        V.addTo err diff val   -- err := diff + val
  where
    n = V.dim s

-- | Sets the target vector to the weighted mean of the vectors.
weightedMeanTo :: (RVector v, BLAS1 e, VNum e)
               => [(Double, v e)] -> STVector s e -> ST s ()
weightedMeanTo wxs m = let
    go _ _ [] = return ()
    go diff w_sum ((w,x):wxs') | w == 0    = go diff w_sum wxs'
                               | otherwise = let w_sum' = w_sum + w
                                             in do
                                    V.subTo diff x m
                                    V.addWithScaleM m 
                                        (realToFrac $ w/w_sum') diff
                                    go diff w_sum' wxs'
    in do
        diff <- V.new_ n
        V.setElems m (replicate n 0)
        go diff 0 wxs
  where
    n = V.dim m
