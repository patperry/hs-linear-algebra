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
    addSumTo,
    meanTo,
    addWeightedSumTo,
    weightedMeanTo,

    ) where

import Prelude hiding ( sum )

import Control.Monad( forM_ )
import Control.Monad.ST( ST )

import Numeric.LinearAlgebra.Types

import Numeric.LinearAlgebra.Vector.Base( Vector )
import Numeric.LinearAlgebra.Vector.STBase( RVector, STVector )
import qualified Numeric.LinearAlgebra.Vector.STBase as V


-- | Returns the sum of the vectors.  The first argument gives the dimension
-- of the vectors.
sum :: (BLAS1 e) => Int -> [Vector e] -> Vector e
sum p xs = V.create $ do
    s <- V.new_ p
    V.clear s
    addSumTo s xs
    return s

-- | Returns the mean of the vectors.  The first argument gives the dimension
-- of the vectors.
mean :: (BLAS1 e) => Int -> [Vector e] -> Vector e
mean p xs = V.create $ do
      m <- V.new_ p
      meanTo m xs
      return m

-- | Returns the weighted sum of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedSum :: (BLAS1 e) => Int -> [(e, Vector e)] -> Vector e
weightedSum p wxs = V.create $ do
    s <- V.new_ p
    V.clear s
    addWeightedSumTo s wxs
    return s

-- | Returns the weighted mean of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedMean :: (BLAS1 e)
             => Int -> [(Double, Vector e)] -> Vector e
weightedMean p wxs = V.create $ do
    m <- V.new_ p
    weightedMeanTo m wxs
    return m

-- | Adds the sum of the vectors to the target vector.
addSumTo :: (RVector v, BLAS1 e) => STVector s e -> [v e] -> ST s ()
addSumTo dst = addWeightedSumTo dst . zip (repeat 1)

-- | Sets the target vector to the mean of the vectors.
meanTo :: (RVector v, BLAS1 e)
       => STVector s e -> [v e] -> ST s()
meanTo dst = weightedMeanTo dst . zip (repeat 1)

-- | Adds the weigthed sum of the vectors to the target vector.
addWeightedSumTo :: (RVector v, BLAS1 e)
                 => STVector s e ->  [(e, v e)] -> ST s ()
addWeightedSumTo s wxs = do
    n <- V.getDim s
    err <- V.new n 0
    old_s <- V.new_ n
    diff <- V.new_ n
    val <- V.new_ n
    
    forM_ wxs $ \(w,x) -> do
        V.unsafeCopyTo old_s s -- old_s := s
        
        V.unsafeCopyTo val x   -- val := w * x
        V.scaleByM_ w val

        V.addTo err err val    -- err := err + val
        V.addTo s s err        -- s := s + err
        
        V.subTo diff old_s s   -- diff := old_s - s
        V.addTo err diff val   -- err := diff + val

-- | Sets the target vector to the weighted mean of the vectors.
weightedMeanTo :: (RVector v, BLAS1 e)
               =>  STVector s e -> [(Double, v e)] -> ST s ()
weightedMeanTo m wxs = let
    go _ _ [] = return ()
    go diff w_sum ((w,x):wxs') | w == 0    = go diff w_sum wxs'
                               | otherwise = let w_sum' = w_sum + w
                                             in do
                                    V.subTo diff x m
                                    V.addWithScaleM_ 
                                        (realToFrac $ w/w_sum') diff m
                                    go diff w_sum' wxs'
    in do
        n <- V.getDim m
        diff <- V.new_ n
        V.clear m
        go diff 0 wxs
