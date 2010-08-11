-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Statistics
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Basic multivariate statistics.
--

module Numeric.LinearAlgebra.Statistics (
    -- * Sums and means
    -- ** Immutable interface
    sumVector,
    meanVector,
    weightedSumVector,
    weightedMeanVector,

    -- ** Mutable interface
    sumToVector,
    meanToVector,
    weightedSumToVector,
    weightedMeanToVector,
    ) where

import Control.Monad( forM_ )
import Control.Monad.ST( ST )

import Foreign.VMath( VNum )
import Numeric.LinearAlgebra.Vector
import Numeric.LinearAlgebra.Vector.ST


-- | Returns the sum of the vectors.  The list must be nonempty.
sumVector :: (VNum e) => [Vector e] -> Vector e
sumVector xs | null xs = error "sumVector of empty list"
             | otherwise = runVector $ do
                 s <- newVector_ (dimVector $ head xs)
                 sumToVector xs s
                 return s
{-# INLINE sumVector #-}

-- | Returns the mean of the vectors.  The list must be nonempty.
meanVector :: (VNum e, Fractional e) => [Vector e] -> Vector e
meanVector xs | null xs = error "meanVector of empty list"
              | otherwise = runVector $ do
                  m <- newVector_ (dimVector $ head xs)
                  meanToVector xs m
                  return m
{-# INLINE meanVector #-}

-- | Returns the weighted sum of the vectors.  The list must be nonempty.
weightedSumVector :: (VNum e) => [(e,Vector e)] -> Vector e
weightedSumVector wxs | null wxs = error "weightedSumVector of empty list"
                      | otherwise = runVector $ do
                          s <- newVector_ (dimVector $ snd $ head wxs)
                          weightedSumToVector wxs s
                          return s
{-# INLINE weightedSumVector #-}

-- | Returns the weighted mean of the vectors.  The list must be nonempty.
weightedMeanVector :: (VNum e, Fractional e) => [(e,Vector e)] -> Vector e
weightedMeanVector wxs | null wxs = error "weightedMeanVector of empty list"
                       | otherwise = runVector $ do
                           s <- newVector_ (dimVector $ snd $ head wxs)
                           weightedMeanToVector wxs s
                           return s
{-# INLINE weightedMeanVector #-}

-- | Sets the target vector to the sum of the vectors.
sumToVector :: (RVector v, VNum e) => [v e] -> STVector s e -> ST s ()
sumToVector = weightedSumToVector . zip (repeat 1)
{-# INLINE sumToVector #-}

-- | Sets the target vector to the mean of the vectors.
meanToVector :: (RVector v, VNum e, Fractional e)
             => [v e] -> STVector s e -> ST s()
meanToVector = weightedMeanToVector . zip (repeat 1)
{-# INLINE meanToVector #-}

-- | Sets the target vector to the weigthed sum of the vectors.
weightedSumToVector :: (RVector v, VNum e) => [(e, v e)] -> STVector s e -> ST s ()
weightedSumToVector wxs s = do
    err <- newVector n 0
    old_s <- newVector_ n
    diff <- newVector_ n
    
    setElemsVector s (replicate n 0)
    forM_ wxs $ \(w,x) -> do
        unsafeCopyToVector s old_s
        addToVectorWithScales 1 err w x err  -- err  := err + val
        addToVector s err s                  -- sum  := sum + err
        
        subToVector old_s s diff             
        addToVectorWithScales 1 diff w x err -- err  := (old_sum - sum) + val    
  where
    n = dimVector s
{-# INLINE weightedSumToVector #-}

-- | Sets the target vector to the weighted mean of the vectors.
weightedMeanToVector :: (RVector v, VNum e, Fractional e)
                     => [(e, v e)] -> STVector s e -> ST s ()
weightedMeanToVector wxs m = let
    go _ _ [] = return ()
    go diff w_sum ((w,x):wxs') | w == 0    = go diff w_sum wxs'
                               | otherwise = let w_sum' = w_sum + w
                                             in do
                                    subToVector x m diff
                                    addToVectorWithScales (w/w_sum') diff 1 m m
                                    go diff w_sum' wxs'
    in do
        diff <- newVector_ n
        setElemsVector m (replicate n 0)
        go diff 0 wxs
  where
    n = dimVector m
{-# INLINE weightedMeanToVector #-}
