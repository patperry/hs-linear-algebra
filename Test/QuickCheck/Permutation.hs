-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Permutation
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Permutation
    where

import Data.List ( sortBy )
import Data.Ord ( comparing )

import Test.QuickCheck

import Data.Permutation ( Permutation )
import qualified Data.Permutation as P

permutation :: Int -> Gen Permutation
permutation n = do
    xs <- vector n :: Gen [Int]
    let is = (snd . unzip) $ sortBy (comparing fst) $ zip xs [0..]
    return $ P.permutation n is
    
newtype TestPermutation = TestPermutation Permutation deriving Show

instance Arbitrary TestPermutation where
    arbitrary = do
        n <- arbitrary >>= return . abs
        p <- permutation n
        return $ TestPermutation p
        
    coarbitrary (TestPermutation p) = 
            coarbitrary $ P.toList p
