{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.Vector.Dense
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.Vector.Dense (
    SubVector(..),
    VectorPair(..),
    VectorTriple(..),
    
    vector,
    ) where

import Test.QuickCheck hiding ( vector )
import Test.QuickCheck.BLASBase( SubVector(..) )
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense hiding ( vector )
import Data.Elem.BLAS ( BLAS1 )

vector :: (BLAS1 e, Arbitrary e) => Int -> Gen (Vector n e)
vector = Test.vector
    
data VectorPair n e = 
    VectorPair (Vector n e) 
               (Vector n e) 
    deriving (Show)
    
data VectorTriple n e = 
    VectorTriple (Vector n e) 
                 (Vector n e) 
                 (Vector n e) 
    deriving (Show)

instance (Arbitrary e, BLAS1 e) => Arbitrary (Vector n e) where
    arbitrary = sized $ \m ->
        choose (0,m) >>= vector
        
    coarbitrary = undefined

instance (Arbitrary e, BLAS1 e) => Arbitrary (VectorPair n e) where
    arbitrary = sized $ \m -> do
        n <- choose (0,m)
        x <- vector n
        y <- vector n
        return $ VectorPair x y
        
    coarbitrary = undefined
        
instance (Arbitrary e, BLAS1 e) => Arbitrary (VectorTriple n e) where
    arbitrary = sized $ \m -> do
        n <- choose (0,m)
        x <- vector n
        y <- vector n
        z <- vector n
        return $ VectorTriple x y z
    
    coarbitrary = undefined
