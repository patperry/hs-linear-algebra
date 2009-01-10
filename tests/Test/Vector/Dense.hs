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
import Test.QuickCheck.BLAS ( TestElem )
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense hiding ( vector )
import Data.Elem.BLAS ( BLAS1 )

vector :: (TestElem e) => Int -> Gen (Vector n e)
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

instance (TestElem e) => Arbitrary (Vector n e) where
    arbitrary = do
        n <- Test.dim
        Test.vector n
        
    coarbitrary = undefined

instance (TestElem e) => Arbitrary (VectorPair n e) where
    arbitrary = do
        n <- Test.dim
        x <- vector n
        y <- vector n
        return $ VectorPair x y
        
    coarbitrary = undefined
        
instance (TestElem e) => Arbitrary (VectorTriple n e) where
    arbitrary = do
        n <- Test.dim
        x <- vector n
        y <- vector n
        z <- vector n
        return $ VectorTriple x y z
    
    coarbitrary = undefined

data SubVector n e = 
    SubVector Int 
              (Vector n e) 
              Int 
              Int 
    deriving (Show)

instance (TestElem e) => Arbitrary (SubVector n e) where
    arbitrary = do
        n <- Test.dim
        o <- choose (0,5)
        s <- choose (1,5)
        e <- choose (0,5)
        x <- Test.vector (o + s*n + e)
        return (SubVector s x o n)
        
    coarbitrary = undefined
