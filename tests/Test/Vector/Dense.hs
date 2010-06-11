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

vector :: (TestElem e) => Int -> Gen (Vector e)
vector = Test.vector
    
data VectorPair e = 
    VectorPair (Vector e) 
               (Vector e) 
    deriving (Show)
    
data VectorTriple e = 
    VectorTriple (Vector e) 
                 (Vector e) 
                 (Vector e) 
    deriving (Show)

instance (TestElem e) => Arbitrary (Vector e) where
    arbitrary = do
        n <- Test.dim
        Test.vector n

instance (TestElem e) => Arbitrary (VectorPair e) where
    arbitrary = do
        n <- Test.dim
        x <- vector n
        y <- vector n
        return $ VectorPair x y
        
instance (TestElem e) => Arbitrary (VectorTriple e) where
    arbitrary = do
        n <- Test.dim
        x <- vector n
        y <- vector n
        z <- vector n
        return $ VectorTriple x y z

data SubVector e = 
    SubVector Int 
              (Vector e) 
              Int 
              Int 
    deriving (Show)

instance (TestElem e) => Arbitrary (SubVector e) where
    arbitrary = do
        n <- Test.dim
        o <- choose (0,5)
        s <- choose (1,5)
        e <- choose (0,5)
        x <- Test.vector (o + s*n + e)
        return (SubVector s x o n)
