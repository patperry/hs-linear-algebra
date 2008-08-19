-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Vector.Dense
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Generators.Vector.Dense (
    SubVector(..),
    VectorPair(..),
    VectorTriple(..),
    
    vector,
    rawVector,
    conjVector,
    subVector
    ) where

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Data.Vector.Dense hiding ( vector )
import BLAS.Elem ( Elem, BLAS1 )

data SubVector n e = SubVector Int (Vector n e) Int Int deriving (Show)
data VectorPair n e = Pair (Vector n e) (Vector n e) deriving (Show)
data VectorTriple n e = Triple (Vector n e) (Vector n e) (Vector n e) deriving (Show)


vector :: (Elem e, Arbitrary e) => Int -> Gen (Vector n e)
vector n =
    frequency [ (3, rawVector n)  
              , (2, conjVector n)
              , (1, subVector n    >>= \(SubVector s x o _) -> 
                    return $ subvectorWithStride s x o n)
              ]    

rawVector :: (Elem e, Arbitrary e) => Int -> Gen (Vector n e)
rawVector n = do
    es <- QC.vector n
    return $ listVector n es

conjVector :: (Elem e, Arbitrary e) => Int -> Gen (Vector n e)
conjVector n = do
    x <- vector n
    return $ (conj x)

subVector :: (Elem e, Arbitrary e) => Int -> Gen (SubVector n e)
subVector n = do
    o <- choose (0,5)
    s <- choose (1,5)
    e <- choose (0,5)
    x <- vector (o + s*n + e)
    return (SubVector s x o n)

instance (Arbitrary e, BLAS1 e) => Arbitrary (Vector n e) where
    arbitrary = sized $ \m ->
        choose (0,m) >>= vector
    coarbitrary x =
        coarbitrary (elems x)

instance (Arbitrary e, BLAS1 e) => Arbitrary (SubVector n e) where
    arbitrary = sized $ \m -> 
        choose (0,m) >>= subVector
    coarbitrary (SubVector s x o n) = 
        coarbitrary (s,x,o,n)

instance (Arbitrary e, BLAS1 e) => Arbitrary (VectorPair n e) where
    arbitrary = sized $ \m -> do
        n <- choose (0,m)
        x <- vector n
        y <- vector n
        return $ Pair x y
        
    coarbitrary (Pair x y) = 
        coarbitrary (x,y)
        
instance (Arbitrary e, BLAS1 e) => Arbitrary (VectorTriple n e) where
    arbitrary = sized $ \m -> do
        n <- choose (0,m)
        x <- vector n
        y <- vector n
        z <- vector n
        return $ Triple x y z
        
    coarbitrary (Triple x y z) = 
        coarbitrary (x,y,z)
        