-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Vector.Dense
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Vector.Dense (
    TestVector(..),
    SubVector(..),
    VectorPair(..),
    VectorTriple(..),
    
    dvector,
    rawVector,
    conjVector,
    subVector
    ) where

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Data.Vector.Dense
import BLAS.Elem ( BLAS1 )

newtype TestVector n e = TestVector (Vector n e) deriving (Eq, Show)
data SubVector n e = SubVector Int (Vector n e) Int Int deriving (Eq, Show)
data VectorPair n e = Pair (Vector n e) (Vector n e) deriving (Eq, Show)
data VectorTriple n e = Triple (Vector n e) (Vector n e) (Vector n e) deriving (Eq, Show)


dvector :: (BLAS1 e, Arbitrary e) => Int -> Gen (Vector n e)
dvector n =
    frequency [ (3, rawVector n)  
              , (2, conjVector n)
              , (1, subVector n    >>= \(SubVector s x o _) -> 
                    return $ subvectorWithStride s x o n)
              ]    

rawVector :: (BLAS1 e, Arbitrary e) => Int -> Gen (Vector n e)
rawVector n = do
    es <- QC.vector n
    return $ listVector n es

conjVector :: (BLAS1 e, Arbitrary e) => Int -> Gen (Vector n e)
conjVector n = do
    x <- dvector n
    return $ (conj x)

subVector :: (BLAS1 e, Arbitrary e) => Int -> Gen (SubVector n e)
subVector n = do
    o <- choose (0,5)
    s <- choose (1,5)
    e <- choose (0,5)
    x <- dvector (o + s*n + e)
    return (SubVector s x o n)

instance (Arbitrary e, BLAS1 e) => Arbitrary (TestVector n e) where
    arbitrary = sized $ \m ->
        choose (0,m) >>= dvector >>= return . TestVector
    coarbitrary (TestVector x) =
        coarbitrary (elems x)

instance (Arbitrary e, BLAS1 e) => Arbitrary (SubVector n e) where
    arbitrary = sized $ \m -> 
        choose (0,m) >>= subVector
    coarbitrary (SubVector s x o n) = 
        coarbitrary (s,TestVector x,o,n)

instance (Arbitrary e, BLAS1 e) => Arbitrary (VectorPair n e) where
    arbitrary = sized $ \m -> do
        n <- choose (0,m)
        x <- dvector n
        y <- dvector n
        return $ Pair x y
        
    coarbitrary (Pair x y) = 
        coarbitrary (TestVector x, TestVector y)
        
instance (Arbitrary e, BLAS1 e) => Arbitrary (VectorTriple n e) where
    arbitrary = sized $ \m -> do
        n <- choose (0,m)
        x <- dvector n
        y <- dvector n
        z <- dvector n
        return $ Triple x y z
        
    coarbitrary (Triple x y z) = 
        coarbitrary (TestVector x, TestVector y, TestVector z)
        