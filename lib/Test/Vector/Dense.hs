{-# OPTIONS_GHC -fno-warn-orphans #-}
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
import qualified Test.QuickCheck as QC

import Data.Vector.Dense hiding ( vector )
import Data.Elem.BLAS ( BLAS1 )

data SubVector n e = 
    SubVector Int 
              (Vector n e) 
              Int 
              Int 
    deriving (Show)
    
data VectorPair n e = 
    VectorPair (Vector n e) 
               (Vector n e) 
    deriving (Show)
    
data VectorTriple n e = 
    VectorTriple (Vector n e) 
                 (Vector n e) 
                 (Vector n e) 
    deriving (Show)

vector :: (BLAS1 e, Arbitrary e) => Int -> Gen (Vector n e)
vector n =
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
    x <- vector n
    return $ (conj x)

subVector :: (BLAS1 e, Arbitrary e) => Int -> Gen (SubVector n e)
subVector n = do
    o <- choose (0,5)
    s <- choose (1,5)
    e <- choose (0,5)
    x <- vector (o + s*n + e)
    return (SubVector s x o n)

instance (Arbitrary e, BLAS1 e) => Arbitrary (Vector n e) where
    arbitrary = sized $ \m ->
        choose (0,m) >>= vector
        
    coarbitrary = undefined

instance (Arbitrary e, BLAS1 e) => Arbitrary (SubVector n e) where
    arbitrary = sized $ \m -> 
        choose (0,m) >>= subVector
        
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
