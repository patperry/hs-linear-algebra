
-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Vector
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Vector
    where

import Data.List ( nub )
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC


newtype Index = Index Int deriving (Eq, Show)
instance Arbitrary Index where
    arbitrary = sized $ \n -> do
        i <- if n == 0 
                then return 0
                else choose (0,n-1)
        return (Index i)
        
    coarbitrary (Index i) = coarbitrary i
    
data Basis = Basis Int Int deriving (Eq, Show)
instance Arbitrary Basis where
    arbitrary = do
        n <- arbitrary >>= (\(Index x) -> return (x+1))
        i <- choose (0,n-1)
        return $ Basis n i

    coarbitrary (Basis n i) = coarbitrary (n,i)

data Assocs e = Assocs Int [(Int,e)] deriving (Eq, Show)
instance Arbitrary e => Arbitrary (Assocs e) where
    arbitrary = sized $ \n -> do
        m <- choose (0,n)
        is <- QC.vector m 
              >>= mapM (\(Index i) -> return i) 
              >>= return . nub 
              >>= return . filter (<n)
        es <- QC.vector m
        return $ Assocs n (zip is es)
        
    coarbitrary (Assocs n ies) = coarbitrary (n,ies)
    