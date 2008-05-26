-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Complex
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Complex
    where

import Test.QuickCheck 

import Data.Complex

newtype TestComplex e = TestComplex (Complex e)

instance (Arbitrary e, RealFloat e) => Arbitrary (TestComplex e) where
    arbitrary = do
        (x,y) <- arbitrary
        return $ TestComplex (x :+ y)
    coarbitrary (TestComplex (x :+ y)) = 
        coarbitrary (x,y)
