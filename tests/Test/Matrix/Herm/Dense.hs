{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.Matrix.Herm.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.Matrix.Herm.Dense( 
    HermMatrix(..),
    HermMatrixMV(..),
    HermMatrixMM(..),
    ) where

import Control.Monad( liftM )

import Test.QuickCheck hiding ( Test.vector )
import qualified Test.QuickCheck.BLAS as Test

import Data.Elem.BLAS ( BLAS3 )

import Data.Vector.Dense hiding ( Test.vector )
import Data.Matrix.Dense hiding ( Test.matrix )
import Data.Matrix.Herm



hermMatrix :: (BLAS3 e, Arbitrary e) => Int -> Gen (Matrix (n,n) e)
hermMatrix n  = do
    a <- Test.matrix (n,n)
    let h = (a + herm a)
    elements [ h, herm h ]


data HermMatrix n e = 
    HermMatrix (Herm Matrix (n,n) e)
               (Matrix (n,n) e)
    deriving Show

instance (Arbitrary e, BLAS3 e) => Arbitrary (HermMatrix n e) where
    arbitrary = do
        n <- liftM fst Test.shape
        a <- hermMatrix n
        
        junk <- Test.elements (n*n)
        let (u ,b ) = (Upper, a // zip (filter (uncurry (>)) $ indices a) junk)
            (u',b') = (Lower, a // zip (filter (uncurry (<)) $ indices a) junk)

        h <- elements [ hermFromBase u             b
                      , hermFromBase (flipUpLo u)  (herm b)
                      , hermFromBase u'            b'
                      , hermFromBase (flipUpLo u') (herm b')
                      ]
        return $ HermMatrix h a
        
    coarbitrary = undefined


data HermMatrixMV n e = 
    HermMatrixMV (Herm Matrix (n,n) e) 
                 (Matrix (n,n) e) 
                 (Vector n e) 
    deriving Show

instance (Arbitrary e, BLAS3 e) => Arbitrary (HermMatrixMV n e) where
    arbitrary = do
        (HermMatrix h a) <- arbitrary
        x <- Test.vector (numCols a)
        return $ HermMatrixMV h a x
        
    coarbitrary = undefined

    
data HermMatrixMM m n e = 
    HermMatrixMM (Herm Matrix (m,m) e) 
                 (Matrix (m,m) e) 
                 (Matrix (m,n) e) 
    deriving Show
    
instance (Arbitrary e, BLAS3 e) => Arbitrary (HermMatrixMM m n e) where
    arbitrary = do
        (HermMatrix h a) <- arbitrary
        n <- liftM fst Test.shape
        b <- Test.matrix (numCols a,n)
        return $ HermMatrixMM h a b
            
    coarbitrary = undefined
        