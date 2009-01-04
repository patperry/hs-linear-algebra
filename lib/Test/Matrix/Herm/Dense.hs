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

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Test.Vector.Dense ( vector )
import Test.Matrix ( matrixSized )
import Test.Matrix.Dense ( matrix )

import Data.Elem.BLAS ( BLAS3 )
import BLAS.Types ( flipUpLo )

import Data.Vector.Dense hiding ( vector )
import Data.Matrix.Dense hiding ( matrix )
import Data.Matrix.Herm



hermMatrix :: (BLAS3 e, Arbitrary e) => Int -> Gen (Matrix (n,n) e)
hermMatrix n  = do
    a <- matrix (n,n)
    let h = (a + herm a)
    elements [ h, herm h ]


data HermMatrix n e = 
    HermMatrix (Herm Matrix (n,n) e)
               (Matrix (n,n) e)
    deriving Show

instance (Arbitrary e, BLAS3 e) => Arbitrary (HermMatrix n e) where
    arbitrary = matrixSized $ \k -> do
        n <- choose (0,k)
        a <- hermMatrix n
        
        junk <- QC.vector (n*n)
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
        x <- vector (numCols a)
        return $ HermMatrixMV h a x
        
    coarbitrary = undefined

    
data HermMatrixMM m n e = 
    HermMatrixMM (Herm Matrix (m,m) e) 
                 (Matrix (m,m) e) 
                 (Matrix (m,n) e) 
    deriving Show
    
instance (Arbitrary e, BLAS3 e) => Arbitrary (HermMatrixMM m n e) where
    arbitrary = matrixSized $ \k -> do
        (HermMatrix h a) <- arbitrary
        n <- choose (0,k)
        b <- matrix (numCols a,n)
        return $ HermMatrixMM h a b
            
    coarbitrary = undefined
        