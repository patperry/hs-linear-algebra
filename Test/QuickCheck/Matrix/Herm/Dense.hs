-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Matrix.Herm.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Matrix.Herm.Dense 
    where

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Test.QuickCheck.Vector.Dense ( dvector )
import Test.QuickCheck.Matrix.Dense ( dmatrix, rawMatrix )

import BLAS.Elem ( BLAS2 )

import Data.Vector.Dense
import Data.Matrix.Dense
import Data.Matrix.Herm



hermMatrix :: (BLAS2 e, Arbitrary e) => Int -> Gen (Matrix (n,n) e)
hermMatrix n  = do
    a <- rawMatrix (n,n)
    return $ (a + herm a)

data HermMatrixMV n e = 
    HermMatrixMV UpLo (Matrix (n,n) e) (Matrix (n,n) e) (Vector n e) deriving (Show)

instance (Arbitrary e, BLAS2 e) => Arbitrary (HermMatrixMV n e) where
    arbitrary = sized $ \k ->
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            u <- elements [ Upper, Lower ]
            n <- choose (0,k')
            h <- hermMatrix n
            let f  = case u of
                         Upper -> \(i,j) -> i > j
                         Lower -> \(i,j) -> i < j
                zs = zip (filter f (indices h)) (repeat 0)
                a  = h // zs
                
            x <- dvector n
            return $ HermMatrixMV u a h x
    coarbitrary = undefined
    
data HermMatrixMM m n e = 
    HermMatrixMM UpLo (Matrix (m,m) e) (Matrix (m,m) e) (Matrix (m,n) e) deriving (Show)
    
instance (Arbitrary e, BLAS2 e) => Arbitrary (HermMatrixMM m n e) where
    arbitrary = sized $ \k ->
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            u <- elements [ Upper, Lower ]
            m <- choose (0,k')
            n <- choose (0,k')
            h <- hermMatrix m
            let f  = case u of
                         Upper -> \(i,j) -> i > j
                         Lower -> \(i,j) -> i < j
                zs = zip (filter f (indices h)) (repeat 0)
                a  = h // zs

            b <- dmatrix (m,n)
            return $ HermMatrixMM u a h b
            
    coarbitrary = undefined
        