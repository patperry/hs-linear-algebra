{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Matrix.Tri.Dense
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Matrix.Tri.Dense
    where

import BLAS.Access 
import Data.Matrix.Dense.Internal ( DMatrix )
import Data.Ix ( range )

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Test.QuickCheck.Vector.Dense ( TestVector(..), dvector )
import Test.QuickCheck.Matrix.Dense ( dmatrix )

import Data.Vector.Dense
import Data.Matrix.Dense
import BLAS.Elem ( BLAS1, BLAS3 )

import Data.Matrix.Tri.Dense ( Tri, UpLo(..), Diag(..), upper, lower, upperU, lowerU )

triMatrix :: (BLAS1 e, Arbitrary e) => UpLo -> Diag -> Int -> Gen (Matrix (n,n) e)
triMatrix u d n = 
    let nz = case d of
                 NonUnit -> n * (n + 1) `div` 2
                 Unit    -> n * (n - 1) `div` 2
    in do
        h <- arbitrary
        let f = case (h,u,d) of
                 (False, Upper, NonUnit) -> \(i,j) -> i <= j
                 (False, Upper,    Unit) -> \(i,j) -> i <  j
                 (False, Lower, NonUnit) -> \(i,j) -> i >= j
                 (False, Lower,    Unit) -> \(i,j) -> i >  j
                 (True,  Upper, NonUnit) -> \(i,j) -> i >= j
                 (True,  Upper,    Unit) -> \(i,j) -> i >  j
                 (True,  Lower, NonUnit) -> \(i,j) -> i <= j
                 (True,  Lower,    Unit) -> \(i,j) -> i <  j
            ijs = filter f $ range ((0,0), (n-1,n-1))
        es <- QC.vector nz
        let a = matrix (n,n) $ zip ijs es

        a' <- case h of
                  False -> return a
                  True  -> return (herm a)
        return a'


data TriMatrix n e = TriMatrix UpLo Diag (Matrix (n,n) e) deriving (Eq, Show)

instance (Arbitrary e, BLAS1 e) => Arbitrary (TriMatrix n e) where
    arbitrary = sized $ \k ->
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            u <- elements [ Upper, Lower  ]
            d <- elements [ Unit, NonUnit ]
            n <- choose (0,k')
            a <- triMatrix u d n
            return $ TriMatrix u d a
            
    coarbitrary = undefined

data TriMatrixMV n e = 
    TriMatrixMV UpLo Diag (Matrix (n,n) e) (Vector n e) deriving (Eq, Show)

instance (Arbitrary e, BLAS1 e) => Arbitrary (TriMatrixMV n e) where
    arbitrary = do
        (TriMatrix u d a) <- arbitrary
        x <- dvector (numCols a)
        return $ TriMatrixMV u d a x
        
    coarbitrary (TriMatrixMV u d a x) =
        coarbitrary (TriMatrix u d a, TestVector x)
        
data TriMatrixMM m n e = 
    TriMatrixMM UpLo Diag (Matrix (m,m) e) (Matrix (m,n) e) deriving (Eq, Show)

instance (Arbitrary e, BLAS1 e) => Arbitrary (TriMatrixMM m n e) where
    arbitrary = sized $ \k ->
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            (TriMatrix u d a) <- arbitrary
            n <- choose (0,k')
            b <- dmatrix (numCols a, n)
            return $ TriMatrixMM u d a b
            
    coarbitrary = undefined
        
data TriMatrixSV n e = 
    TriMatrixSV (Tri (DMatrix Imm) (n,n) e) (Vector n e) deriving (Show)
    
instance (Arbitrary e, BLAS3 e) => Arbitrary (TriMatrixSV n e) where
    arbitrary = do
        (TriMatrix u d a) <- arbitrary
        let t = case (u,d) of
                    (Lower,NonUnit) -> lower  a
                    (Lower,Unit)    -> lowerU a
                    (Upper,NonUnit) -> upper  a
                    (Upper,Unit)    -> upperU a
        k <- arbitrary
        t' <- elements [ t
                       , k *> t
                       ]
        x <- dvector (numCols t')
        let y = t' <*> x
        return (TriMatrixSV t' y)
        
    coarbitrary = undefined


data TriMatrixSM m n e = 
    TriMatrixSM (Tri (DMatrix Imm) (m,m) e) (Matrix (m,n) e) 
    deriving (Show)
    
instance (Arbitrary e, BLAS3 e) => Arbitrary (TriMatrixSM m n e) where
    arbitrary = sized $ \k ->
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            (TriMatrixSV t _) <- arbitrary
            n <- choose (0, k')
            a <- dmatrix (numCols t, n)
            
            let b = t <**> a
            return (TriMatrixSM t b)
        
    coarbitrary = undefined

    