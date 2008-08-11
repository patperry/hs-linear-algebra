-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Matrix.Diag
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Matrix.Diag
    where

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Test.QuickCheck.Vector.Dense ( dvector )
import Test.QuickCheck.Matrix.Dense ( dmatrix )
import Test.QuickCheck.Matrix ( matrixSized )

import BLAS.Elem ( BLAS1, BLAS2 )

import Data.Vector.Dense ( Vector )
import Data.Matrix.Dense ( Matrix, matrix, diag )
import Data.Matrix.Diag 

diagMatrix :: (Arbitrary e, BLAS1 e) => Int -> Gen (Matrix (n,n) e)
diagMatrix n = do
    es <- QC.vector n
    let a = matrix (n,n) [ ((i,i),e) | (i,e) <- zip [0..(n-1)] es ]
    elements [ a, herm a ]


data TestDiag n e = 
    TestDiag (Diag (n,n) e) 
             (Matrix (n,n) e) 
    deriving Show

instance (Arbitrary e, BLAS1 e) => Arbitrary (TestDiag n e) where
    arbitrary = matrixSized $ \s -> do
        n <- choose (0,s)
        a <- diagMatrix n
        return $ TestDiag (fromVector $ diag a 0) a
        
    coarbitrary = undefined


data DiagMV n e = 
    DiagMV (Diag (n,n) e)
           (Matrix (n,n) e)
           (Vector n e)
    deriving Show
    
instance (Arbitrary e, BLAS1 e) => Arbitrary (DiagMV n e) where
    arbitrary = do
        (TestDiag d a) <- arbitrary
        x <- dvector (numCols a)
        return $ DiagMV d a x

    coarbitrary = undefined


data DiagMM m n e = 
    DiagMM (Diag (m,m) e)
           (Matrix (m,m) e)
           (Matrix (m,n) e)
    deriving Show

instance (Arbitrary e, BLAS1 e) => Arbitrary (DiagMM m n e) where
    arbitrary = matrixSized $ \s -> do
        (TestDiag d a) <- arbitrary
        n <- choose (0,s)
        b <- dmatrix (numCols a, n)
        return $ DiagMM d a b

    coarbitrary = undefined
    
    
data DiagSV n e =
    DiagSV (Diag (n,n) e)
           (Vector n e)
    deriving Show

instance (Arbitrary e, BLAS2 e) => Arbitrary (DiagSV n e) where
    arbitrary = do
        (DiagMV d _ x) <- arbitrary
        return $ DiagSV d (d <*> x)
        
    coarbitrary = undefined

data DiagSM m n e =
    DiagSM (Diag (m,m) e)
           (Matrix (m,n) e)
    deriving Show
    
instance (Arbitrary e, BLAS2 e) => Arbitrary (DiagSM m n e) where
    arbitrary = do
        (DiagMM d _ b) <- arbitrary
        return $ DiagSM d (d <**> b)

    coarbitrary = undefined
