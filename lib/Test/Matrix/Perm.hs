-----------------------------------------------------------------------------
-- |
-- Module     : Test.Matrix.Perm
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.Matrix.Perm
    where

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Test.Permutation
import Test.Vector.Dense ( vector )
import Test.Matrix.Dense ( matrix )
import Test.Matrix ( matrixSized )

import BLAS.Elem ( Elem, BLAS1, BLAS3 )

import Data.Vector.Dense ( Vector )
import Data.Matrix.Dense ( Matrix, shape )
import Data.Matrix.Perm

pmatrix :: (Elem e) => Int -> Gen (Perm (n,n) e)
pmatrix n = frequency [ (2, rawPerm n)
                      , (2, hermedPerm n)
                      , (1, idPerm n)
                      ]

rawPerm :: Int -> Gen (Perm (n,n) e)
rawPerm n = permutation n >>= return . fromPermutation

hermedPerm :: (Elem e) => Int -> Gen (Perm (n,n) e)
hermedPerm n = rawPerm n >>= return . herm

idPerm :: Int -> Gen (Perm (n,n) e)
idPerm = return . identity


newtype TestPerm n e = TestPerm (Perm (n,n) e) deriving Show

instance (Elem e) => Arbitrary (TestPerm n e) where
    arbitrary = matrixSized $ \s -> do
        n <- choose (0,s)
        p <- pmatrix n
        return $ TestPerm p
        
    coarbitrary = undefined

data PermMBasis n e = PermMBasis (Perm (n,n) e) Int deriving Show

instance (Elem e) => Arbitrary (PermMBasis n e) where
    arbitrary = do
        (TestPerm p) <- arbitrary
        i <- choose (0, numCols p - 1)
        return $ PermMBasis p i
        
    coarbitrary = undefined
    
data PermMV n e = PermMV (Perm (n,n) e) (Vector n e) deriving Show

instance (BLAS1 e, Arbitrary e) => Arbitrary (PermMV n e) where
    arbitrary = do
        (TestPerm p) <- arbitrary
        x <- vector (numCols p)
        return $ PermMV p x
    
    coarbitrary = undefined
    
data PermMVPair n e = 
    PermMVPair (Perm (n,n) e) (Vector n e) (Vector n e) deriving Show
    
instance (BLAS1 e, Arbitrary e) => Arbitrary (PermMVPair n e) where
    arbitrary = do
        (PermMV p x) <- arbitrary
        y <- vector (numCols p)
        return $ PermMVPair p x y
        
    coarbitrary = undefined
    
data PermMM m n e = PermMM (Perm (m,m) e) (Matrix (m,n) e) deriving Show

instance (BLAS3 e, Arbitrary e) => Arbitrary (PermMM m n e) where
    arbitrary = matrixSized $ \s -> do
        (TestPerm p) <- arbitrary
        n <- choose (0,s)
        a <- matrix (numCols p, n)
        return $ PermMM p a
        
    coarbitrary = undefined
    
data PermMMPair m n e = 
    PermMMPair (Perm (m,m) e) (Matrix (m,n) e) (Matrix (m,n) e) deriving Show
    
instance (BLAS3 e, Arbitrary e) => Arbitrary (PermMMPair m n e) where
    arbitrary = do
        (PermMM p a)<- arbitrary
        b <- matrix (shape a)
        return $ PermMMPair p a b
        
    coarbitrary = undefined
