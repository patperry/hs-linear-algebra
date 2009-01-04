{-# LANGUAGE FlexibleInstances #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.Matrix.Dense
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.Matrix.Dense (
    matrix,
    
    MatrixAt(..),
    SubMatrix(..),
    MatrixPair(..),
    MatrixMV(..),
    MatrixMVPair(..),
    MatrixMM(..),
    MatrixMMPair(..)
    ) where

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Test.Vector.Dense ( vector )
import Test.Matrix

import Data.Vector.Dense ( Vector, dim )
import Data.Matrix.Dense hiding ( matrix )
import Data.Elem.BLAS ( BLAS3 )

data SubMatrix m n e = 
    SubMatrix (Matrix (m,n) e) 
              (Int,Int) 
              (Int,Int) 
    deriving (Show)

matrix :: (BLAS3 e, Arbitrary e) => (Int,Int) -> Gen (Matrix (m,n) e)
matrix mn = frequency [ (3, rawMatrix mn)  
                      , (2, hermedMatrix mn)
                      , (1, subMatrix mn >>= \(SubMatrix a ij _) -> 
                                 return $ submatrix a ij mn)
                      ]

rawMatrix :: (BLAS3 e, Arbitrary e) => (Int,Int) -> Gen (Matrix (m,n) e)
rawMatrix (m,n) = do
    es <- QC.vector (m*n)
    return $ listMatrix (m,n) es

hermedMatrix :: (BLAS3 e, Arbitrary e) => (Int,Int) -> Gen (Matrix (m,n) e)
hermedMatrix (m,n) = do
    x <- matrix (n,m)
    return $ (herm x)

subMatrix :: (BLAS3 e, Arbitrary e) => (Int,Int) -> Gen (SubMatrix m n e)
subMatrix (m,n) = do
    i <- choose (0,5)
    j <- choose (0,5)
    e <- choose (0,5)
    f <- choose (0,5)
    x <- matrix (i+m+e, j+n+f)

    return $ SubMatrix x (i,j) (m,n)

instance (Arbitrary e, BLAS3 e) => Arbitrary (Matrix (m,n) e) where
    arbitrary = matrixSized $ \k -> do
        m <- choose (0,k)
        n <- choose (0,k)
        matrix (m,n)
        
    coarbitrary = undefined

data MatrixAt m n e =
    MatrixAt (Matrix (m,n) e)
             (Int,Int)
    deriving (Show)
instance (Arbitrary e, BLAS3 e) => Arbitrary (MatrixAt m n e) where
    arbitrary = matrixSized $ \k -> do
        m <- choose (0,k) >>= return . (+1)
        n <- choose (0,k) >>= return . (+1)
        i <- choose (0,m-1)
        j <- choose (0,n-1)
        a <- matrix (m,n)
        return $ MatrixAt a (i,j)
        
    coarbitrary = undefined


instance (Arbitrary e, BLAS3 e) => Arbitrary (SubMatrix m n e) where
    arbitrary = matrixSized $ \k -> do
        m <- choose (0,k)
        n <- choose (0,k)
        (SubMatrix a ij mn) <- subMatrix (m,n)
        return $ SubMatrix a ij mn
        
    coarbitrary = undefined
        
data MatrixPair m n e = 
    MatrixPair (Matrix (m,n) e) 
               (Matrix (m,n) e) 
    deriving (Show)

instance (Arbitrary e, BLAS3 e) => Arbitrary (MatrixPair m n e) where
    arbitrary = do
        a <- arbitrary
        b <- matrix (shape a)
        return $ MatrixPair a b
        
    coarbitrary = undefined
  
data MatrixMV m n e = 
    MatrixMV (Matrix (m,n) e) 
             (Vector n e) 
    deriving (Show)

instance (Arbitrary e, BLAS3 e) => Arbitrary (MatrixMV m n e) where
    arbitrary = do
        a <- arbitrary
        x <- vector (numCols a)
        return $ MatrixMV a x
            
    coarbitrary = undefined

data MatrixMVPair m n e = 
    MatrixMVPair (Matrix (m,n) e) 
                 (Vector n e) 
                 (Vector n e) 
    deriving (Show)
    
instance (Arbitrary e, BLAS3 e) => Arbitrary (MatrixMVPair m n e) where
    arbitrary = do
        (MatrixMV a x) <- arbitrary
        y <- vector (dim x)
        return $ MatrixMVPair a x y
        
    coarbitrary = undefined

        
data MatrixMM m n k e = 
    MatrixMM (Matrix (m,k) e) 
             (Matrix (k,n) e) 
    deriving (Show)

instance (Arbitrary e, BLAS3 e) => Arbitrary (MatrixMM m n k e) where
    arbitrary = matrixSized $ \s -> do
        a <- arbitrary
        n <- choose (0,s)
        b <- matrix (numCols a, n)
        return $ MatrixMM a b
            
    coarbitrary = undefined
        
data MatrixMMPair m n k e = 
    MatrixMMPair (Matrix (m,k) e) 
                 (Matrix (k,n) e) 
                 (Matrix (k,n) e)
    deriving (Show)
    
instance (Arbitrary e, BLAS3 e) => Arbitrary (MatrixMMPair m n k e) where
    arbitrary = do
        (MatrixMM a b) <- arbitrary
        c <- matrix (shape b)
        return $ MatrixMMPair a b c
        
    coarbitrary = undefined
