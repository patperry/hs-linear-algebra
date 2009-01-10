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

import Test.QuickCheck hiding ( Test.vector )
import Test.QuickCheck.BLAS ( TestElem )
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense ( Vector, dim )
import Data.Matrix.Dense hiding ( Test.matrix )
import Data.Elem.BLAS ( BLAS3 )

matrix :: (TestElem e) => (Int,Int) -> Gen (Matrix (m,n) e)
matrix = Test.matrix


instance (TestElem e) => Arbitrary (Matrix (m,n) e) where
    arbitrary   = Test.matrix =<< Test.shape
    coarbitrary = undefined

data SubMatrix m n e = 
    SubMatrix (Matrix (m,n) e) 
              (Int,Int) 
              (Int,Int) 
    deriving (Show)

instance (TestElem e) => Arbitrary (SubMatrix m n e) where
    arbitrary = do
        (m,n) <- Test.shape
        i <- choose (0,5)
        j <- choose (0,5)
        e <- choose (0,5)
        f <- choose (0,5)
        a <- Test.matrix (i+m+e, j+n+f)
        return $ SubMatrix a (i,j) (m,n)
        
    coarbitrary = undefined


data MatrixAt m n e =
    MatrixAt (Matrix (m,n) e)
             (Int,Int)
    deriving (Show)
instance (TestElem e) => Arbitrary (MatrixAt m n e) where
    arbitrary = do
        (m',n') <- Test.shape
        i <- choose (0,m')
        j <- choose (0,n')
        a <- Test.matrix (m'+1,n'+1)
        return $ MatrixAt a (i,j)
        
    coarbitrary = undefined



        
data MatrixPair m n e = 
    MatrixPair (Matrix (m,n) e) 
               (Matrix (m,n) e) 
    deriving (Show)

instance (TestElem e) => Arbitrary (MatrixPair m n e) where
    arbitrary = do
        a <- arbitrary
        b <- Test.matrix (shape a)
        return $ MatrixPair a b
        
    coarbitrary = undefined
  
data MatrixMV m n e = 
    MatrixMV (Matrix (m,n) e) 
             (Vector n e) 
    deriving (Show)

instance (TestElem e) => Arbitrary (MatrixMV m n e) where
    arbitrary = do
        a <- arbitrary
        x <- Test.vector (numCols a)
        return $ MatrixMV a x
            
    coarbitrary = undefined

data MatrixMVPair m n e = 
    MatrixMVPair (Matrix (m,n) e) 
                 (Vector n e) 
                 (Vector n e) 
    deriving (Show)
    
instance (TestElem e) => Arbitrary (MatrixMVPair m n e) where
    arbitrary = do
        (MatrixMV a x) <- arbitrary
        y <- Test.vector (dim x)
        return $ MatrixMVPair a x y
        
    coarbitrary = undefined

        
data MatrixMM m n k e = 
    MatrixMM (Matrix (m,k) e) 
             (Matrix (k,n) e) 
    deriving (Show)

instance (TestElem e) => Arbitrary (MatrixMM m n k e) where
    arbitrary = do
        a <- arbitrary
        (_,n) <- Test.shape
        b <- Test.matrix (numCols a, n)
        return $ MatrixMM a b
            
    coarbitrary = undefined
        
data MatrixMMPair m n k e = 
    MatrixMMPair (Matrix (m,k) e) 
                 (Matrix (k,n) e) 
                 (Matrix (k,n) e)
    deriving (Show)
    
instance (TestElem e) => Arbitrary (MatrixMMPair m n k e) where
    arbitrary = do
        (MatrixMM a b) <- arbitrary
        c <- Test.matrix (shape b)
        return $ MatrixMMPair a b c
        
    coarbitrary = undefined
