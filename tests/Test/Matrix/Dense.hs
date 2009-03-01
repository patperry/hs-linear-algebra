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

matrix :: (TestElem e) => (Int,Int) -> Gen (Matrix e)
matrix = Test.matrix


instance (TestElem e) => Arbitrary (Matrix e) where
    arbitrary   = Test.matrix =<< Test.shape
    coarbitrary = undefined

data SubMatrix e = 
    SubMatrix (Matrix e) 
              (Int,Int) 
              (Int,Int) 
    deriving (Show)

instance (TestElem e) => Arbitrary (SubMatrix e) where
    arbitrary = do
        (m,n) <- Test.shape
        i <- choose (0,5)
        j <- choose (0,5)
        e <- choose (0,5)
        f <- choose (0,5)
        a <- Test.matrix (i+m+e, j+n+f)
        return $ SubMatrix a (i,j) (m,n)
        
    coarbitrary = undefined


data MatrixAt e =
    MatrixAt (Matrix e)
             (Int,Int)
    deriving (Show)
instance (TestElem e) => Arbitrary (MatrixAt e) where
    arbitrary = do
        (m',n') <- Test.shape
        i <- choose (0,m')
        j <- choose (0,n')
        a <- Test.matrix (m'+1,n'+1)
        return $ MatrixAt a (i,j)
        
    coarbitrary = undefined



        
data MatrixPair e = 
    MatrixPair (Matrix e) 
               (Matrix e) 
    deriving (Show)

instance (TestElem e) => Arbitrary (MatrixPair e) where
    arbitrary = do
        a <- arbitrary
        b <- Test.matrix (shape a)
        return $ MatrixPair a b
        
    coarbitrary = undefined
  
data MatrixMV e = 
    MatrixMV (Matrix e) 
             (Vector e) 
    deriving (Show)

instance (TestElem e) => Arbitrary (MatrixMV e) where
    arbitrary = do
        a <- arbitrary
        x <- Test.vector (numCols a)
        return $ MatrixMV a x
            
    coarbitrary = undefined

data MatrixMVPair e = 
    MatrixMVPair (Matrix e) 
                 (Vector e) 
                 (Vector e) 
    deriving (Show)
    
instance (TestElem e) => Arbitrary (MatrixMVPair e) where
    arbitrary = do
        (MatrixMV a x) <- arbitrary
        y <- Test.vector (dim x)
        return $ MatrixMVPair a x y
        
    coarbitrary = undefined

        
data MatrixMM e = 
    MatrixMM (Matrix e) 
             (Matrix e) 
    deriving (Show)

instance (TestElem e) => Arbitrary (MatrixMM e) where
    arbitrary = do
        a <- arbitrary
        (_,n) <- Test.shape
        b <- Test.matrix (numCols a, n)
        return $ MatrixMM a b
            
    coarbitrary = undefined
        
data MatrixMMPair e = 
    MatrixMMPair (Matrix e) 
                 (Matrix e) 
                 (Matrix e)
    deriving (Show)
    
instance (TestElem e) => Arbitrary (MatrixMMPair e) where
    arbitrary = do
        (MatrixMM a b) <- arbitrary
        c <- Test.matrix (shape b)
        return $ MatrixMMPair a b c
        
    coarbitrary = undefined
