{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Matrix.Dense
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Matrix.Dense
    where

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Test.QuickCheck.Vector.Dense ( TestVector(..), dvector )

import Data.Vector.Dense ( Vector, dim )
import Data.Matrix.Dense
import BLAS.Elem ( Elem, BLAS1 )

newtype TestMatrix m n e = TestMatrix (Matrix (m,n) e) deriving (Eq, Show)
data SubMatrix m n e = SubMatrix (Matrix (m,n) e) (Int,Int) (Int,Int) deriving (Eq, Show)

dmatrix :: (Elem e, Arbitrary e) => (Int,Int) -> Gen (Matrix (m,n) e)
dmatrix mn = frequency [ (3, rawMatrix mn)  
                       , (2, hermedMatrix mn)
                       , (1, subMatrix mn 
                             >>= \(SubMatrix a ij _) -> 
                                 return $ submatrix a ij mn)
                       ]

rawMatrix :: (Elem e, Arbitrary e) => (Int,Int) -> Gen (Matrix (m,n) e)
rawMatrix (m,n) = do
    es <- QC.vector (m*n)
    return $ listMatrix (m,n) es

hermedMatrix :: (Elem e, Arbitrary e) => (Int,Int) -> Gen (Matrix (m,n) e)
hermedMatrix (m,n) = do
    x <- dmatrix (n,m)
    return $ (herm x)

subMatrix :: (Elem e, Arbitrary e) => (Int,Int) -> Gen (SubMatrix m n e)
subMatrix (m,n) = do
    i <- choose (0,5)
    j <- choose (0,5)
    e <- choose (0,5)
    f <- choose (0,5)
    x <- dmatrix (i+m+e, j+n+f)

    return $ SubMatrix x (i,j) (m,n)



instance (Arbitrary e, BLAS1 e) => Arbitrary (TestMatrix m n e) where
    arbitrary = sized $ \k -> 
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            m <- choose (0,k')
            n <- choose (0,k')
            dmatrix (m,n) >>= return . TestMatrix
        
    coarbitrary (TestMatrix x) =
        coarbitrary (elems x)

data MatAt m n e = MatAt (Matrix (m,n) e) (Int,Int) deriving (Eq, Show)
instance (Arbitrary e, Elem e) => Arbitrary (MatAt m n e) where
    arbitrary = sized $ \k ->
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            m <- choose (1,k'+1)
            n <- choose (1,k'+1)
            i <- choose (0,m-1)
            j <- choose (0,n-1)
            a <- dmatrix (m,n)
            return $ MatAt a (i,j)

    coarbitrary = undefined
            

instance (Arbitrary e, BLAS1 e) => Arbitrary (SubMatrix m n e) where
    arbitrary = sized $ \k -> 
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double) 
        in do
            m <- choose (0,k')
            n <- choose (0,k')
            (SubMatrix a ij mn) <- subMatrix (m,n)
            return $ SubMatrix a ij mn
        
    coarbitrary (SubMatrix a ij mn) =
        coarbitrary (elems a, ij, mn)
        
data MatrixPair m n e = Pair (Matrix (m,n) e) (Matrix (m,n) e) deriving (Eq, Show)

instance (Arbitrary e, Elem e) => Arbitrary (MatrixPair m n e) where
    arbitrary = sized $ \k -> 
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            m <- choose (0,k')
            n <- choose (0,k')
            a <- dmatrix (m,n)
            b <- dmatrix (m,n)
            return $ Pair a b
        
    coarbitrary = undefined
  
data MultMV m n e = MultMV (Matrix (m,n) e) (Vector n e) deriving (Eq, Show)

instance (Arbitrary e, BLAS1 e) => Arbitrary (MultMV m n e) where
    arbitrary = sized $ \k -> 
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            m <- choose (0,k')
            n <- choose (0,k')
            a <- dmatrix (m,n)
            x <- dvector n
            return $ MultMV a x
            
    coarbitrary = undefined

data MultMVPair m n e = MultMVPair (Matrix (m,n) e) (Vector n e) (Vector n e) 
    deriving (Eq, Show)
    
instance (Arbitrary e, BLAS1 e) => Arbitrary (MultMVPair m n e) where
    arbitrary = do
        (MultMV a x) <- arbitrary
        y <- dvector (dim x)
        return $ MultMVPair a x y
        
    coarbitrary (MultMVPair a x y) =
        coarbitrary (MultMV a x, TestVector y)
        
data MultMM m n k e = MultMM (Matrix (m,k) e) (Matrix (k,n) e) deriving (Eq, Show)

instance (Arbitrary e, Elem e) => Arbitrary (MultMM m n k e) where
    arbitrary = sized $ \s ->
        let s' = ceiling (sqrt $ fromInteger $ toInteger s :: Double)
        in do
            m <- choose (0,s')
            n <- choose (0,s')
            k <- choose (0,s')
            a <- dmatrix (m,k)
            b <- dmatrix (k,n)
            return $ MultMM a b
            
    coarbitrary = undefined
        
data MultMMPair m n k e = MultMMPair (Matrix (m,k) e) (Matrix (k,n) e) (Matrix (k,n) e)
    deriving (Eq, Show)
    
instance (Arbitrary e, Elem e) => Arbitrary (MultMMPair m n k e) where
    arbitrary = do
        (MultMM a b) <- arbitrary
        c <- dmatrix (shape b)
        return $ MultMMPair a b c
        
    coarbitrary = undefined
