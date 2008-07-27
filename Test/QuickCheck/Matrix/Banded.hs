{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Matrix.Banded
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Matrix.Banded
    where

import Control.Monad( forM )

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Test.QuickCheck.Vector.Dense ( TestVector(..), dvector )
import Test.QuickCheck.Matrix.Dense ( dmatrix )

import Data.Vector.Dense ( Vector, dim )
import Data.Matrix.Dense ( Matrix )
import Data.Matrix.Banded
import BLAS.Elem ( Elem, BLAS1 )

newtype TestBanded m n e = TestBanded (Banded (m,n) e) deriving (Eq, Show)

bmatrix :: (BLAS1 e, Arbitrary e) => 
    (Int,Int) -> (Int,Int) -> Gen (Banded (m,n) e)
bmatrix mn kl = frequency [ (3, rawBanded mn kl)  
                          , (2, hermedBanded mn kl)
                          ]

rawBanded :: (BLAS1 e, Arbitrary e) => 
    (Int,Int) -> (Int,Int) -> Gen (Banded (m,n) e)
rawBanded (m,n) (k,l) = do
    xs <- QC.vector ((k+1+l)*len)
    return $ listsBanded (m,n) (k,l) (splitDiags xs)
  where
    len = min m n
    
    splitDiags [] = [[]]
    splitDiags xs = let (d,xs') = splitAt len xs
                    in d:(splitDiags xs')


hermedBanded :: (BLAS1 e, Arbitrary e) => 
    (Int,Int) -> (Int,Int) -> Gen (Banded (m,n) e)
hermedBanded (m,n) (l,u) = do
    x <- rawBanded (n,m) (u,l)
    return $ (herm x)

instance (Arbitrary e, BLAS1 e) => Arbitrary (TestBanded m n e) where
    arbitrary = sized $ \k -> 
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            m <- choose (0,k')
            n <- choose (0,k')
            kl <- if m == 0 then return 0 else choose (0,m-1)
            ku <- if n == 0 then return 0 else choose (0,n-1)
            bmatrix (m,n) (kl,ku) >>= return . TestBanded
        
    coarbitrary (TestBanded x) =
        coarbitrary (elems x)

data BandedAt m n e = BandedAt (Banded (m,n) e) (Int,Int) deriving (Eq, Show)
instance (Arbitrary e, BLAS1 e) => Arbitrary (BandedAt m n e) where
    arbitrary = sized $ \k ->
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            m  <- choose (1,k'+1)
            n  <- choose (1,k'+1)
            kl <- choose (0,m-1)
            ku <- choose (0,n-1)
            i  <- choose (0,m-1)
            j  <- choose (0,n-1)
            a  <- bmatrix (m,n) (kl,ku)
            
            return $ BandedAt a (i,j)

    coarbitrary = undefined


data ListsBanded e = ListsBanded (Int,Int) (Int,Int) [[e]] deriving (Eq,Show)
instance (Arbitrary e, Elem e) => Arbitrary (ListsBanded e) where
    arbitrary = sized $ \k ->
       let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
       in do
           m <- choose (1,k'+1)
           n <- choose (1,k'+1)
           kl <- choose (0,m-1)
           ku <- choose (0,n-1)
          
           ds <- forM [(-kl)..ku] $ \i ->
                     let beginPad = max (-i)    0
                         endPad   = max (m-n+i) 0
                         len      = m - (beginPad+endPad)
                     in do
                         xs <- QC.vector len
                         return $ replicate beginPad 0 ++ xs ++ replicate endPad 0
         
           return $ ListsBanded (m,n) (kl,ku) ds
          
    coarbitrary (ListsBanded mn lu ds) = coarbitrary (mn,lu,ds)
   
{-
                    
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
-}  
  
data BandedMV m n e = BandedMV (Banded (m,n) e) (Vector n e) deriving (Eq, Show)

instance (Arbitrary e, BLAS1 e) => Arbitrary (BandedMV m n e) where
    arbitrary = sized $ \k -> 
        let k' = ceiling (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            m <- choose (0,k')
            n <- choose (0,k')
            kl <- if m == 0 then return 0 else choose (0,m-1)
            ku <- if n == 0 then return 0 else choose (0,n-1)
            a <- bmatrix (m,n) (kl,ku)             
            x <- dvector n
            return $ BandedMV a x
            
    coarbitrary = undefined

data BandedMVPair m n e = BandedMVPair (Banded (m,n) e) (Vector n e) (Vector n e) 
    deriving (Eq, Show)
    
instance (Arbitrary e, BLAS1 e) => Arbitrary (BandedMVPair m n e) where
    arbitrary = do
        (BandedMV a x) <- arbitrary
        y <- dvector (dim x)
        return $ BandedMVPair a x y
        
    coarbitrary (BandedMVPair a x y) =
        coarbitrary (BandedMV a x, TestVector y)
        
data BandedMM m n k e = BandedMM (Banded (m,k) e) (Matrix (k,n) e) deriving (Eq, Show)

instance (Arbitrary e, BLAS1 e) => Arbitrary (BandedMM m n k e) where
    arbitrary = sized $ \s ->
        let s' = ceiling (sqrt $ fromInteger $ toInteger s :: Double)
        in do
            m <- choose (0,s')
            k <- choose (0,s')
            n <- choose (0,s')
            kl <- if m == 0 then return 0 else choose (0,m-1)
            ku <- if k == 0 then return 0 else choose (0,k-1)
            a <- bmatrix (m,k) (kl,ku)             
            b <- dmatrix (k,n)
            return $ BandedMM a b
            
    coarbitrary = undefined
        
data BandedMMPair m n k e = BandedMMPair (Banded (m,k) e) (Matrix (k,n) e) (Matrix (k,n) e)
    deriving (Eq, Show)
    
instance (Arbitrary e, BLAS1 e) => Arbitrary (BandedMMPair m n k e) where
    arbitrary = do
        (BandedMM a b) <- arbitrary
        c <- dmatrix (shape b)
        return $ BandedMMPair a b c
        
    coarbitrary = undefined
