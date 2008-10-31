{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Matrix.Tri.Banded
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Matrix.Tri.Banded
    where

import Control.Monad ( replicateM )

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Test.QuickCheck.Vector.Dense ( TestVector(..), dvector )
import Test.QuickCheck.Matrix.Dense ( dmatrix )
import Test.QuickCheck.Matrix ( matrixSized )

import Data.Vector.Dense ( Vector )
import Data.Matrix.Dense ( Matrix )
import Data.Matrix.Banded
import BLAS.Elem ( BLAS1, BLAS2 )

import Data.Matrix.Tri.Banded ( Tri, UpLo(..), Diag(..), fromBase )

triBanded :: (BLAS1 e, Arbitrary e) => UpLo -> Diag -> Int -> Int -> Gen (Banded (n,n) e)
triBanded Upper NonUnit n k = do
    a <- triBanded Upper Unit n k
    d <- QC.vector n
    let (_,_,(_:ds)) = listsFromBanded a
    return $ listsBanded (n,n) (0,k) (d:ds)

triBanded Lower NonUnit n k = do
    a <- triBanded Lower Unit n k
    d <- QC.vector n
    let (_,_,ds) = listsFromBanded a
        ds' = (init ds) ++ [d]
    return $ listsBanded (n,n) (k,0) ds'
    
triBanded _ Unit n 0 = do
    return $ listsBanded (n,n) (0,0) [replicate n 1]
    
triBanded Upper Unit n k = do
    a <- triBanded Upper Unit n (k-1)
    let (_,_,ds) = listsFromBanded a
    
    d <- QC.vector (n-k) >>= \xs -> return $ xs ++ replicate k 0
    
    return $ listsBanded (n,n) (0,k) $ ds ++ [d]
    
triBanded Lower Unit n k = do
    a <- triBanded Lower Unit n (k-1)
    let (_,_,ds) = listsFromBanded a

    d <- QC.vector (n-k) >>= \xs -> return $ replicate k 0 ++ xs
    
    return $ listsBanded (n,n) (k,0) $ [d] ++ ds
    
    
data TriBanded n e = 
    TriBanded (Tri Banded (n,n) e) (Banded (n,n) e) deriving Show

instance (Arbitrary e, BLAS1 e) => Arbitrary (TriBanded n e) where
    arbitrary = matrixSized $ \s -> do
        u <- elements [ Upper, Lower  ]
        d <- elements [ Unit, NonUnit ]
        n <- choose (0,s)
        k <- if n == 0 then return 0 else choose (0,n-1)
        a <- triBanded u d n k

        l    <- if n == 0 then return 0 else choose (0,n-1)
        junk <- replicateM l $ QC.vector n
        diagJunk <- QC.vector n
        let (_,_,ds) = listsFromBanded a
            t = case (u,d) of 
                    (Upper,NonUnit) -> 
                        listsBanded (n,n) (l,k) $ junk ++ ds
                    (Upper,Unit) ->
                        listsBanded (n,n) (l,k) $ junk ++ [diagJunk] ++ tail ds
                    (Lower,NonUnit) -> 
                        listsBanded (n,n) (k,l) $ ds ++ junk
                    (Lower,Unit) -> 
                        listsBanded (n,n) (k,l) $ init ds ++ [diagJunk] ++ junk
            t' = fromBase u d t 

        return $ TriBanded t' a
            
    coarbitrary = undefined

data TriBandedMV n e = 
    TriBandedMV (Tri Banded (n,n) e) (Banded (n,n) e) (Vector n e) deriving Show

instance (Arbitrary e, BLAS1 e) => Arbitrary (TriBandedMV n e) where
    arbitrary = do
        (TriBanded t a) <- arbitrary
        x <- dvector (numCols t)
        return $ TriBandedMV t a x
        
    coarbitrary (TriBandedMV t a x) =
        coarbitrary (TriBanded t a, TestVector x)
        
data TriBandedMM m n e = 
    TriBandedMM (Tri Banded (m,m) e) (Banded (m,m) e) (Matrix (m,n) e) deriving Show

instance (Arbitrary e, BLAS1 e) => Arbitrary (TriBandedMM m n e) where
    arbitrary = matrixSized $ \s -> do
        (TriBanded t a) <- arbitrary
        n <- choose (0,s)
        b <- dmatrix (numCols t, n)
        return $ TriBandedMM t a b
            
    coarbitrary = undefined
        
data TriBandedSV n e = 
    TriBandedSV (Tri Banded (n,n) e) (Vector n e) deriving (Show)
    
instance (Arbitrary e, BLAS2 e) => Arbitrary (TriBandedSV n e) where
    arbitrary = do
        (TriBanded t _) <- arbitrary
        t' <- elements [ t
                       , herm t
                       ]
        x <- dvector (numCols t')
        let y = t' <*> x
        return (TriBandedSV t' y)
        
    coarbitrary = undefined


data TriBandedSM m n e = 
    TriBandedSM (Tri Banded (m,m) e) (Matrix (m,n) e) 
    deriving (Show)
    
instance (Arbitrary e, BLAS2 e) => Arbitrary (TriBandedSM m n e) where
    arbitrary = matrixSized $ \s -> do
        (TriBandedSV t _) <- arbitrary
        t' <- elements [ t
                       , herm t
                       ]

        n <- choose (0, s)
        a <- dmatrix (numCols t, n)
        
        let b = t' <**> a
        return (TriBandedSM t' b)
        
    coarbitrary = undefined
    