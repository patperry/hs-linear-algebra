{-# LANGUAGE FlexibleInstances #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.Matrix.Tri.Banded
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.Matrix.Tri.Banded (
    TriBanded(..),
    TriBandedMV(..),
    TriBandedMM(..),
    TriBandedSV(..),
    TriBandedSM(..),
    ) where

import Control.Monad ( replicateM )

import Test.QuickCheck hiding ( vector )
import Test.QuickCheck.BLAS ( TestElem )
import qualified Test.QuickCheck as QC
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense ( Vector )
import Data.Matrix.Dense ( Matrix )
import Data.Matrix.Banded
import Data.Elem.BLAS ( Elem, BLAS1, BLAS3 )

import Data.Matrix.Tri ( Tri, triFromBase )
import Data.Matrix.Class( UpLoEnum(..), DiagEnum(..) )

listsFromBanded :: (Elem e) => Banded e -> ((Int,Int), (Int,Int),[[e]])
listsFromBanded a = ( (m,n)
            , (kl,ku)
            , map paddedDiag [(-kl)..ku]
            )
  where
    (m,n)   = shape a
    (kl,ku) = bandwidths a
    
    padBegin i   = replicate (max (-i) 0)    0
    padEnd   i   = replicate (max (m-n+i) 0) 0
    paddedDiag i = (  padBegin i
                   ++ elems (diagBanded a i)
                   ++ padEnd i 
                   )
                   
triBanded :: (TestElem e) => UpLoEnum -> DiagEnum -> Int -> Int -> Gen (Banded e)
triBanded Upper NonUnit n k = do
    a <- triBanded Upper Unit n k
    d <- Test.elems n
    let (_,_,(_:ds)) = listsFromBanded a
    return $ listsBanded (n,n) (0,k) (d:ds)

triBanded Lower NonUnit n k = do
    a <- triBanded Lower Unit n k
    d <- Test.elems n
    let (_,_,ds) = listsFromBanded a
        ds' = (init ds) ++ [d]
    return $ listsBanded (n,n) (k,0) ds'
    
triBanded _ Unit n 0 = do
    return $ listsBanded (n,n) (0,0) [replicate n 1]
    
triBanded Upper Unit n k = do
    a <- triBanded Upper Unit n (k-1)
    let (_,_,ds) = listsFromBanded a
    
    d <- Test.elems (n-k) >>= \xs -> return $ xs ++ replicate k 0
    
    return $ listsBanded (n,n) (0,k) $ ds ++ [d]
    
triBanded Lower Unit n k = do
    a <- triBanded Lower Unit n (k-1)
    let (_,_,ds) = listsFromBanded a

    d <- Test.elems (n-k) >>= \xs -> return $ replicate k 0 ++ xs
    
    return $ listsBanded (n,n) (k,0) $ [d] ++ ds
    
    
data TriBanded e = 
    TriBanded (Tri Banded e) (Banded e) deriving Show

instance (TestElem e) => Arbitrary (TriBanded e) where
    arbitrary = do
        u <- elements [ Upper, Lower  ]
        d <- elements [ Unit, NonUnit ]
        (m,n) <- Test.shape
        (_,k) <- Test.bandwidths (m,n)
        a <- triBanded u d n k

        l    <- if n == 0 then return 0 else choose (0,n-1)
        junk <- replicateM l $ Test.elems n
        diagJunk <- Test.elems n
        let (_,_,ds) = listsFromBanded a
            t = triFromBase u d $ case (u,d) of 
                    (Upper,NonUnit) -> 
                        listsBanded (n,n) (l,k) $ junk ++ ds
                    (Upper,Unit) ->
                        listsBanded (n,n) (l,k) $ junk ++ [diagJunk] ++ tail ds
                    (Lower,NonUnit) -> 
                        listsBanded (n,n) (k,l) $ ds ++ junk
                    (Lower,Unit) -> 
                        listsBanded (n,n) (k,l) $ init ds ++ [diagJunk] ++ junk

        (t',a') <- elements [ (t,a), (herm t, herm a)]
        return $ TriBanded t' a'

data TriBandedMV e = 
    TriBandedMV (Tri Banded e) (Banded e) (Vector e) deriving Show

instance (TestElem e) => Arbitrary (TriBandedMV e) where
    arbitrary = do
        (TriBanded t a) <- arbitrary
        x <- Test.vector (numCols t)
        return $ TriBandedMV t a x
        
data TriBandedMM e = 
    TriBandedMM (Tri Banded e) (Banded e) (Matrix e) deriving Show

instance (TestElem e, BLAS3 e) => Arbitrary (TriBandedMM e) where
    arbitrary = do
        (TriBanded t a) <- arbitrary
        (_,n) <- Test.shape
        b <- Test.matrix (numCols t, n)
        return $ TriBandedMM t a b
        
data TriBandedSV e = 
    TriBandedSV (Tri Banded e) (Vector e) deriving (Show)
    
instance (TestElem e, BLAS3 e) => Arbitrary (TriBandedSV e) where
    arbitrary = do
        (TriBanded t a) <- arbitrary
        if any (== 0) (elems $ diagBanded a 0)
            then arbitrary
            else do
                x <- Test.vector (numCols t)
                let y = t <*> x
                return (TriBandedSV t y)


data TriBandedSM e = 
    TriBandedSM (Tri Banded e) (Matrix e) 
    deriving (Show)
    
instance (TestElem e, BLAS3 e) => Arbitrary (TriBandedSM e) where
    arbitrary = do
        (TriBandedSV t _) <- arbitrary
        (_,n) <- Test.shape
        a <- Test.matrix (numCols t, n)
        
        let b = t <**> a
        return (TriBandedSM t b)
    