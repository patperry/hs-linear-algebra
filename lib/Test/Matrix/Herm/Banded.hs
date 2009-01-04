{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.Matrix.Herm.Banded
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.Matrix.Herm.Banded (
    HermBanded(..),
    HermBandedMV(..),
    HermBandedMM(..),
    ) where

import Control.Monad ( replicateM )

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Test.Vector.Dense ( vector )
import Test.Matrix ( matrixSized )
import Test.Matrix.Dense ( matrix )

import Data.Elem.BLAS ( Elem, BLAS2, BLAS3, toReal, fromReal, conj )
import BLAS.Types ( flipUpLo )

import Data.Vector.Dense ( Vector )
import Data.Matrix.Banded
import Data.Matrix.Dense ( Matrix )
import Data.Matrix.Herm



hermBanded :: (BLAS2 e, Arbitrary e) => Int -> Int -> Gen (Banded (n,n) e)
hermBanded n k 
    | n < 0 = 
        error $ "hermBanded: n must be non-negative"
    | n == 0 =
        return $ listsBanded (0,0) (0,0) []
    | k >= n =
        error $ "hermBanded: k must be less than n"
    | k < 0 = 
        error $ "hermBanded: k must be non-negative"
    | k == 0 = do
        d <- QC.vector n >>= return . realPart
        return $ listsBanded (n,n) (0,0) [d]
    | otherwise = do
        a <- hermBanded n (k-1)
        let (_,_,ds) = listsFromBanded a
        
        d <- QC.vector (n-k)
        let d'  = map conj d
            pad = replicate k 0
            ds' = [pad ++ d] ++ ds ++ [d' ++ pad]

        return $ listsBanded (n,n) (k,k) ds'
        
  where
    realPart :: Elem e => [e] -> [e]
    realPart = map (fromReal . toReal)

data HermBanded n e =
    HermBanded (Herm Banded (n,n) e)
               (Banded (n,n) e)
    deriving Show
    
instance (Arbitrary e, BLAS2 e) => Arbitrary (HermBanded n e) where
    arbitrary = matrixSized $ \s -> do
        n <- choose (0,s)
        k <- if n == 0 then return 0 else choose (0,n-1)
        l <- if n == 0 then return 0 else choose (0,n-1)
            
        a <- hermBanded n k
            
        junk <- replicateM l $ QC.vector n
        let (_,_,ds) = listsFromBanded a
            (u ,b ) = (Upper, listsBanded (n,n) (l,k) $ junk ++ (drop k ds))
            (u',b') = (Lower, listsBanded (n,n) (k,l) $ (take (k+1) ds) ++ junk)
        
        h <- elements [ hermFromBase u             b
                      , hermFromBase (flipUpLo u)  (herm b)
                      , hermFromBase u'            b'
                      , hermFromBase (flipUpLo u') (herm b')
                      ]
            
        return $ HermBanded h a

    coarbitrary = undefined

data HermBandedMV n e = 
    HermBandedMV (Herm Banded (n,n) e) 
                 (Banded (n,n) e) 
                 (Vector n e) 
    deriving Show

instance (Arbitrary e, BLAS2 e) => Arbitrary (HermBandedMV n e) where
    arbitrary = do
        (HermBanded h a) <- arbitrary
        x <- vector (numCols a)
        return $ HermBandedMV h a x

    coarbitrary = undefined
    
    
data HermBandedMM m n e = 
    HermBandedMM (Herm Banded (m,m) e) 
                 (Banded (m,m) e) 
                 (Matrix (m,n) e) 
    deriving Show
    
instance (Arbitrary e, BLAS3 e) => Arbitrary (HermBandedMM m n e) where
    arbitrary = matrixSized $ \s -> do
        (HermBanded a h) <- arbitrary
        n <- choose (0,s)
        b <- matrix (numCols h,n)

        return $ HermBandedMM a h b
            
    coarbitrary = undefined
