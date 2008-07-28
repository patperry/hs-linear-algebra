-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Matrix.Herm.Banded
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Matrix.Herm.Banded 
    where

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Test.QuickCheck.Vector.Dense ( dvector )
import Test.QuickCheck.Matrix.Dense ( dmatrix )

import BLAS.Elem ( Elem, BLAS2, toReal, fromReal, conj )

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
        let (_,_,ds) = toLists a
        
        d <- QC.vector (n-k)
        let d'  = map conj d
            pad = replicate k 0
            ds' = [pad ++ d] ++ ds ++ [d' ++ pad]

        return $ listsBanded (n,n) (k,k) ds'
        
  where
    realPart :: Elem e => [e] -> [e]
    realPart = map (fromReal . toReal)


data HermBandedMV n e = 
    HermBandedMV UpLo (Banded (n,n) e) (Banded (n,n) e) (Vector n e) deriving (Show)

instance (Arbitrary e, BLAS2 e) => Arbitrary (HermBandedMV n e) where
    arbitrary = sized $ \s ->
        let s' = ceiling (sqrt $ fromInteger $ toInteger s :: Double)
        in do
            u <- elements [ Upper, Lower ]
            n <- choose (0,s')
            k <- if n == 0 then return 0 else choose (0,n-1)
            
            h <- hermBanded n k
            
            let (_,_,ds) = toLists h
                a = case u of
                        Upper -> listsBanded (n,n) (0,k) (drop k ds)
                        Lower -> listsBanded (n,n) (k,0) (take (k+1) ds)
                        
            x <- dvector n
            
            return $ HermBandedMV u a h x

    coarbitrary = undefined

    
data HermBandedMM m n e = 
    HermBandedMM UpLo (Banded (m,m) e) (Banded (m,m) e) (Matrix (m,n) e) deriving (Show)
    
instance (Arbitrary e, BLAS2 e) => Arbitrary (HermBandedMM m n e) where
    arbitrary = sized $ \s ->
        let s' = ceiling (sqrt $ fromInteger $ toInteger s :: Double)
        in do
            (HermBandedMV u a h _) <- arbitrary
            n <- choose (0,s')
            
            let m = numCols h
            b <- dmatrix (m,n)

            return $ HermBandedMM u a h b
            
    coarbitrary = undefined
