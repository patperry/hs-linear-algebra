{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.Matrix.Tri.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.Matrix.Tri.Dense (
    TriMatrix(..),
    TriMatrixMV(..),
    TriMatrixMM(..),
    TriMatrixSV(..),
    TriMatrixSM(..),
    ) where

import Data.Ix ( range )
import Control.Monad( liftM )

import Test.QuickCheck hiding ( vector )
import Test.QuickCheck.BLAS ( TestElem )
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense hiding ( vector )
import Data.Matrix.Dense
import Data.Elem.BLAS ( BLAS3 )

import Data.Matrix.Tri ( Tri, triFromBase )
import Data.Matrix.Class( UpLoEnum(..), DiagEnum(..) )

triMatrix :: (TestElem e) => UpLoEnum -> DiagEnum -> (Int,Int) -> Gen (Matrix e)
triMatrix u d (m,n) =
    let ijs = filter (isTriIndex u d) $ range ((0,0), (m-1,n-1))
    in do
        es <- Test.elems (m*n)
        let a  = matrix (m,n) $ zip ijs es
            a' = case d of
                    NonUnit -> a
                    Unit    -> a // [ ((i,i),1) | i <- [0..(mn-1)] ]
        return $ a'
  where
    mn = min m n

isTriIndex :: UpLoEnum -> DiagEnum -> (Int,Int) -> Bool
isTriIndex Upper NonUnit (i,j) = i <= j
isTriIndex Upper Unit    (i,j) = i <  j
isTriIndex Lower NonUnit (i,j) = i >= j
isTriIndex Lower Unit    (i,j) = i >  j

-- | A triangular Test.matrix and an equivalent dense Test.matrix
data TriMatrix e = 
    TriMatrix (Tri Matrix e) 
              (Matrix e) 
    deriving Show

instance (TestElem e) => Arbitrary (TriMatrix e) where
    arbitrary = do
        u <- elements [ Upper, Lower  ]
        d <- elements [ Unit, NonUnit ]
        (m,n) <- Test.shape
        a <- triMatrix u d (m,n)
        
        junk <- Test.elems (m*n)
        let ijs = [ (i,j) | i <- [0..(m-1)]
                          , j <- [0..(n-1)]
                          , (not . (isTriIndex u d)) (i,j) ]
            t   = triFromBase u d $ a // zip ijs junk
            
        (t',a') <- elements [ (t,a), (herm t, herm a) ]
            
        return $ TriMatrix t' $ submatrix a' (0,0) (shape t')


-- | A triangular Test.matrix, and equivalent dense Test.matrix, and a Test.vector in
-- their domain.        
data TriMatrixMV e = 
    TriMatrixMV (Tri Matrix e) 
                (Matrix e) 
                (Vector e) 
    deriving Show

instance (TestElem e) => Arbitrary (TriMatrixMV e) where
    arbitrary = do
        (TriMatrix t a) <- arbitrary
        x <- Test.vector (numCols a)
        return $ TriMatrixMV t a x

-- | A triangular Test.matrix, and equivalent dense Test.matrix, and a Test.matrix in
-- their domain.        
data TriMatrixMM e = 
    TriMatrixMM (Tri Matrix e) 
                (Matrix e) 
                (Matrix e) 
    deriving Show

instance (TestElem e) => Arbitrary (TriMatrixMM e) where
    arbitrary =  do
        (TriMatrix t a) <- arbitrary
        n <- liftM fst Test.shape
        b <- Test.matrix (numCols a, n)
        return $ TriMatrixMM t a b

-- | A triangular Test.matrix and a Test.vector in its range
data TriMatrixSV e = 
    TriMatrixSV (Tri Matrix e) 
                (Vector e) 
    deriving Show
    
instance (BLAS3 e, TestElem e) => Arbitrary (TriMatrixSV e) where
    arbitrary = do
        (TriMatrix t a) <- arbitrary
        if any (== 0) (elems $ diag a 0)
            then arbitrary
            else do
                x <- Test.vector (numCols a)
                let y  = a <*> x
                return (TriMatrixSV t y)

-- | A triangular Test.matrix and a Test.matrix in its range
data TriMatrixSM e = 
    TriMatrixSM (Tri Matrix e) 
                (Matrix e) 
    deriving Show
    
instance (BLAS3 e, TestElem e) => Arbitrary (TriMatrixSM e) where
    arbitrary = do
        (TriMatrix t a) <- arbitrary
        if any (== 0) (elems $ diag a 0)
            then arbitrary
            else do
                n <- liftM fst Test.shape
                b <- Test.matrix (numCols a, n)
                let c  = a <**> b
                return (TriMatrixSM t c)

normF :: (BLAS3 e) => Matrix e -> Double
normF = sum . map norm2 . cols
