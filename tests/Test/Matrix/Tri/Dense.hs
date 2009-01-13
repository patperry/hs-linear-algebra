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

import Test.QuickCheck hiding ( Test.vector )
import Test.QuickCheck.BLAS ( TestElem )
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense hiding ( Test.vector )
import Data.Matrix.Dense
import Data.Elem.BLAS ( BLAS3 )

import Data.Matrix.Tri ( Tri, triFromBase )
import Data.Matrix.Class( UpLoEnum(..), DiagEnum(..) )

import Unsafe.Coerce

triMatrix :: (TestElem e) => UpLoEnum -> DiagEnum -> (Int,Int) -> Gen (Matrix (m,n) e)
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
data TriMatrix m n e = 
    TriMatrix (Tri Matrix (m,n) e) 
              (Matrix (m,n) e) 
    deriving Show

instance (TestElem e) => Arbitrary (TriMatrix m n e) where
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
            
        (t',a') <- elements [ (t,a), unsafeCoerce (herm t, herm a) ]
            
        return $ TriMatrix t' $ submatrix a' (0,0) (shape t')
            
    coarbitrary = undefined


-- | A triangular Test.matrix, and equivalent dense Test.matrix, and a Test.vector in
-- their domain.        
data TriMatrixMV m n e = 
    TriMatrixMV (Tri Matrix (m,n) e) 
                (Matrix (m,n) e) 
                (Vector n e) 
    deriving Show

instance (TestElem e) => Arbitrary (TriMatrixMV m n e) where
    arbitrary = do
        (TriMatrix t a) <- arbitrary
        x <- Test.vector (numCols a)
        return $ TriMatrixMV t a x

    coarbitrary = undefined

-- | A triangular Test.matrix, and equivalent dense Test.matrix, and a Test.matrix in
-- their domain.        
data TriMatrixMM m k n e = 
    TriMatrixMM (Tri Matrix (m,k) e) 
                (Matrix (m,k) e) 
                (Matrix (k,n) e) 
    deriving Show

instance (TestElem e) => Arbitrary (TriMatrixMM m k n e) where
    arbitrary =  do
        (TriMatrix t a) <- arbitrary
        n <- liftM fst Test.shape
        b <- Test.matrix (numCols a, n)
        return $ TriMatrixMM t a b
            
    coarbitrary = undefined

-- | A triangular Test.matrix and a Test.vector in its range
data TriMatrixSV m n e = 
    TriMatrixSV (Tri Matrix (m,n) e) 
                (Vector m e) 
    deriving Show
    
instance (BLAS3 e, TestElem e) => Arbitrary (TriMatrixSV m n e) where
    arbitrary = do
        (TriMatrix t a) <- arbitrary
        if any (== 0) (elems $ diag a 0)
            then arbitrary
            else do
                x <- Test.vector (numCols a)
                let y  = a <*> x
                return (TriMatrixSV t y)
        
    coarbitrary = undefined

-- | A triangular Test.matrix and a Test.matrix in its range
data TriMatrixSM m k n e = 
    TriMatrixSM (Tri Matrix (m,k) e) 
                (Matrix (m,n) e) 
    deriving Show
    
instance (BLAS3 e, TestElem e) => Arbitrary (TriMatrixSM m k n e) where
    arbitrary = do
        (TriMatrix t a) <- arbitrary
        if any (== 0) (elems $ diag a 0)
            then arbitrary
            else do
                n <- liftM fst Test.shape
                b <- Test.matrix (numCols a, n)
                let c  = a <**> b
                return (TriMatrixSM t c)
        
    coarbitrary = undefined

normF :: (BLAS3 e) => Matrix (n,p) e -> Double
normF = sum . map norm2 . cols
