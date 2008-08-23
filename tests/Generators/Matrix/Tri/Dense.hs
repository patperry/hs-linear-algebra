-----------------------------------------------------------------------------
-- |
-- Module     : Generators.Matrix.Tri.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Generators.Matrix.Tri.Dense (
    TriMatrix(..),
    TriMatrixMV(..),
    TriMatrixMM(..),
    TriMatrixSV(..),
    TriMatrixSM(..),
    ) where

import Data.Ix ( range )

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Generators.Vector.Dense ( vector )
import qualified Generators.Matrix.Dense as Test
import Generators.Matrix ( matrixSized )

import Data.Vector.Dense hiding ( vector )
import Data.Matrix.Dense
import BLAS.Elem ( BLAS1, BLAS3 )

import Data.Matrix.Tri ( Tri, UpLo(..), Diag(..), fromBase )

import Unsafe.Coerce

triMatrix :: (BLAS1 e, Arbitrary e) => UpLo -> Diag -> (Int,Int) -> Gen (Matrix (m,n) e)
triMatrix u d (m,n) =
    let ijs = filter (isTriIndex u d) $ range ((0,0), (m-1,n-1))
    in do
        es <- QC.vector (m*n)
        let a = matrix (m,n) $ zip ijs es
            a' = case d of
                    NonUnit -> a
                    Unit    -> a // [ ((i,i),1) | i <- [0..(mn-1)] ]
        return $ a'
  where
    mn = min m n

isTriIndex :: UpLo -> Diag -> (Int,Int) -> Bool
isTriIndex Upper NonUnit (i,j) = i <= j
isTriIndex Upper Unit    (i,j) = i <  j
isTriIndex Lower NonUnit (i,j) = i >= j
isTriIndex Lower Unit    (i,j) = i >  j

-- | A triangular matrix and an equivalent dense matrix
data TriMatrix m n e = 
    TriMatrix (Tri Matrix (m,n) e) 
              (Matrix (m,n) e) 
    deriving Show

instance (Arbitrary e, BLAS1 e) => Arbitrary (TriMatrix m n e) where
    arbitrary = matrixSized $ \k -> do
        u <- elements [ Upper, Lower  ]
        d <- elements [ Unit, NonUnit ]
        m <- choose (0,k)
        n <- choose (0,k)
        a <- triMatrix u d (m,n)
        
        junk <- QC.vector (m*n)
        let ijs = [ (i,j) | i <- [0..(m-1)]
                          , j <- [0..(n-1)]
                          , (not . (isTriIndex u d)) (i,j) ]
            t   = fromBase u d $ a // zip ijs junk
            
        (t',a') <- elements [ (t,a), unsafeCoerce (herm t, herm a) ]
            
        return $ TriMatrix t' $ submatrix a' (0,0) (shape t')
            
    coarbitrary = undefined


-- | A triangular matrix, and equivalent dense matrix, and a vector in
-- their domain.        
data TriMatrixMV m n e = 
    TriMatrixMV (Tri Matrix (m,n) e) 
                (Matrix (m,n) e) 
                (Vector n e) 
    deriving Show

instance (Arbitrary e, BLAS1 e) => Arbitrary (TriMatrixMV m n e) where
    arbitrary = do
        (TriMatrix t a) <- arbitrary
        x <- vector (numCols a)
        return $ TriMatrixMV t a x

    coarbitrary = undefined

-- | A triangular matrix, and equivalent dense matrix, and a matrix in
-- their domain.        
data TriMatrixMM m k n e = 
    TriMatrixMM (Tri Matrix (m,k) e) 
                (Matrix (m,k) e) 
                (Matrix (k,n) e) 
    deriving Show

instance (Arbitrary e, BLAS1 e) => Arbitrary (TriMatrixMM m k n e) where
    arbitrary = matrixSized $ \s -> do
        (TriMatrix t a) <- arbitrary
        n <- choose (0,s)
        b <- Test.matrix (numCols a, n)
        return $ TriMatrixMM t a b
            
    coarbitrary = undefined

-- | A triangular matrix and a vector in its range
data TriMatrixSV m n e = 
    TriMatrixSV (Tri Matrix (m,n) e) 
                (Vector m e) 
    deriving Show
    
instance (Arbitrary e, BLAS3 e) => Arbitrary (TriMatrixSV m n e) where
    arbitrary = do
        (TriMatrix t a) <- arbitrary
        x <- vector (numCols a)
        let y = a <*> x
        return (TriMatrixSV t y)
        
    coarbitrary = undefined

-- | A triangular matrix and a matrix in its range
data TriMatrixSM m k n e = 
    TriMatrixSM (Tri Matrix (m,k) e) 
                (Matrix (m,n) e) 
    deriving Show
    
instance (Arbitrary e, BLAS3 e) => Arbitrary (TriMatrixSM m k n e) where
    arbitrary = matrixSized $ \s -> do
        (TriMatrix t a) <- arbitrary
        n <- choose (0, s)
        b <- Test.matrix (numCols a, n)
        let c = a <**> b
        return (TriMatrixSM t c)
        
    coarbitrary = undefined
    