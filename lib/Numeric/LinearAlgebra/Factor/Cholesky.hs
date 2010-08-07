-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Factor.Cholesky
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Cholesky decompositions of symmetric (Hermitian) positive-definite
-- matrices.
--

-- | A Cholesky decomposition view of a matrix.
data Chol m e = Chol Uplo (m e) deriving (Show)

class LAPACK e where

cholFactorHermMatrix :: (LAPACK e)
                     => Herm Matrix e
                     -> ST s (Either Int (Chol Matrix e))


cholFactorToHermMatrix :: (LAPACK e, RMatrix m)
                       => Herm m e
                       => Herm (STMatrix s) e
                       => ST s (Either Int (Chol (STMatrix s) e))

cholMatrixSolveVector :: (LAPACK e)
                      => Chol Matrix e
                      -> Vector e
                      -> Vector e

cholMatrixSolveMatrix :: (LAPACK e)
                      => Chol Matrix e
                      -> Matrix e
                      -> Matrix e

cholMatrixSolveToVector :: (LAPACK e, RMatrix m, RVector v)
                        -> Chol m e
                        -> v e
                        -> STVector s e
                        -> ST s ()

cholMatrixSolveToMatrix :: (LAPACK e, RMatrix m1 RMatrix m2)
                        -> Chol m1 e
                        -> m2 e
                        -> STMatrix s e
                        -> ST s ()
