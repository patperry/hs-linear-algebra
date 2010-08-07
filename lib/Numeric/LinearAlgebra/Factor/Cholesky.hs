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

module Numeric.LinearAlgebra.Factor.Cholesky (
    Chol(..),
    cholFactorHermMatrix,
    cholFactorToHermMatrix,
    cholMatrixSolveVector,
    cholMatrixSolveToVector,
    cholMatrixSolveMatrix,
    cholMatrixSolveToMatrix,
    ) where

import Control.Monad.ST( ST, runST )

-- import qualified Foreign.LAPACK as LAPACK

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Matrix
import Numeric.LinearAlgebra.Matrix.Herm
import Numeric.LinearAlgebra.Matrix.ST
import Numeric.LinearAlgebra.Vector
import Numeric.LinearAlgebra.Vector.ST

-- | A Cholesky decomposition view of a matrix.
data Chol m e = Chol Uplo (m e) deriving (Show)

cholFactorHermMatrix :: (LAPACK e)
                     => Herm Matrix e
                     -> Either Int (Chol Matrix e)
cholFactorHermMatrix (Herm uplo a) = runST $ do
    ma <- newCopyMatrix a
    cholFactorToHermMatrix (Herm uplo ma)
        >>= either (return . Left) (\(Chol uplo' ma') -> do
                a' <- unsafeFreezeMatrix ma'
                return $ Right (Chol uplo' a')
                )

cholMatrixSolveVector :: (LAPACK e)
                      => Chol Matrix e
                      -> Vector e
                      -> Vector e
cholMatrixSolveVector a x = runVector $ do
    y <- newVector_ (dimVector x)
    cholMatrixSolveToVector a x y
    return y

cholMatrixSolveMatrix :: (LAPACK e)
                      => Chol Matrix e
                      -> Matrix e
                      -> Matrix e
cholMatrixSolveMatrix a c = runMatrix $ do
    c' <- newMatrix_ (dimMatrix c)
    cholMatrixSolveToMatrix a c c'
    return c'

cholFactorToHermMatrix :: (LAPACK e)
                       => Herm (STMatrix s) e
                       -> ST s (Either Int (Chol (STMatrix s) e))
cholFactorToHermMatrix = undefined

cholMatrixSolveToVector :: (LAPACK e, RMatrix m, RVector v)
                        => Chol m e
                        -> v e
                        -> STVector s e
                        -> ST s ()
cholMatrixSolveToVector a x y =
    withMatrixViewColVector x $ \x' ->
        cholMatrixSolveToMatrix a x' (matrixViewColVector y)

cholMatrixSolveToMatrix :: (LAPACK e, RMatrix m1, RMatrix m2)
                        => Chol m1 e
                        -> m2 e
                        -> STMatrix s e
                        -> ST s ()
cholMatrixSolveToMatrix = undefined
