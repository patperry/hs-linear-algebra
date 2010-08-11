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
    
    -- * Immutable interface
    cholFactorMatrix,
    cholMatrixSolveVector,
    cholMatrixSolveMatrix,

    -- * Mutable interface
    cholFactorToMatrix,
    cholMatrixSolveToVector,
    cholMatrixSolveToMatrix,
    ) where

import Control.Monad.ST( ST, runST, unsafeIOToST )
import Text.Printf( printf )

import qualified Foreign.LAPACK as LAPACK

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Matrix
import Numeric.LinearAlgebra.Matrix.Herm
import Numeric.LinearAlgebra.Matrix.ST
import Numeric.LinearAlgebra.Matrix.STBase
import Numeric.LinearAlgebra.Vector
import Numeric.LinearAlgebra.Vector.ST

-- | A Cholesky decomposition view of a matrix.
data Chol m e = Chol Uplo (m e) deriving (Show)

-- | @cholFactorMatrix a@ tries to compute the Cholesky
-- factorization of @a@.  If @a@ is positive-definite then the routine
-- returns @Right@ with the factorization.  If the leading minor of order @i@
-- is not positive-definite then the routine returns @Left i@.
cholFactorMatrix :: (LAPACK e)
                 => Herm Matrix e
                 -> Either Int (Chol Matrix e)
cholFactorMatrix (Herm uplo a) = runST $ do
    ma <- newCopyMatrix a
    cholFactorToMatrix (Herm uplo ma)
        >>= either (return . Left) (\(Chol uplo' ma') -> do
                a' <- unsafeFreezeMatrix ma'
                return $ Right (Chol uplo' a')
                )

-- | @cholMatrixSolveVector a x@ returns @a \\ x@.
cholMatrixSolveVector :: (LAPACK e)
                      => Chol Matrix e
                      -> Vector e
                      -> Vector e
cholMatrixSolveVector a x = runVector $ do
    y <- newVector_ (dimVector x)
    cholMatrixSolveToVector a x y
    return y

-- | @cholMatrixSolveMatrix a b@ returns @a \\ b@.
cholMatrixSolveMatrix :: (LAPACK e)
                      => Chol Matrix e
                      -> Matrix e
                      -> Matrix e
cholMatrixSolveMatrix a c = runMatrix $ do
    c' <- newMatrix_ (dimMatrix c)
    cholMatrixSolveToMatrix a c c'
    return c'

-- | @cholFactorToMatrix a@ tries to compute the Cholesky
-- factorization of @a@ in place.  If @a@ is positive-definite then the
-- routine returns @Right@ with the factorization, stored in the same
-- memory as @a@.  If the leading minor of order @i@ is not
-- positive-definite then the routine returns @Left i@.
-- In either case, the original storage of @a@ is destroyed.
cholFactorToMatrix :: (LAPACK e)
                   => Herm (STMatrix s) e
                   -> ST s (Either Int (Chol (STMatrix s) e))
cholFactorToMatrix (Herm uplo a)
    | not $ (ma,na) == (n,n) = error $
        printf ("cholFactorToMatrix"
                ++ " (Herm _ <matrix with dim (%d,%d)>): nonsquare matrix"
               ) ma na
    | otherwise =
        unsafeIOToST $
            unsafeWithMatrix a $ \pa lda -> do
                info <- LAPACK.potrf uplo n pa lda
                return $ if info > 0 then Left info
                                     else Right (Chol uplo a)
  where
    (ma,na) = dimMatrix a
    n = na

-- | @cholMatrixSolveToVector a x x'@ sets
-- @x' := a \\ x@.  Arguments @x@ and @x'@ can be the same.
cholMatrixSolveToVector :: (LAPACK e, RMatrix m, RVector v)
                        => Chol m e
                        -> v e
                        -> STVector s e
                        -> ST s ()
cholMatrixSolveToVector a x y =
    withMatrixViewColVector x $ \x' ->
        cholMatrixSolveToMatrix a x' (matrixViewColVector y)

-- | @cholMatrixSolveToMatrix a b b'@ sets
-- @b' := a \\ b@.  Arguments @b@ and @b'@ can be the same.
cholMatrixSolveToMatrix :: (LAPACK e, RMatrix m1, RMatrix m2)
                        => Chol m1 e
                        -> m2 e
                        -> STMatrix s e
                        -> ST s ()
cholMatrixSolveToMatrix (Chol uplo a) b b'
    | (not . and) [ (ma,na) == (n,n)
                  , (mb,nb) == (n,nrhs)
                  , (mb',nb') == (n,nrhs)
                  ] = error $
        printf ("cholMatrixSolveToMatrix"
                ++ " (Chol _ <matrix with dim (%d,%d)>)"
                ++ " <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
               ma na mb nb mb' nb'
    | otherwise = do
        unsafeCopyToMatrix b b'
        unsafeIOToST $
            unsafeWithMatrix a $ \pa lda ->
            unsafeWithMatrix b' $ \pb ldb ->
                LAPACK.potrs uplo n nrhs pa lda pb ldb
  where
    (ma,na) = dimMatrix a
    (mb,nb) = dimMatrix b
    (mb',nb') = dimMatrix b
    (n,nrhs) = (mb',nb')
