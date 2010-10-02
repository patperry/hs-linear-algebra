-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.Cholesky
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Cholesky factorizations of symmetric (Hermitian) positive-definite
-- matrices.
--

module Numeric.LinearAlgebra.Matrix.Cholesky (
    -- * Immutable interface
    cholFactor,
    cholSolveVector,
    cholSolveMatrix,
    
    -- * Mutable interface
    cholFactorTo,
    cholSolveToVector,
    cholSolveToMatrix,
    
    ) where

import Control.Monad.ST( ST, runST, unsafeIOToST )
import Text.Printf( printf )

import qualified Foreign.LAPACK as LAPACK

import Numeric.LinearAlgebra.Types

import Numeric.LinearAlgebra.Matrix.Base( Matrix )
import Numeric.LinearAlgebra.Matrix.STBase( RMatrix, STMatrix )
import qualified Numeric.LinearAlgebra.Matrix.Base as M
import qualified Numeric.LinearAlgebra.Matrix.STBase as M

import Numeric.LinearAlgebra.Vector( Vector, RVector, STVector )
import qualified Numeric.LinearAlgebra.Vector as V


-- | @cholFactor a@ tries to compute the Cholesky
-- factorization of @a@.  If @a@ is positive-definite then the routine
-- returns @Right@ with the factorization.  If the leading minor of order @i@
-- is not positive-definite then the routine returns @Left i@.
cholFactor :: (LAPACK e)
                 => Herm Matrix e
                 -> Either Int (Chol Matrix e)
cholFactor (Herm uplo a) = runST $ do
    ma <- M.newCopy a
    cholFactorTo (Herm uplo ma)
        >>= either (return . Left) (\(Chol uplo' ma') -> do
                a' <- M.unsafeFreeze ma'
                return $ Right (Chol uplo' a')
                )

-- | @cholSolveVector a x@ returns @a \\ x@.
cholSolveVector :: (LAPACK e)
                      => Chol Matrix e
                      -> Vector e
                      -> Vector e
cholSolveVector a x = V.create $ do
    y <- V.new_ (V.dim x)
    cholSolveToVector a x y
    return y

-- | @cholSolveMatrix a b@ returns @a \\ b@.
cholSolveMatrix :: (LAPACK e)
                      => Chol Matrix e
                      -> Matrix e
                      -> Matrix e
cholSolveMatrix a c = M.create $ do
    c' <- M.new_ (M.dim c)
    cholSolveToMatrix a c c'
    return c'

-- | @cholFactorTo a@ tries to compute the Cholesky
-- factorization of @a@ in place.  If @a@ is positive-definite then the
-- routine returns @Right@ with the factorization, stored in the same
-- memory as @a@.  If the leading minor of order @i@ is not
-- positive-definite then the routine returns @Left i@.
-- In either case, the original storage of @a@ is destroyed.
cholFactorTo :: (LAPACK e)
                   => Herm (STMatrix s) e
                   -> ST s (Either Int (Chol (STMatrix s) e))
cholFactorTo (Herm uplo a)
    | not $ (ma,na) == (n,n) = error $
        printf ("cholFactorTo"
                ++ " (Herm _ <matrix with dim (%d,%d)>): nonsquare matrix"
               ) ma na
    | otherwise =
        unsafeIOToST $
            M.unsafeWith a $ \pa lda -> do
                info <- LAPACK.potrf uplo n pa lda
                return $ if info > 0 then Left info
                                     else Right (Chol uplo a)
  where
    (ma,na) = M.dim a
    n = na

-- | @cholSolveToVector a x x'@ sets
-- @x' := a \\ x@.  Arguments @x@ and @x'@ can be the same.
cholSolveToVector :: (LAPACK e, RMatrix m, RVector v)
                        => Chol m e
                        -> v e
                        -> STVector s e
                        -> ST s ()
cholSolveToVector a x y =
    M.withViewFromCol x $ \x' ->
    M.withViewFromSTCol y $ \y' ->
        cholSolveToMatrix a x' y'

-- | @cholSolveToMatrix a b b'@ sets
-- @b' := a \\ b@.  Arguments @b@ and @b'@ can be the same.
cholSolveToMatrix :: (LAPACK e, RMatrix m1, RMatrix m2)
                        => Chol m1 e
                        -> m2 e
                        -> STMatrix s e
                        -> ST s ()
cholSolveToMatrix (Chol uplo a) b b'
    | (not . and) [ (ma,na) == (n,n)
                  , (mb,nb) == (n,nrhs)
                  , (mb',nb') == (n,nrhs)
                  ] = error $
        printf ("cholSolveToMatrix"
                ++ " (Chol _ <matrix with dim (%d,%d)>)"
                ++ " <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
               ma na mb nb mb' nb'
    | otherwise = do
        M.unsafeCopyTo b b'
        unsafeIOToST $
            M.unsafeWith a $ \pa lda ->
            M.unsafeWith b' $ \pb ldb ->
                LAPACK.potrs uplo n nrhs pa lda pb ldb
  where
    (ma,na) = M.dim a
    (mb,nb) = M.dim b
    (mb',nb') = M.dim b
    (n,nrhs) = (mb',nb')
