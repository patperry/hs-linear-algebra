-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Factor.Cholesky
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Cholesky factorizations of symmetric (Hermitian) positive-definite
-- matrices.
--

module Numeric.LinearAlgebra.Factor.Cholesky (
    Chol(..),
    
    -- * Immutable interface
    -- ** Matrix
    cholFactorMatrix,
    cholMatrixSolveVector,
    cholMatrixSolveMatrix,
    
    -- ** Packed
    cholFactorPacked,
    cholPackedSolveVector,
    cholPackedSolveMatrix,

    -- * Mutable interface
    -- ** Matrix
    cholFactorToMatrix,
    cholMatrixSolveToVector,
    cholMatrixSolveToMatrix,
    
    -- ** Packed
    cholFactorToPacked,
    cholPackedSolveToVector,
    cholPackedSolveToMatrix,
    ) where

import Control.Monad.ST( ST, runST, unsafeIOToST )
import Text.Printf( printf )

import qualified Foreign.LAPACK as LAPACK

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Matrix.Herm

import Numeric.LinearAlgebra.Packed( Packed, RPacked, STPacked )
import qualified Numeric.LinearAlgebra.Packed as P

import Numeric.LinearAlgebra.Matrix.Base( Matrix )
import qualified Numeric.LinearAlgebra.Matrix.Base as M

import Numeric.LinearAlgebra.Matrix.STBase( STMatrix, RMatrix )
import qualified Numeric.LinearAlgebra.Matrix.STBase as M

import Numeric.LinearAlgebra.Vector( Vector, RVector, STVector )
import qualified Numeric.LinearAlgebra.Vector as V


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
    ma <- M.newCopy a
    cholFactorToMatrix (Herm uplo ma)
        >>= either (return . Left) (\(Chol uplo' ma') -> do
                a' <- M.unsafeFreeze ma'
                return $ Right (Chol uplo' a')
                )

-- | @cholMatrixSolveVector a x@ returns @a \\ x@.
cholMatrixSolveVector :: (LAPACK e)
                      => Chol Matrix e
                      -> Vector e
                      -> Vector e
cholMatrixSolveVector a x = V.create $ do
    y <- V.new_ (V.dim x)
    cholMatrixSolveToVector a x y
    return y

-- | @cholMatrixSolveMatrix a b@ returns @a \\ b@.
cholMatrixSolveMatrix :: (LAPACK e)
                      => Chol Matrix e
                      -> Matrix e
                      -> Matrix e
cholMatrixSolveMatrix a c = M.create $ do
    c' <- M.new_ (M.dim c)
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
            M.unsafeWith a $ \pa lda -> do
                info <- LAPACK.potrf uplo n pa lda
                return $ if info > 0 then Left info
                                     else Right (Chol uplo a)
  where
    (ma,na) = M.dim a
    n = na

-- | @cholMatrixSolveToVector a x x'@ sets
-- @x' := a \\ x@.  Arguments @x@ and @x'@ can be the same.
cholMatrixSolveToVector :: (LAPACK e, RMatrix m, RVector v)
                        => Chol m e
                        -> v e
                        -> STVector s e
                        -> ST s ()
cholMatrixSolveToVector a x y =
    M.withViewFromCol x $ \x' ->
        cholMatrixSolveToMatrix a x' (M.fromSTCol y)

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

-- | @cholFactorPacked a@ tries to compute the Cholesky
-- factorization of @a@.  If @a@ is positive-definite then the routine
-- returns @Right@ with the factorization.  If the leading minor of order @i@
-- is not positive-definite then the routine returns @Left i@.
cholFactorPacked :: (LAPACK e)
                 => Herm Packed e
                 -> Either Int (Chol Packed e)
cholFactorPacked (Herm uplo a) = runST $ do
    ma <- P.newCopy a
    cholFactorToPacked (Herm uplo ma)
        >>= either (return . Left) (\(Chol uplo' ma') -> do
                a' <- P.unsafeFreeze ma'
                return $ Right (Chol uplo' a')
                )

-- | @cholPackedSolveVector a x@ returns @a \\ x@.
cholPackedSolveVector :: (LAPACK e)
                      => Chol Packed e
                      -> Vector e
                      -> Vector e
cholPackedSolveVector a x = V.create $ do
    y <- V.new_ (V.dim x)
    cholPackedSolveToVector a x y
    return y

-- | @cholPackedSolveMatrix a b@ returns @a \\ b@.
cholPackedSolveMatrix :: (LAPACK e)
                      => Chol Packed e
                      -> Matrix e
                      -> Matrix e
cholPackedSolveMatrix a c = M.create $ do
    c' <- M.new_ (M.dim c)
    cholPackedSolveToMatrix a c c'
    return c'

-- | @cholFactorToPacked a@ tries to compute the Cholesky
-- factorization of @a@ in place.  If @a@ is positive-definite then the
-- routine returns @Right@ with the factorization, stored in the same
-- memory as @a@.  If the leading minor of order @i@ is not
-- positive-definite then the routine returns @Left i@.
-- In either case, the original storage of @a@ is destroyed.
cholFactorToPacked :: (LAPACK e)
                   => Herm (STPacked s) e
                   -> ST s (Either Int (Chol (STPacked s) e))
cholFactorToPacked (Herm uplo a) =
    unsafeIOToST $
        P.unsafeWith a $ \pa -> do
            info <- LAPACK.pptrf uplo n pa
            return $ if info > 0 then Left info
                                 else Right (Chol uplo a)
  where
    n = P.dim a

-- | @cholPackedSolveToVector a x x'@ sets
-- @x' := a \\ x@.  Arguments @x@ and @x'@ can be the same.
cholPackedSolveToVector :: (LAPACK e, RPacked p, RVector v)
                        => Chol p e
                        -> v e
                        -> STVector s e
                        -> ST s ()
cholPackedSolveToVector a x y =
    M.withViewFromCol x $ \x' ->
        cholPackedSolveToMatrix a x' (M.fromSTCol y)

-- | @cholPackedSolveToMatrix a b b'@ sets
-- @b' := a \\ b@.  Arguments @b@ and @b'@ can be the same.
cholPackedSolveToMatrix :: (LAPACK e, RPacked p, RMatrix m)
                        => Chol p e
                        -> m e
                        -> STMatrix s e
                        -> ST s ()
cholPackedSolveToMatrix (Chol uplo a) b b'
    | (not . and) [ na == n
                  , (mb,nb) == (n,nrhs)
                  , (mb',nb') == (n,nrhs)
                  ] = error $
        printf ("cholPackedSolveToMatrix"
                ++ " (Chol _ <packed matrix with dim %d>)"
                ++ " <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
               na mb nb mb' nb'
    | otherwise = do
        M.unsafeCopyTo b b'
        unsafeIOToST $
            P.unsafeWith a $ \pa ->
            M.unsafeWith b' $ \pb ldb ->
                LAPACK.pptrs uplo n nrhs pa pb ldb
  where
    na = P.dim a
    (mb,nb) = M.dim b
    (mb',nb') = M.dim b
    (n,nrhs) = (mb',nb')
