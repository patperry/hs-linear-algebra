-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Packed.Cholesky
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Cholesky factorizations of symmetric (Hermitian) positive-definite
-- matrices.
--

module Numeric.LinearAlgebra.Packed.Cholesky (
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

import Numeric.LinearAlgebra.Packed.Base( Packed, RPacked, STPacked )
import qualified Numeric.LinearAlgebra.Packed.Base as P

import Numeric.LinearAlgebra.Matrix( Matrix, RMatrix, STMatrix )
import qualified Numeric.LinearAlgebra.Matrix as M

import Numeric.LinearAlgebra.Vector( Vector, RVector, STVector )
import qualified Numeric.LinearAlgebra.Vector as V

-- | @cholFactor a@ tries to compute the Cholesky
-- factorization of @a@.  If @a@ is positive-definite then the routine
-- returns @Right@ with the factorization.  If the leading minor of order @i@
-- is not positive-definite then the routine returns @Left i@.
cholFactor :: (LAPACK e)
                 => Herm Packed e
                 -> Either Int (Chol Packed e)
cholFactor (Herm uplo a) = runST $ do
    ma <- P.newCopy a
    cholFactorTo (Herm uplo ma)
        >>= either (return . Left) (\(Chol uplo' ma') -> do
                a' <- P.unsafeFreeze ma'
                return $ Right (Chol uplo' a')
                )

-- | @cholSolveVector a x@ returns @a \\ x@.
cholSolveVector :: (LAPACK e)
                      => Chol Packed e
                      -> Vector e
                      -> Vector e
cholSolveVector a x = V.create $ do
    y <- V.new_ (V.dim x)
    cholSolveToVector a x y
    return y

-- | @cholSolveMatrix a b@ returns @a \\ b@.
cholSolveMatrix :: (LAPACK e)
                      => Chol Packed e
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
                   => Herm (STPacked s) e
                   -> ST s (Either Int (Chol (STPacked s) e))
cholFactorTo (Herm uplo a) =
    unsafeIOToST $
        P.unsafeWith a $ \pa -> do
            info <- LAPACK.pptrf uplo n pa
            return $ if info > 0 then Left info
                                 else Right (Chol uplo a)
  where
    n = P.dim a

-- | @cholSolveToVector a x x'@ sets
-- @x' := a \\ x@.  Arguments @x@ and @x'@ can be the same.
cholSolveToVector :: (LAPACK e, RPacked p, RVector v)
                        => Chol p e
                        -> v e
                        -> STVector s e
                        -> ST s ()
cholSolveToVector a x y =
    M.withViewFromCol x $ \x' ->
    M.withViewFromSTCol y $ \y' ->
        cholSolveToMatrix a x' y'

-- | @cholSolveToMatrix a b b'@ sets
-- @b' := a \\ b@.  Arguments @b@ and @b'@ can be the same.
cholSolveToMatrix :: (LAPACK e, RPacked p, RMatrix m)
                        => Chol p e
                        -> m e
                        -> STMatrix s e
                        -> ST s ()
cholSolveToMatrix (Chol uplo a) b b'
    | (not . and) [ na == n
                  , (mb,nb) == (n,nrhs)
                  , (mb',nb') == (n,nrhs)
                  ] = error $
        printf ("cholSolveToMatrix"
                ++ " (Chol _ <packed matrix with dim %d>)"
                ++ " <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
               na mb nb mb' nb'
    | otherwise = do
        M.unsafeCopyTo b' b
        unsafeIOToST $
            P.unsafeWith a $ \pa ->
            M.unsafeWith b' $ \pb ldb ->
                LAPACK.pptrs uplo n nrhs pa pb ldb
  where
    na = P.dim a
    (mb,nb) = M.dim b
    (mb',nb') = M.dim b
    (n,nrhs) = (mb',nb')
