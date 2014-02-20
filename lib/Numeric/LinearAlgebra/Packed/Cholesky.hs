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
    cholFactorM,
    cholSolveVectorM_,
    cholSolveMatrixM_,
    ) where

import Control.Monad( when )
import Control.Monad.ST( ST, runST )
import Control.Monad.ST.Unsafe( unsafeIOToST )
import Text.Printf( printf )

import qualified Foreign.LAPACK as LAPACK

import Numeric.LinearAlgebra.Types

import Numeric.LinearAlgebra.Packed.Base( Packed, RPacked, STPacked )
import qualified Numeric.LinearAlgebra.Packed.Base as P

import Numeric.LinearAlgebra.Matrix( Matrix, STMatrix )
import qualified Numeric.LinearAlgebra.Matrix as M

import Numeric.LinearAlgebra.Vector( Vector, STVector )
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
    cholFactorM (Herm uplo ma)
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
    x' <- V.newCopy x
    cholSolveVectorM_ a x'
    return x'

-- | @cholSolveMatrix a b@ returns @a \\ b@.
cholSolveMatrix :: (LAPACK e)
                      => Chol Packed e
                      -> Matrix e
                      -> Matrix e
cholSolveMatrix a c = M.create $ do
    c' <- M.newCopy c
    cholSolveMatrixM_ a c'
    return c'

-- | @cholFactorM a@ tries to compute the Cholesky
-- factorization of @a@ in place.  If @a@ is positive-definite then the
-- routine returns @Right@ with the factorization, stored in the same
-- memory as @a@.  If the leading minor of order @i@ is not
-- positive-definite then the routine returns @Left i@.
-- In either case, the original storage of @a@ is destroyed.
cholFactorM :: (LAPACK e)
             => Herm (STPacked s) e
             -> ST s (Either Int (Chol (STPacked s) e))
cholFactorM (Herm uplo a) = do
    n <- P.getDim a
    unsafeIOToST $
        P.unsafeWith a $ \pa -> do
            info <- LAPACK.pptrf uplo n pa
            return $ if info > 0 then Left info
                                 else Right (Chol uplo a)

-- | @cholSolveVectorM_ a x@ sets @x := a \\ x@.
cholSolveVectorM_ :: (LAPACK e, RPacked p)
                  => Chol p e
                  -> STVector s e
                  -> ST s ()
cholSolveVectorM_ a x =
    M.withFromColM x $ \x' ->
        cholSolveMatrixM_ a x'

-- | @cholSolveMatrixM_ a b@ sets @b := a \\ b@.
cholSolveMatrixM_ :: (LAPACK e, RPacked p)
                  => Chol p e
                  -> STMatrix s e
                  -> ST s ()
cholSolveMatrixM_ (Chol uplo a) b = do
    na <- P.getDim a
    (mb,nb) <- M.getDim b
    let (n,nrhs) = (mb,nb)

    when ((not . and) [ na == n
                      , (mb,nb) == (n,nrhs)
                      ]) $ error $
        printf ("cholSolveMatrixM_"
                ++ " (Chol _ <packed matrix with dim %d>)"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
               na mb nb

    unsafeIOToST $
        P.unsafeWith a $ \pa ->
        M.unsafeWith b $ \pb ldb ->
            LAPACK.pptrs uplo n nrhs pa pb ldb
