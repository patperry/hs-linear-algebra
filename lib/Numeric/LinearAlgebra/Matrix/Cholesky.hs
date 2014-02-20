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

import Numeric.LinearAlgebra.Matrix.Base( Matrix )
import Numeric.LinearAlgebra.Matrix.STBase( RMatrix, STMatrix )
import qualified Numeric.LinearAlgebra.Matrix.STBase as M

import Numeric.LinearAlgebra.Vector( Vector, STVector )
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
    cholFactorM (Herm uplo ma)
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
    x' <- V.newCopy x
    cholSolveVectorM_ a x'
    return x'

-- | @cholSolveMatrix a b@ returns @a \\ b@.
cholSolveMatrix :: (LAPACK e)
                => Chol Matrix e
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
            => Herm (STMatrix s) e
            -> ST s (Either Int (Chol (STMatrix s) e))
cholFactorM (Herm uplo a) = do
    (ma,na) <- M.getDim a
    let n = na

    when (not $ (ma,na) == (n,n)) $ error $
        printf ("cholFactorM"
                ++ " (Herm _ <matrix with dim (%d,%d)>): nonsquare matrix"
               ) ma na

    unsafeIOToST $
        M.unsafeWith a $ \pa lda -> do
            info <- LAPACK.potrf uplo n pa lda
            return $ if info > 0 then Left info
                                 else Right (Chol uplo a)


-- | @cholSolveVectorM_ a x@ sets @x := a \\ x@.
cholSolveVectorM_ :: (LAPACK e, RMatrix m)
                  => Chol m e
                  -> STVector s e
                  -> ST s ()
cholSolveVectorM_ a x =
    M.withFromColM x $ \x' ->
        cholSolveMatrixM_ a x'

-- | @cholSolveMatrixM_ a b@ sets @b := a \\ b@.
cholSolveMatrixM_ :: (LAPACK e, RMatrix m)
                  => Chol m e
                  -> STMatrix s e
                  -> ST s ()
cholSolveMatrixM_ (Chol uplo a) b = do
    (ma,na) <- M.getDim a
    (mb,nb) <- M.getDim b
    let (n,nrhs) = (mb,nb)
    
    when ((not . and) [ (ma,na) == (n,n)
                      , (mb,nb) == (n,nrhs)
                      ]) $ error $
        printf ("cholSolveMatrixM_"
                ++ " (Chol _ <matrix with dim (%d,%d)>)"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
               ma na mb nb

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        M.unsafeWith b $ \pb ldb ->
            LAPACK.potrs uplo n nrhs pa lda pb ldb
