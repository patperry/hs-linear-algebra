-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.Eigen
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Eigenvalue decompositions.
--
module Numeric.LinearAlgebra.Matrix.Eigen (
    hermEigen,
    hermEigenvalues,
    hermEigenM_,
    hermEigenvaluesM_,
    ) where

import Control.Monad( when )
import Control.Monad.ST( ST, runST )
import Control.Monad.ST.Unsafe( unsafeIOToST )
import Foreign( allocaArray, nullPtr )
import Text.Printf( printf )

import Foreign.LAPACK( LAPACK, EigJob(..), EigRange(..) )
import qualified Foreign.LAPACK as LAPACK

import Numeric.LinearAlgebra.Types( Herm(..) )
import Numeric.LinearAlgebra.Matrix.Base( Matrix )
import Numeric.LinearAlgebra.Matrix.STBase( STMatrix )
import qualified Numeric.LinearAlgebra.Matrix.STBase as M
import Numeric.LinearAlgebra.Vector( Vector, STVector )
import qualified Numeric.LinearAlgebra.Vector as V

-- | Compute the eigenvalues and eigenvectors of a Hermitian matrix.
-- Return the eigenvalues are in ascending order in the result vector; 
-- store the corresponding eigenvectors are in the columns of the result
-- matrix.
hermEigen :: (LAPACK e) => Herm Matrix e -> (Vector Double, Matrix e)
hermEigen (Herm uplo a) = runST $ do
    (_,n) <- M.getDim a
    ma <- M.newCopy a
    mw <- V.new_ n
    mz <- M.new_ (n,n)
    hermEigenM_ (Herm uplo ma) mw mz
    w <- V.unsafeFreeze mw
    z <- M.unsafeFreeze mz
    return (w,z)


-- | Return the eigenvalues of a Hermitian matrix in ascending order.
hermEigenvalues :: (LAPACK e) => Herm Matrix e -> Vector Double
hermEigenvalues (Herm uplo a) = V.create $ do
    (_,n) <- M.getDim a
    ma <- M.newCopy a
    w <- V.new_ n
    hermEigenvaluesM_ (Herm uplo ma) w
    return w


-- | Compute and copy the eigenvalues and eigenvectors of a mutable
-- Hermitian matrix.  This destroys the original Hermitian matrix.
hermEigenM_ :: (LAPACK e)
                  => Herm (STMatrix s) e
                  -> STVector s Double
                  -> STMatrix s e
                  -> ST s ()
hermEigenM_ (Herm uplo a) w z = do
    (ma,na) <- M.getDim a
    (mz,nz) <- M.getDim z
    nw <- V.getDim w
    let n = na
    
    when (ma /= na) $ error $
        printf ("hermEigenM_"
                ++ " (Herm _ <matrix with dim (%d,%d)>): nonsquare matrix"
               ) ma na
    when ((not . and) [ (ma,na) == (n,n)
                      , nw == n
                      , (mz,nz) == (n,n)
                      ]) $ error $
        printf ("hermEigenM_"
                ++ " (Herm _ <matrix with dim (%d,%d)>)"
                ++ " <vector with dim %d>:"
                ++ " <matrix with dim (%d,%d)>"
                ++ " dimension mismatch")
               ma na nw mz nz

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        V.unsafeWith w $ \pw ->
        M.unsafeWith z $ \pz ldz ->
        allocaArray (2*n) $ \psuppz -> do
            _m <- LAPACK.heevr jobz range uplo n pa lda abstol pw pz ldz
                               psuppz
            return ()
  where
    jobz = EigVectors
    range = AllEigs
    abstol = 0
      

-- | Compute and copy the eigenvalues of a mutable Hermitian matrix.  This
-- destroys the original Hermitian matrix.
hermEigenvaluesM_ :: (LAPACK e)
                        => Herm (STMatrix s) e
                        -> STVector s Double
                        -> ST s ()
hermEigenvaluesM_ (Herm uplo a) w = do
    (ma,na) <- M.getDim a
    nw <- V.getDim w
    let n = na
    
    when (ma /= na) $  error $
        printf ("hermEigenvaluesM_"
                ++ " (Herm _ <matrix with dim (%d,%d)>): nonsquare matrix"
               ) ma na

    when ((not . and) [ (ma,na) == (n,n)
                      , nw == n
                      ]) $ error $
        printf ("hermEigenvaluesM_"
                ++ " (Herm _ <matrix with dim (%d,%d)>)"
                ++ " <vector with dim %d>: dimension mismatch")
               ma na nw
               
    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        V.unsafeWith w $ \pw -> do
            _m <- LAPACK.heevr jobz range uplo n pa lda abstol pw pz ldz
                               psuppz
            return ()
    
  where
    jobz = NoEigVectors
    range = AllEigs
    abstol = 0
    pz = nullPtr
    ldz = 1
    psuppz = nullPtr
