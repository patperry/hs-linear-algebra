-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Factor.Eigen
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Eigenvalue decompositions.
--
module Numeric.LinearAlgebra.Factor.Eigen (
    eigenHermMatrix,
    eigenvaluesHermMatrix,
    eigenToHermMatrix,
    eigenvaluesToHermMatrix,
    ) where

import Control.Monad.ST( ST, runST, unsafeIOToST )
import Foreign( allocaArray, nullPtr )
import Text.Printf( printf )

import Foreign.LAPACK( LAPACK, EigJob(..), EigRange(..) )
import qualified Foreign.LAPACK as LAPACK

import Numeric.LinearAlgebra.Matrix
import Numeric.LinearAlgebra.Matrix.Herm
import Numeric.LinearAlgebra.Matrix.ST
import Numeric.LinearAlgebra.Vector( Vector )
import qualified Numeric.LinearAlgebra.Vector as V
import Numeric.LinearAlgebra.Vector.ST( STVector )
import qualified Numeric.LinearAlgebra.Vector.ST as V

-- | Compute the eigenvalues and eigenvectors of a Hermitian matrix.
-- Return the eigenvalues are in ascending order in the result vector; 
-- store the corresponding eigenvectors are in the columns of the result
-- matrix.
eigenHermMatrix :: (LAPACK e) => Herm Matrix e -> (Vector Double, Matrix e)
eigenHermMatrix (Herm uplo a) = runST $ do
    ma <- newCopyMatrix a
    mw <- V.new_ n
    mz <- newMatrix_ (n,n)
    eigenToHermMatrix (Herm uplo ma) mw mz
    w <- V.unsafeFreeze mw
    z <- unsafeFreezeMatrix mz
    return (w,z)
  where
    (_,n) = dimMatrix a

-- | Return the eigenvalues of a Hermitian matrix in ascending order.
eigenvaluesHermMatrix :: (LAPACK e) => Herm Matrix e -> Vector Double
eigenvaluesHermMatrix (Herm uplo a) = V.create $ do
    ma <- newCopyMatrix a
    w <- V.new_ n
    eigenvaluesToHermMatrix (Herm uplo ma) w
    return w
  where
    (_,n) = dimMatrix a

-- | Compute and copy the eigenvalues and eigenvectors of a mutable
-- Hermitian matrix.  This destroys the original Hermitian matrix.
eigenToHermMatrix :: (LAPACK e)
                  => Herm (STMatrix s) e
                  -> STVector s Double
                  -> STMatrix s e
                  -> ST s ()
eigenToHermMatrix (Herm uplo a) w z
    | ma /= na =  error $
        printf ("eigenToHermMatrix"
                ++ " (Herm _ <matrix with dim (%d,%d)>): nonsquare matrix"
               ) ma na
    | (not . and) [ (ma,na) == (n,n)
                  , nw == n
                  , (mz,nz) == (n,n)
                  ] = error $
        printf ("eigenToHermMatrix"
                ++ " (Herm _ <matrix with dim (%d,%d)>)"
                ++ " <vector with dim %d>:"
                ++ " <matrix with dim (%d,%d)>"
                ++ " dimension mismatch")
               ma na nw mz nz
    | otherwise =
        unsafeIOToST $
            unsafeWithMatrix a $ \pa lda ->
            V.unsafeWith w $ \pw ->
            unsafeWithMatrix z $ \pz ldz ->
            allocaArray (2*n) $ \psuppz -> do
                _m <- LAPACK.heevr jobz range uplo n pa lda abstol pw pz ldz
                                   psuppz
                return ()
    
  where
    (ma,na) = dimMatrix a
    (mz,nz) = dimMatrix z
    nw = V.dim w
    n = na
    jobz = EigVectors
    range = AllEigs
    abstol = 0

-- | Compute and copy the eigenvalues of a mutable Hermitian matrix.  This
-- destroys the original Hermitian matrix.
eigenvaluesToHermMatrix :: (LAPACK e)
                        => Herm (STMatrix s) e
                        -> STVector s Double
                        -> ST s ()
eigenvaluesToHermMatrix (Herm uplo a) w
    | ma /= na =  error $
        printf ("eigenvaluesToHermMatrix"
                ++ " (Herm _ <matrix with dim (%d,%d)>): nonsquare matrix"
               ) ma na
    | (not . and) [ (ma,na) == (n,n)
                  , nw == n
                  ] = error $
        printf ("eigenvaluesToHermMatrix"
                ++ " (Herm _ <matrix with dim (%d,%d)>)"
                ++ " <vector with dim %d>: dimension mismatch")
               ma na nw
    | otherwise =
        unsafeIOToST $
            unsafeWithMatrix a $ \pa lda ->
            V.unsafeWith w $ \pw -> do
                _m <- LAPACK.heevr jobz range uplo n pa lda abstol pw pz ldz
                                   psuppz
                return ()
    
  where
    (ma,na) = dimMatrix a
    nw = V.dim w
    n = na
    jobz = NoEigVectors
    range = AllEigs
    abstol = 0
    pz = nullPtr
    ldz = 1
    psuppz = nullPtr
