{-# LANGUAGE Rank2Types #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.Tri
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Triangular views of matrices.
--

module Numeric.LinearAlgebra.Matrix.Tri (
    -- * Immutable interface
    
    -- ** Vector multiplication
    triMulVector,
    
    -- ** Matrix multiplication
    triMulMatrix,
    triMulMatrixWithScale,

    -- ** Vector solving
    triSolvVector,
    
    -- ** Matrix solving
    triSolvMatrix,
    triSolvMatrixWithScale,

    -- * Mutable interface
    triCreate,
    
    -- ** Vector multiplication
    triMulVectorM_,
    
    -- ** Matrix multiplication
    triMulMatrixM_,
    triMulMatrixWithScaleM_,

    -- ** Vector solving
    triSolvVectorM_,

    -- ** Matrix solving
    triSolvMatrixM_,
    triSolvMatrixWithScaleM_,

    ) where

import Control.Monad( when )
import Control.Monad.ST( ST, runST )
import Control.Monad.ST.Unsafe( unsafeIOToST )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Vector( Vector, STVector )
import qualified Numeric.LinearAlgebra.Vector as V
import Numeric.LinearAlgebra.Matrix.Base( Matrix )
import Numeric.LinearAlgebra.Matrix.STBase( STMatrix, RMatrix )
import qualified Numeric.LinearAlgebra.Matrix.STBase as M
import Numeric.LinearAlgebra.Types
import qualified Foreign.BLAS as BLAS


-- | A safe way to create and work with a mutable Tri Matrix before returning 
-- an immutable one for later perusal.
triCreate :: (Storable e)
           => (forall s. ST s (Tri (STMatrix s) e))
           -> Tri Matrix e
triCreate mt = runST $ do
    (Tri u d ma) <- mt
    a <- M.unsafeFreeze ma
    return $ Tri u d a

-- | @triMulVector trans a x@ returns @op(a) * x@, where @op(a)@ is
-- determined by @trans@.
triMulVector :: (BLAS2 e)
             => Trans
             -> Tri Matrix e
             -> Vector e
             -> Vector e
triMulVector trans a x =
    V.create $ do
        x' <- V.newCopy x
        triMulVectorM_ trans a x'
        return x'

-- | @triMulMatrix side a b@
-- returns @alpha * op(a) * b@ when @side@ is @LeftSide@ and
-- @alpha * b * op(a)@ when @side@ is @RightSide@.  Operation
-- @op(a)@ is determined by @trans@.
triMulMatrix :: (BLAS3 e)
              => Side
              -> Trans -> Tri Matrix e
              -> Matrix e
              -> Matrix e
triMulMatrix side trans a b = 
    M.create $ do
        b' <- M.newCopy b
        triMulMatrixM_ side trans a b'
        return b'

-- | @triMulMatrixWithScale alpha side trans a b@
-- returns @alpha * op(a) * b@ when @side@ is @LeftSide@ and
-- @alpha * b * op(a)@ when @side@ is @RightSide@.  Operation
-- @op(a)@ is determined by @trans@.
triMulMatrixWithScale :: (BLAS3 e)
                       => e
                       -> Side
                       -> Trans -> Tri Matrix e
                       -> Matrix e
                       -> Matrix e
triMulMatrixWithScale alpha side trans a b =
    M.create $ do
        b' <- M.newCopy b
        triMulMatrixWithScaleM_ alpha side trans a b'
        return b'

-- | @triMulVectorM_ a x@ sets @x := op(a) * x@, where @op(a)@ is determined
-- by @trans@.
triMulVectorM_ :: (RMatrix m, BLAS2 e)
               => Trans -> Tri m e
               -> STVector s e
               -> ST s ()
triMulVectorM_ trans (Tri uplo diag a) x = do
    (ma,na) <- M.getDim a
    nx <- V.getDim x
    let n = nx
    
    when (ma /= na) $ error $
        printf ("triMulVectorM_"
                ++ " _"
                ++ " (Tri _ _ <matrix with dim (%d,%d)>)"
                ++ " _"
                ++ ": matrix is not square")
               ma na
               
    when ((not . and) [ (ma,na) == (n,n)
                      , nx == n
                      ]) $ error $
        printf ("triMulVectorM_"
                ++ " _"
                ++ " (Tri _ _ <matrix with dim (%d,%d)>)"
                ++ " <vector with dim %d>"
                ++ ": dimension mismatch")
               ma na
               nx

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        V.unsafeWith x $ \px ->
            BLAS.trmv uplo trans diag n pa lda px 1


-- | @triMulMatrixM_ side trans a b@
-- sets @b := op(a) * b@ when @side@ is @LeftSide@ and
-- @b := b * op(a)@ when @side@ is @RightSide@.  Operation
-- @op(a)@ is determined by @trans@.
triMulMatrixM_ :: (RMatrix m, BLAS3 e)
               => Side 
               -> Trans -> Tri m e
               -> STMatrix s e
               -> ST s ()
triMulMatrixM_ = triMulMatrixWithScaleM_ 1

-- | @triMulMatrixWithScaleM_ alpha side trans a b@
-- sets @b := alpha * op(a) * b@ when @side@ is @LeftSide@ and
-- @b := alpha * b * op(a)@ when @side@ is @RightSide@.  Operation
-- @op(a)@ is determined by @trans@.
triMulMatrixWithScaleM_ :: (RMatrix m, BLAS3 e)
                         => e
                         -> Side
                         -> Trans -> Tri m e
                         -> STMatrix s e
                         -> ST s ()
triMulMatrixWithScaleM_ alpha side trans (Tri uplo diag a) b = do
    (ma,na) <- M.getDim a
    (mb,nb) <- M.getDim b
    let (m,n) = (mb,nb)
    
    when (ma /= na) $ error $
        printf ("triMulMatrixWithScaleM_"
                ++ " _"
                ++ " _"
                ++ " _"
                ++ " (Tri _ _ <matrix with dim (%d,%d)>)"
                ++ " _"
                ++ ": matrix is not square")
               ma na

    when ((not . and) [ case side of LeftSide  -> (ma,na) == (m,m)
                                     RightSide -> (ma,na) == (n,n)
                      , (mb, nb ) == (m,n)
                      ]) $ error $
        printf ("triMulMatrixWithScaleM_"
                ++ " _"
                ++ " %s"
                ++ " _"
                ++ " (Tri _ _ <matrix with dim (%d,%d)>)"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
               (show side)
               ma na
               mb nb

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        M.unsafeWith b $ \pb ldb ->
            BLAS.trmm side uplo trans diag m n alpha pa lda pb ldb


-- | @triSolvVector trans a x@ returns @op(a) \\ x@, where @op(a)@ is
-- determined by @trans@.
triSolvVector :: (BLAS2 e)
             => Trans
             -> Tri Matrix e
             -> Vector e
             -> Vector e
triSolvVector trans a x =
    V.create $ do
        x' <- V.newCopy x
        triSolvVectorM_ trans a x'
        return x'

-- | @triSolvMatrix side a b@
-- returns @alpha * op(a) \\ b@ when @side@ is @LeftSide@ and
-- @alpha * b * op(a)@ when @side@ is @RightSide@.  Operation
-- @op(a)@ is determined by @trans@.
triSolvMatrix :: (BLAS3 e)
              => Side
              -> Trans -> Tri Matrix e
              -> Matrix e
              -> Matrix e
triSolvMatrix side trans a b = 
    M.create $ do
        b' <- M.newCopy b
        triSolvMatrixM_ side trans a b'
        return b'

-- | @triSolvMatrixWithScale alpha side trans a b@
-- returns @alpha * op(a) \\ b@ when @side@ is @LeftSide@ and
-- @alpha * b * op(a)@ when @side@ is @RightSide@.  Operation
-- @op(a)@ is determined by @trans@.
triSolvMatrixWithScale :: (BLAS3 e)
                       => e
                       -> Side
                       -> Trans -> Tri Matrix e
                       -> Matrix e
                       -> Matrix e
triSolvMatrixWithScale alpha side trans a b =
    M.create $ do
        b' <- M.newCopy b
        triSolvMatrixWithScaleM_ alpha side trans a b'
        return b'

-- | @triSolvVectorM_ a x@ sets @x := op(a) \\ x@, where @op(a)@ is determined
-- by @trans@.
triSolvVectorM_ :: (RMatrix m, BLAS2 e)
               => Trans -> Tri m e
               -> STVector s e
               -> ST s ()
triSolvVectorM_ trans (Tri uplo diag a) x = do
    (ma,na) <- M.getDim a
    nx <- V.getDim x
    let n = nx
    
    when (ma /= na) $ error $
        printf ("triSolvVectorM_"
                ++ " _"
                ++ " (Tri _ _ <matrix with dim (%d,%d)>)"
                ++ " _"
                ++ ": matrix is not square")
               ma na
               
    when ((not . and) [ (ma,na) == (n,n)
                      , nx == n
                      ]) $ error $
        printf ("triSolvVectorM_"
                ++ " _"
                ++ " (Tri _ _ <matrix with dim (%d,%d)>)"
                ++ " <vector with dim %d>"
                ++ ": dimension mismatch")
               ma na
               nx

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        V.unsafeWith x $ \px ->
            BLAS.trsv uplo trans diag n pa lda px 1


-- | @triSolvMatrixM_ side trans a b@
-- sets @b := op(a) \\ b@ when @side@ is @LeftSide@ and
-- @b := b * op(a)@ when @side@ is @RightSide@.  Operation
-- @op(a)@ is determined by @trans@.
triSolvMatrixM_ :: (RMatrix m, BLAS3 e)
               => Side 
               -> Trans -> Tri m e
               -> STMatrix s e
               -> ST s ()
triSolvMatrixM_ = triSolvMatrixWithScaleM_ 1

-- | @triSolvMatrixWithScaleM_ alpha side trans a b@
-- sets @b := alpha * op(a) \\ b@ when @side@ is @LeftSide@ and
-- @b := alpha * b * op(a)@ when @side@ is @RightSide@.  Operation
-- @op(a)@ is determined by @trans@.
triSolvMatrixWithScaleM_ :: (RMatrix m, BLAS3 e)
                         => e
                         -> Side
                         -> Trans -> Tri m e
                         -> STMatrix s e
                         -> ST s ()
triSolvMatrixWithScaleM_ alpha side trans (Tri uplo diag a) b = do
    (ma,na) <- M.getDim a
    (mb,nb) <- M.getDim b
    let (m,n) = (mb,nb)
    
    when (ma /= na) $ error $
        printf ("triSolvMatrixWithScaleM_"
                ++ " _"
                ++ " _"
                ++ " _"
                ++ " (Tri _ _ <matrix with dim (%d,%d)>)"
                ++ " _"
                ++ ": matrix is not square")
               ma na

    when ((not . and) [ case side of LeftSide  -> (ma,na) == (m,m)
                                     RightSide -> (ma,na) == (n,n)
                      , (mb, nb ) == (m,n)
                      ]) $ error $
        printf ("triSolvMatrixWithScaleM_"
                ++ " _"
                ++ " %s"
                ++ " _"
                ++ " (Tri _ _ <matrix with dim (%d,%d)>)"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
               (show side)
               ma na
               mb nb

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        M.unsafeWith b $ \pb ldb ->
            BLAS.trsm side uplo trans diag m n alpha pa lda pb ldb
