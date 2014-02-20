{-# LANGUAGE Rank2Types #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Packed.Tri
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Triangular views of packed matrices.
--

module Numeric.LinearAlgebra.Packed.Tri (
    -- * Immutable interface
    
    -- ** Vector multiplication
    triMulVector,
    
    -- ** Vector solving
    triSolvVector,
    
    -- * Mutable interface
    triCreate,
    
    -- ** Vector multiplication
    triMulVectorM_,
    
    -- ** Vector solving
    triSolvVectorM_,

    ) where

import Control.Monad( when )
import Control.Monad.ST( ST, runST )
import Control.Monad.ST.Unsafe( unsafeIOToST )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Vector( Vector, STVector )
import qualified Numeric.LinearAlgebra.Vector as V
import Numeric.LinearAlgebra.Packed.Base( Packed, RPacked, STPacked )
import qualified Numeric.LinearAlgebra.Packed.Base as P
import Numeric.LinearAlgebra.Types
import qualified Foreign.BLAS as BLAS


-- | A safe way to create and work with a mutable Tri Packed before returning 
-- an immutable one for later perusal.
triCreate :: (Storable e)
          => (forall s. ST s (Tri (STPacked s) e))
          -> Tri Packed e
triCreate mt = runST $ do
    (Tri u d ma) <- mt
    a <- P.unsafeFreeze ma
    return $ Tri u d a

-- | @triMulVector trans a x@ returns @op(a) * x@, where @op(a)@ is
-- determined by @trans@.
triMulVector :: (BLAS2 e)
             => Trans
             -> Tri Packed e
             -> Vector e
             -> Vector e
triMulVector trans a x =
    V.create $ do
        x' <- V.newCopy x
        triMulVectorM_ trans a x'
        return x'

-- | @triMulVectorM_ a x@ sets @x := op(a) * x@, where @op(a)@ is determined
-- by @trans@.
triMulVectorM_ :: (RPacked p, BLAS2 e)
               => Trans -> Tri p e
               -> STVector s e
               -> ST s ()
triMulVectorM_ trans (Tri uplo diag a) x = do
    na <- P.getDim a
    nx <- V.getDim x
    let n = nx
    
    when (nx /= n) $ error $
        printf ("triMulVectorM_"
                ++ " _"
                ++ " (Tri _ _ <packed matrix with dim %d>)"
                ++ " <vector with dim %d>"
                ++ ": dimension mismatch")
               na
               nx

    unsafeIOToST $
        P.unsafeWith a $ \pa ->
        V.unsafeWith x $ \px ->
            BLAS.tpmv uplo trans diag n pa px 1


-- | @triSolvVector trans a x@ returns @op(a) \\ x@, where @op(a)@ is
-- determined by @trans@.
triSolvVector :: (BLAS2 e)
             => Trans
             -> Tri Packed e
             -> Vector e
             -> Vector e
triSolvVector trans a x =
    V.create $ do
        x' <- V.newCopy x
        triSolvVectorM_ trans a x'
        return x'

-- | @triSolvVectorM_ a x@ sets @x := op(a) \\ x@, where @op(a)@ is determined
-- by @trans@.
triSolvVectorM_ :: (RPacked p, BLAS2 e)
               => Trans -> Tri p e
               -> STVector s e
               -> ST s ()
triSolvVectorM_ trans (Tri uplo diag a) x = do
    na <- P.getDim a
    nx <- V.getDim x
    let n = nx
    
    when (nx /= n) $ error $
        printf ("triMulVectorM_"
                ++ " _"
                ++ " (Tri _ _ <packed matrix with dim %d>)"
                ++ " <vector with dim %d>"
                ++ ": dimension mismatch")
               na
               nx

    unsafeIOToST $
        P.unsafeWith a $ \pa ->
        V.unsafeWith x $ \px ->
            BLAS.tpsv uplo trans diag n pa px 1
