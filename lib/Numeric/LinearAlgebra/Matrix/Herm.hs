{-# LANGUAGE Rank2Types #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.Herm
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Hermitian views of matrices.
--

module Numeric.LinearAlgebra.Matrix.Herm (
    -- * Hermitian views of matrices
    Herm(..),
    withHerm,
    
    -- * Immutable interface
    
    -- ** Matrix-Vector multiplication
    hermMulVector,
    hermMulVectorWithScale,
    hermMulAddVectorWithScales,
    
    -- ** Matrix-Matrix  multiplication
    hermMulMatrix,
    hermMulMatrixWithScale,
    hermMulAddMatrixWithScales,

    -- ** Updates
    hermRank1Update,
    hermRank2Update,
    hermRankKUpdate,   
    hermRank2KUpdate,    

    
    -- * Mutable interface
    hermCreate,
    
    -- ** Matrix-Vector multiplication
    hermMulToVector,
    hermMulToVectorWithScale,
    hermMulAddToVectorWithScales,
    
    -- ** Matrix-Matrix multiplication
    hermMulToMatrix,
    hermMulToMatrixWithScale,
    hermMulAddToMatrixWithScales,

    -- ** Updates
    hermRank1UpdateTo,
    hermRank2UpdateTo,
    hermRankKUpdateTo,  
    hermRank2KUpdateTo,
    
    ) where

import Control.Monad.ST( ST, runST, unsafeIOToST )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Vector( Vector, RVector, STVector )
import qualified Numeric.LinearAlgebra.Vector as V
import Numeric.LinearAlgebra.Matrix.Base( Matrix )
import qualified Numeric.LinearAlgebra.Matrix.Base as M
import Numeric.LinearAlgebra.Matrix.STBase( STMatrix, RMatrix )
import qualified Numeric.LinearAlgebra.Matrix.STBase as M
import Numeric.LinearAlgebra.Types
import qualified Foreign.BLAS as BLAS

-- | A hermitian view of an underlying matrix.  The view can either be
-- of the upper or lower triangular part of the matrix.  The type arguments
-- are as follows:
--
--     * @m@: the underlyting matrix type.
--
--     * @e@: the element type of the matrix.
--
data Herm m e = Herm Uplo (m e) deriving (Show)

-- | Apply a function to the unerlying 'Uplo' and matrix.
withHerm :: Herm m e -> (Uplo -> m e -> a) -> a
withHerm (Herm u m) f = f u m

-- | A safe way to V.create and work with a mutable Herm Matrix before returning 
-- an immutable one for later perusal.
hermCreate :: (Storable e)
              => (forall s. ST s (Herm (STMatrix s) e))
              -> Herm Matrix e
hermCreate mh = runST $ do
    (Herm u ma) <- mh
    a <- M.unsafeFreeze ma
    return $ Herm u a

-- | @hermRank1Update alpha x a@ returns
-- @alpha * x * x^H + a@.
hermRank1Update :: (BLAS2 e)
                      => Double -> Vector e -> Herm Matrix e -> Herm Matrix e
hermRank1Update alpha x (Herm uplo a) = runST $ do
    ma' <- M.newCopy a
    hermRank1UpdateTo alpha x (Herm uplo ma')
    a' <- M.unsafeFreeze ma'
    return $ Herm uplo a'

-- | @hermRank2Update alpha x y a@ returns
-- @alpha * x * y^H + conj(alpha) * y * x^H + a@.
hermRank2Update :: (BLAS2 e)
                      => e -> Vector e -> Vector e -> Herm Matrix e
                      -> Herm Matrix e
hermRank2Update alpha x y (Herm uplo a) = runST $ do
    ma' <- M.newCopy a
    hermRank2UpdateTo alpha x y (Herm uplo ma')
    a' <- M.unsafeFreeze ma'
    return $ Herm uplo a'

-- | @hermRankKUpdate alpha trans a beta c@ returns
-- @c := alpha * a * a^H + beta * c@ when @trans@ is @NoTrans@ and
-- @c := alpha * a^H * a + beta * c@ when @trans@ is @ConjTrans@.  The
-- function signals an error when @trans@ is @Trans@.
hermRankKUpdate :: (BLAS3 e)
                      => e -> Trans -> Matrix e -> e -> Herm Matrix e
                      -> Herm Matrix e
hermRankKUpdate alpha trans a beta (Herm uplo c) = runST $ do
    mc' <- M.newCopy c
    hermRankKUpdateTo alpha trans a beta (Herm uplo mc')
    c' <- M.unsafeFreeze mc'
    return $ Herm uplo c'

-- | @hermRank2KUpdate alpha trans a b beta c@ returns
-- @c := alpha * a * b^H + conj(alpha) * b * a^H + beta * c@ when @trans@ is
-- @NoTrans@ and @c := alpha * b^H * a + conj(alpha) * a^H * b + beta * c@
-- when @trans@ is @ConjTrans@.  The function signals an error when @trans@
-- is @Trans@.
hermRank2KUpdate :: (BLAS3 e)
                       => e -> Trans -> Matrix e -> Matrix e -> e -> Herm Matrix e
                       -> Herm Matrix e
hermRank2KUpdate alpha trans a b beta (Herm uplo c) = runST $ do
    mc' <- M.newCopy c
    hermRank2KUpdateTo alpha trans a b beta (Herm uplo mc')
    c' <- M.unsafeFreeze mc'
    return $ Herm uplo c'

-- | @hermRank1UpdateTo alpha x a@ sets
-- @a := alpha * x * x^H + a@.
hermRank1UpdateTo :: (RVector v, BLAS2 e)
                        => Double -> v e -> Herm (STMatrix s) e -> ST s ()
hermRank1UpdateTo alpha x (Herm uplo a)
    | (not . and) [ nx == n, (ma,na) == (n,n) ] = error $
        printf ("hermRank1UpdateTo _ <vector with dim %d>"
                 ++ " (Herm _ <matrix with dim (%d,%d)>):"
                 ++ " invalid dimensions") nx ma na
    | otherwise =
        unsafeIOToST $
        V.unsafeWith x $ \px ->
        M.unsafeWith a $ \pa lda ->
            BLAS.her uplo n alpha px 1 pa lda
  where
    nx = V.dim x
    (ma,na) = M.dim a
    n = nx

-- | @hermRank2UpdateTo alpha x y a@ sets
-- @a := alpha * x * y^H + conj(alpha) * y * x^H + a@.
hermRank2UpdateTo :: (RVector v1, RVector v2, BLAS2 e)
                        => e -> v1 e -> v2 e -> Herm (STMatrix s) e -> ST s ()
hermRank2UpdateTo alpha x y (Herm uplo a)
    | (not . and) [ nx == n, ny == n, (ma,na) == (n,n) ] = error $
        printf ("hermRank2UpdateTo _ <vector with dim %d>"
                 ++ " <vector with dim %d>"
                 ++ " (Herm _ <matrix with dim (%d,%d)>):"
                 ++ " invalid dimensions") nx ny ma na
    | otherwise =
        unsafeIOToST $
        V.unsafeWith x $ \px ->
        V.unsafeWith y $ \py ->
        M.unsafeWith a $ \pa lda ->
            BLAS.her2 uplo n alpha px 1 py 1 pa lda
  where
    nx = V.dim x
    ny = V.dim y
    (ma,na) = M.dim a
    n = nx

-- | @hermRankKUpdateTo alpha trans a beta c@ sets
-- @c := alpha * a * a^H + beta * c@ when @trans@ is @NoTrans@ and
-- @c := alpha * a^H * a + beta * c@ when @trans@ is @ConjTrans@.  The
-- function signals an error when @trans@ is @Trans@.
hermRankKUpdateTo :: (RMatrix m, BLAS3 e)
                        => e -> Trans -> m e -> e -> Herm (STMatrix s) e
                        -> ST s ()
hermRankKUpdateTo alpha trans a beta (Herm uplo c)
    | trans == Trans = error $
        printf ("hermRankKUpdateTo _ %s:"
                 ++ " trans argument must be NoTrans or ConjTrans")
               (show trans)
    | (not . and) [ (mc,nc) == (n,n)
                  , case trans of NoTrans -> (ma,na) == (n,k)
                                  _       -> (ma,na) == (k,n)
                  ] = error $
            printf ("hermRankKUpdateTo _ %s <matrix with dim (%d,%d)> _"
                    ++ " (Herm _ <matrix with dim (%d,%d)>):"
                    ++ " invalid dimensions") (show trans) ma na mc nc
    | otherwise =
        unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        M.unsafeWith c $ \pc ldc ->
            BLAS.herk uplo trans n k alpha pa lda beta pc ldc
  where
    (ma,na) = M.dim a
    (mc,nc) = M.dim c
    (n,k) = if trans == NoTrans then (ma,na) else (na,ma)


-- | @hermRank2KUpdateTo alpha trans a b beta c@ sets
-- @c := alpha * a * b^H + conj(alpha) * b * a^H + beta * c@ when @trans@ is
-- @NoTrans@ and @c := alpha * b^H * a + conj(alpha) * a^H * b + beta * c@
-- when @trans@ is @ConjTrans@.  The function signals an error when @trans@
-- is @Trans@.
hermRank2KUpdateTo :: (RMatrix m1, RMatrix m2, BLAS3 e)
                         => e -> Trans -> m1 e -> m2 e -> e -> Herm (STMatrix s) e
                         -> ST s ()
hermRank2KUpdateTo alpha trans a b beta (Herm uplo c)
    | trans == Trans = error $
        printf ("hermRank2KUpdateTo _ %s:"
                 ++ " trans argument must be NoTrans or ConjTrans")
               (show trans)
    | (not . and) [ (mc,nc) == (n,n)
                  , (mb,nb) == (ma,na)
                  , case trans of NoTrans -> (ma,na) == (n,k)
                                  _       -> (ma,na) == (k,n)
                  ] = error $
            printf ("hermRank2KUpdateTo _ %s <matrix with dim (%d,%d)>"
                    ++ " <matrix with dim (%d,%d)> _"
                    ++ " (Herm _ <matrix with dim (%d,%d)>):"
                    ++ " invalid dimensions") (show trans) ma na mb nb mc nc
    | otherwise =
        unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        M.unsafeWith b $ \pb ldb ->
        M.unsafeWith c $ \pc ldc ->
            BLAS.her2k uplo trans n k alpha pa lda pb ldb beta pc ldc
  where
    (ma,na) = M.dim a
    (mb,nb) = M.dim b
    (mc,nc) = M.dim c
    (n,k) = if trans == NoTrans then (ma,na) else (na,ma)


-- | @hermMulVector a x@ returns @a * x@.
hermMulVector :: (BLAS2 e)
                    => Herm Matrix e
                    -> Vector e
                    -> Vector e
hermMulVector a x =
    V.create $ do
        y <- V.new_ (V.dim x)
        hermMulToVector a x y
        return y

-- | @hermMulVectorWithScale alpha a x@ retunrs @alpha * a * x@.
hermMulVectorWithScale :: (BLAS2 e)
                             => e
                             -> Herm Matrix e
                             -> Vector e
                             -> Vector e
hermMulVectorWithScale alpha a x =
    V.create $ do
        y <- V.new_ (V.dim x)
        hermMulToVectorWithScale alpha a x y
        return y
                       
-- | @hermMulAddVectorWithScales alpha a x y@
-- returns @alpha * a * x + beta * y@.
hermMulAddVectorWithScales :: (BLAS2 e)
                                 => e
                                 -> Herm Matrix e
                                 -> Vector e
                                 -> e
                                 -> Vector e
                                 -> Vector e
hermMulAddVectorWithScales alpha a x beta y =
    V.create $ do
        y' <- V.newCopy y
        hermMulAddToVectorWithScales alpha a x beta y'
        return y'

-- | @hermMulMatrix side a b@
-- returns @alpha * a * b@ when @side@ is @LeftSide@ and
-- @alpha * b * a@ when @side@ is @RightSide@.
hermMulMatrix :: (BLAS3 e)
                    => Side -> Herm Matrix e
                    -> Matrix e
                    -> Matrix e
hermMulMatrix side a b = 
    M.create $ do
        c <- M.new_ (M.dim b)
        hermMulToMatrix side a b c
        return c

-- | @hermMulMatrixWithScale alpha side a b@
-- returns @alpha * a * b@ when @side@ is @LeftSide@ and
-- @alpha * b * a@ when @side@ is @RightSide@.
hermMulMatrixWithScale :: (BLAS3 e)
                             => e
                             -> Side -> Herm Matrix e
                             -> Matrix e
                             -> Matrix e
hermMulMatrixWithScale alpha side a b =
    M.create $ do
        c <- M.new_ (M.dim b)
        hermMulToMatrixWithScale alpha side a b c
        return c

-- | @hermMulAddMatrixWithScales alpha side a b beta c@
-- returns @alpha * a * b + beta * c@ when @side@ is @LeftSide@ and
-- @alpha * b * a + beta * c@ when @side@ is @RightSide@.
hermMulAddMatrixWithScales :: (BLAS3 e)
                                 => e
                                 -> Side -> Herm Matrix e
                                 -> Matrix e
                                 -> e
                                 -> Matrix e
                                 -> Matrix e
hermMulAddMatrixWithScales alpha side a b beta c = 
    M.create $ do
        c' <- M.newCopy c
        hermMulAddToMatrixWithScales alpha side a b beta c'
        return c'

-- | @hermMulToVector a x y@ sets @y := a * x@.
hermMulToVector :: (RMatrix m, RVector v, BLAS2 e)
                      => Herm m e
                      -> v e
                      -> STVector s e
                      -> ST s ()
hermMulToVector = hermMulToVectorWithScale 1

-- | @hermMulToVectorWithScale alpha a x y@
-- sets @y := alpha * a * x@.
hermMulToVectorWithScale :: (RMatrix m, RVector v, BLAS2 e)
                               => e
                               -> Herm m e
                               -> v e
                               -> STVector s e
                               -> ST s ()
hermMulToVectorWithScale alpha a x y =
    hermMulAddToVectorWithScales alpha a x 0 y

-- | @hermMulAddToVectorWithScales alpha a x beta y@
-- sets @y := alpha * a * x + beta * y@.
hermMulAddToVectorWithScales :: (RMatrix m, RVector v, BLAS2 e)
                                   => e
                                   -> Herm m e
                                   -> v e
                                   -> e
                                   -> STVector s e
                                   -> ST s ()
hermMulAddToVectorWithScales alpha (Herm uplo a) x beta y
    | ma /= na = error $
        printf ("hermMulAddToVectorWithScales _"
                ++ " (Herm %s <matrix with dim (%d,%d)>)"
                ++ " %s <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>: Herm matrix is not square")
               (show uplo) ma na
               nx ny
               
    | (not . and) [ (ma,na) == (n,n)
                  , nx == n
                  , ny == n
                  ] = error $
        printf ("hermMulAddToVectorWithScales _"
                ++ " (Herm %s <matrix with dim (%d,%d)>)"
                ++ " %s <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>: dimension mismatch")
               (show uplo) ma na
               nx ny

    | otherwise =
        unsafeIOToST $
            M.unsafeWith a $ \pa lda ->
            V.unsafeWith x $ \px ->
            V.unsafeWith y $ \py ->
                BLAS.hemv uplo n alpha pa lda px 1 beta py 1
  where
    (ma,na) = M.dim a
    nx = V.dim x
    ny = V.dim y
    n = ny

-- | @hermMulToMatrix side a b c@
-- sets @c := a * b@ when @side@ is @LeftSide@ and
-- @c := b * a@ when @side@ is @RightSide@.
hermMulToMatrix :: (RMatrix m1, RMatrix m2, BLAS3 e)
                      => Side -> Herm m1 e
                      -> m2 e
                      -> STMatrix s e
                      -> ST s ()
hermMulToMatrix = hermMulToMatrixWithScale 1

-- | @hermMulToMatrixWithScale alpha side a b c@
-- sets @c := alpha * a * b@ when @side@ is @LeftSide@ and
-- @c := alpha * b * a@ when @side@ is @RightSide@.
hermMulToMatrixWithScale :: (RMatrix m1, RMatrix m2, BLAS3 e)
                               => e
                               -> Side -> Herm m1 e
                               -> m2 e
                               -> STMatrix s e
                               -> ST s ()
hermMulToMatrixWithScale alpha side a b c =
    hermMulAddToMatrixWithScales alpha side a b 0 c

-- | @hermMulAddToMatrixWithScales alpha side a b beta c@
-- sets @c := alpha * a * b + beta * c@ when @side@ is @LeftSide@ and
-- @c := alpha * b * a + beta * c@ when @side@ is @RightSide@.
hermMulAddToMatrixWithScales :: (RMatrix m1, RMatrix m2, BLAS3 e)
                                   => e
                                   -> Side -> Herm m1 e
                                   -> m2 e
                                   -> e
                                   -> STMatrix s e
                                   -> ST s ()
hermMulAddToMatrixWithScales alpha side (Herm uplo a) b beta c
    | ma /= na = error $
        printf ("hermMulAddToMatrixWithScales _"
                ++ " %s (Herm %s <matrix with dim (%d,%d)>)" 
                ++ " <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>: Herm matrix is not square")
               (show side) (show uplo) ma na
               mb nb
               mc nc
    | (not . and) [ case side of LeftSide  -> (ma,na) == (m,m)
                                 RightSide -> (ma,na) == (n,n)
                  , (mb, nb ) == (m,n)
                  , (mc, nc ) == (m,n)
                  ] = error $
        printf ("hermMulAddToMatrixWithScales _"
                ++ " %s (Herm %s <matrix with dim (%d,%d)>)" 
                ++ " <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
               (show side) (show uplo) ma na
               mb nb
               mc nc
    | otherwise =
        unsafeIOToST $
            M.unsafeWith a $ \pa lda ->
            M.unsafeWith b $ \pb ldb ->
            M.unsafeWith c $ \pc ldc ->
                BLAS.hemm side uplo m n alpha pa lda pb ldb beta pc ldc
  where
    (ma,na) = M.dim a
    (mb,nb) = M.dim b
    (mc,nc) = M.dim c
    (m,n) = M.dim c
