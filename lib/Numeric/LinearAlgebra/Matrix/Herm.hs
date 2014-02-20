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
    -- * Immutable interface
    
    -- ** Vector multiplication
    hermMulVector,
    hermMulVectorWithScale,
    addHermMulVectorWithScales,
    
    -- ** Matrix  multiplication
    hermMulMatrix,
    hermMulMatrixWithScale,
    addHermMulMatrixWithScales,

    -- ** Updates
    hermRank1Update,
    hermRank2Update,
    hermRankKUpdate,   
    hermRank2KUpdate,    

    
    -- * Mutable interface
    hermCreate,
    
    -- ** Vector multiplication
    hermMulVectorTo,
    hermMulVectorWithScaleTo,
    addHermMulVectorWithScalesM_,
    
    -- ** Matrix multiplication
    hermMulMatrixTo,
    hermMulMatrixWithScaleTo,
    addHermMulMatrixWithScalesM_,

    -- ** Updates
    hermRank1UpdateM_,
    hermRank2UpdateM_,
    hermRankKUpdateM_,  
    hermRank2KUpdateM_,
    
    ) where

import Control.Monad( when )
import Control.Monad.ST( ST, runST )
import Control.Monad.ST.Unsafe( unsafeIOToST )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Vector( Vector, RVector, STVector )
import qualified Numeric.LinearAlgebra.Vector as V
import Numeric.LinearAlgebra.Matrix.Base( Matrix )
import Numeric.LinearAlgebra.Matrix.STBase( STMatrix, RMatrix )
import qualified Numeric.LinearAlgebra.Matrix.STBase as M
import Numeric.LinearAlgebra.Types
import qualified Foreign.BLAS as BLAS


-- | A safe way to create and work with a mutable Herm Matrix before returning 
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
    hermRank1UpdateM_ alpha x (Herm uplo ma')
    a' <- M.unsafeFreeze ma'
    return $ Herm uplo a'

-- | @hermRank2Update alpha x y a@ returns
-- @alpha * x * y^H + conj(alpha) * y * x^H + a@.
hermRank2Update :: (BLAS2 e)
                => e -> Vector e -> Vector e -> Herm Matrix e
                -> Herm Matrix e
hermRank2Update alpha x y (Herm uplo a) = runST $ do
    ma' <- M.newCopy a
    hermRank2UpdateM_ alpha x y (Herm uplo ma')
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
    hermRankKUpdateM_ alpha trans a beta (Herm uplo mc')
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
    hermRank2KUpdateM_ alpha trans a b beta (Herm uplo mc')
    c' <- M.unsafeFreeze mc'
    return $ Herm uplo c'

-- | @hermRank1UpdateM_ alpha x a@ sets
-- @a := alpha * x * x^H + a@.
hermRank1UpdateM_ :: (RVector v, BLAS2 e)
                  => Double -> v e -> Herm (STMatrix s) e -> ST s ()
hermRank1UpdateM_ alpha x (Herm uplo a) = do
    nx <- V.getDim x
    (ma,na) <- M.getDim a
    let n = nx

    when ((not . and) [ nx == n, (ma,na) == (n,n) ]) $ error $
        printf ("hermRank1UpdateM_ _ <vector with dim %d>"
                 ++ " (Herm _ <matrix with dim (%d,%d)>):"
                 ++ " invalid dimensions") nx ma na

    unsafeIOToST $
        V.unsafeWith x $ \px ->
        M.unsafeWith a $ \pa lda ->
            BLAS.her uplo n alpha px 1 pa lda


-- | @hermRank2UpdateM_ alpha x y a@ sets
-- @a := alpha * x * y^H + conj(alpha) * y * x^H + a@.
hermRank2UpdateM_ :: (RVector v1, RVector v2, BLAS2 e)
                  => e -> v1 e -> v2 e -> Herm (STMatrix s) e -> ST s ()
hermRank2UpdateM_ alpha x y (Herm uplo a) = do
    nx <- V.getDim x
    ny <- V.getDim y
    (ma,na) <- M.getDim a
    let n = nx
    
    when ((not . and) [ nx == n, ny == n, (ma,na) == (n,n) ]) $ error $
        printf ("hermRank2UpdateM_ _ <vector with dim %d>"
                 ++ " <vector with dim %d>"
                 ++ " (Herm _ <matrix with dim (%d,%d)>):"
                 ++ " invalid dimensions") nx ny ma na

    unsafeIOToST $
        V.unsafeWith x $ \px ->
        V.unsafeWith y $ \py ->
        M.unsafeWith a $ \pa lda ->
            BLAS.her2 uplo n alpha px 1 py 1 pa lda


-- | @hermRankKUpdateM_ alpha trans a beta c@ sets
-- @c := alpha * a * a^H + beta * c@ when @trans@ is @NoTrans@ and
-- @c := alpha * a^H * a + beta * c@ when @trans@ is @ConjTrans@.  The
-- function signals an error when @trans@ is @Trans@.
hermRankKUpdateM_ :: (RMatrix m, BLAS3 e)
                  => e -> Trans -> m e -> e -> Herm (STMatrix s) e
                  -> ST s ()
hermRankKUpdateM_ alpha trans a beta (Herm uplo c) = do
    (ma,na) <- M.getDim a
    (mc,nc) <- M.getDim c
    let (n,k) = if trans == NoTrans then (ma,na) else (na,ma)

    when (trans == Trans) $ error $
        printf ("hermRankKUpdateM_ _ %s:"
                 ++ " trans argument must be NoTrans or ConjTrans")
               (show trans)
               
    when ((not . and) [ (mc,nc) == (n,n)
                      , case trans of NoTrans -> (ma,na) == (n,k)
                                      _       -> (ma,na) == (k,n)
                      ]) $ error $
            printf ("hermRankKUpdateM_ _ %s <matrix with dim (%d,%d)> _"
                    ++ " (Herm _ <matrix with dim (%d,%d)>):"
                    ++ " invalid dimensions") (show trans) ma na mc nc

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        M.unsafeWith c $ \pc ldc ->
            BLAS.herk uplo trans n k alpha pa lda beta pc ldc


-- | @hermRank2KUpdateM_ alpha trans a b beta c@ sets
-- @c := alpha * a * b^H + conj(alpha) * b * a^H + beta * c@ when @trans@ is
-- @NoTrans@ and @c := alpha * b^H * a + conj(alpha) * a^H * b + beta * c@
-- when @trans@ is @ConjTrans@.  The function signals an error when @trans@
-- is @Trans@.
hermRank2KUpdateM_ :: (RMatrix m1, RMatrix m2, BLAS3 e)
                   => e -> Trans -> m1 e -> m2 e -> e -> Herm (STMatrix s) e
                   -> ST s ()
hermRank2KUpdateM_ alpha trans a b beta (Herm uplo c) = do
    (ma,na) <- M.getDim a
    (mb,nb) <- M.getDim b
    (mc,nc) <- M.getDim c
    let (n,k) = if trans == NoTrans then (ma,na) else (na,ma)

    when (trans == Trans) $ error $
        printf ("hermRank2KUpdateM_ _ %s:"
                 ++ " trans argument must be NoTrans or ConjTrans")
               (show trans)

    when ((not . and) [ (mc,nc) == (n,n)
                      , (mb,nb) == (ma,na)
                      , case trans of NoTrans -> (ma,na) == (n,k)
                                      _       -> (ma,na) == (k,n)
                      ]) $ error $
            printf ("hermRank2KUpdateM_ _ %s <matrix with dim (%d,%d)>"
                    ++ " <matrix with dim (%d,%d)> _"
                    ++ " (Herm _ <matrix with dim (%d,%d)>):"
                    ++ " invalid dimensions") (show trans) ma na mb nb mc nc

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        M.unsafeWith b $ \pb ldb ->
        M.unsafeWith c $ \pc ldc ->
            BLAS.her2k uplo trans n k alpha pa lda pb ldb beta pc ldc


-- | @hermMulVector a x@ returns @a * x@.
hermMulVector :: (BLAS2 e)
              => Herm Matrix e
              -> Vector e
              -> Vector e
hermMulVector a x =
    V.create $ do
        n <- V.getDim x
        y <- V.new_ n
        hermMulVectorTo y a x
        return y

-- | @hermMulVectorWithScale alpha a x@ retunrs @alpha * a * x@.
hermMulVectorWithScale :: (BLAS2 e)
                       => e
                       -> Herm Matrix e
                       -> Vector e
                       -> Vector e
hermMulVectorWithScale alpha a x =
    V.create $ do
        n <- V.getDim x
        y <- V.new_ n
        hermMulVectorWithScaleTo y alpha a x
        return y
                       
-- | @addHermMulVectorWithScales alpha a x y@
-- returns @alpha * a * x + beta * y@.
addHermMulVectorWithScales :: (BLAS2 e)
                           => e
                           -> Herm Matrix e
                           -> Vector e
                           -> e
                           -> Vector e
                           -> Vector e
addHermMulVectorWithScales alpha a x beta y =
    V.create $ do
        y' <- V.newCopy y
        addHermMulVectorWithScalesM_ alpha a x beta y'
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
        mn <- M.getDim b
        c <- M.new_ mn
        hermMulMatrixTo c side a b
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
        mn <- M.getDim b
        c <- M.new_ mn
        hermMulMatrixWithScaleTo c alpha side a b
        return c

-- | @addHermMulMatrixWithScales alpha side a b beta c@
-- returns @alpha * a * b + beta * c@ when @side@ is @LeftSide@ and
-- @alpha * b * a + beta * c@ when @side@ is @RightSide@.
addHermMulMatrixWithScales :: (BLAS3 e)
                           => e
                           -> Side -> Herm Matrix e
                           -> Matrix e
                           -> e
                           -> Matrix e
                           -> Matrix e
addHermMulMatrixWithScales alpha side a b beta c = 
    M.create $ do
        c' <- M.newCopy c
        addHermMulMatrixWithScalesM_ alpha side a b beta c'
        return c'

-- | @hermMulVectorTo dst a x@ sets @dst := a * x@.
hermMulVectorTo :: (RMatrix m, RVector v, BLAS2 e)
                => STVector s e
                -> Herm m e
                -> v e
                -> ST s ()
hermMulVectorTo dst = hermMulVectorWithScaleTo dst 1

-- | @hermMulVectorWithScaleTo dst alpha a x@
-- sets @dst := alpha * a * x@.
hermMulVectorWithScaleTo :: (RMatrix m, RVector v, BLAS2 e)
                         => STVector s e
                         -> e
                         -> Herm m e
                         -> v e
                         -> ST s ()
hermMulVectorWithScaleTo dst alpha a x =
    addHermMulVectorWithScalesM_ alpha a x 0 dst

-- | @addHermMulVectorWithScalesM_ alpha a x beta y@
-- sets @y := alpha * a * x + beta * y@.
addHermMulVectorWithScalesM_ :: (RMatrix m, RVector v, BLAS2 e)
                             => e
                             -> Herm m e
                             -> v e
                             -> e
                             -> STVector s e
                             -> ST s ()
addHermMulVectorWithScalesM_ alpha (Herm uplo a) x beta y = do
    (ma,na) <- M.getDim a
    nx <- V.getDim x
    ny <- V.getDim y
    let n = ny
    
    when (ma /= na) $ error $
        printf ("addHermMulVectorWithScalesM_ _"
                ++ " (Herm %s <matrix with dim (%d,%d)>)"
                ++ " %s <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>: Herm matrix is not square")
               (show uplo) ma na
               nx ny
               
    when ((not . and) [ (ma,na) == (n,n)
                      , nx == n
                      , ny == n
                      ]) $ error $
        printf ("addHermMulVectorWithScalesM_ _"
                ++ " (Herm %s <matrix with dim (%d,%d)>)"
                ++ " %s <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>: dimension mismatch")
               (show uplo) ma na
               nx ny

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        V.unsafeWith x $ \px ->
        V.unsafeWith y $ \py ->
            BLAS.hemv uplo n alpha pa lda px 1 beta py 1

-- | @hermMulMatrixTo dst side a b@
-- sets @dst := a * b@ when @side@ is @LeftSide@ and
-- @dst := b * a@ when @side@ is @RightSide@.
hermMulMatrixTo :: (RMatrix m1, RMatrix m2, BLAS3 e)
                => STMatrix s e
                -> Side -> Herm m1 e
                -> m2 e
                -> ST s ()
hermMulMatrixTo dst = hermMulMatrixWithScaleTo dst 1

-- | @hermMulMatrixWithScaleTo dst alpha side a b@
-- sets @dst := alpha * a * b@ when @side@ is @LeftSide@ and
-- @dst := alpha * b * a@ when @side@ is @RightSide@.
hermMulMatrixWithScaleTo :: (RMatrix m1, RMatrix m2, BLAS3 e)
                         => STMatrix s e
                         -> e
                         -> Side -> Herm m1 e
                         -> m2 e
                         -> ST s ()
hermMulMatrixWithScaleTo dst alpha side a b =
    addHermMulMatrixWithScalesM_ alpha side a b 0 dst

-- | @addHermMulMatrixWithScalesM_ alpha side a b beta c@
-- sets @c := alpha * a * b + beta * c@ when @side@ is @LeftSide@ and
-- @c := alpha * b * a + beta * c@ when @side@ is @RightSide@.
addHermMulMatrixWithScalesM_ :: (RMatrix m1, RMatrix m2, BLAS3 e)
                             => e
                             -> Side -> Herm m1 e
                             -> m2 e
                             -> e
                             -> STMatrix s e
                             -> ST s ()
addHermMulMatrixWithScalesM_ alpha side (Herm uplo a) b beta c = do
    (ma,na) <- M.getDim a
    (mb,nb) <- M.getDim b
    (mc,nc) <- M.getDim c
    let (m,n) = (mc,nc)
    
    when (ma /= na) $ error $
        printf ("addHermMulMatrixWithScalesM_ _"
                ++ " %s (Herm %s <matrix with dim (%d,%d)>)" 
                ++ " <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>: Herm matrix is not square")
               (show side) (show uplo) ma na
               mb nb
               mc nc
    when ((not . and) [ case side of LeftSide  -> (ma,na) == (m,m)
                                     RightSide -> (ma,na) == (n,n)
                      , (mb, nb ) == (m,n)
                      , (mc, nc ) == (m,n)
                      ]) $ error $
        printf ("addHermMulMatrixWithScalesM_ _"
                ++ " %s (Herm %s <matrix with dim (%d,%d)>)" 
                ++ " <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
               (show side) (show uplo) ma na
               mb nb
               mc nc

    unsafeIOToST $
        M.unsafeWith a $ \pa lda ->
        M.unsafeWith b $ \pb ldb ->
        M.unsafeWith c $ \pc ldc ->
            BLAS.hemm side uplo m n alpha pa lda pb ldb beta pc ldc
