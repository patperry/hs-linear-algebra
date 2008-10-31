{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Herm
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Herm (
    Herm(..),
    UpLo(..),

    fromBase,
    toBase,
    mapHerm,

    hermL,
    hermU,

    coerceHerm,

    ) where

import Control.Monad( zipWithM_ )
import Control.Monad.ST( ST )
import Unsafe.Coerce

import BLAS.C( BLAS2, BLAS3, colMajor, rightSide, leftSide, cblasUpLo )
import qualified BLAS.C as BLAS
import BLAS.UnsafeIOToM

import BLAS.Matrix
import BLAS.Types ( UpLo(..), flipUpLo )

import Data.Matrix.Banded( Banded )
import Data.Matrix.Banded.Class 
import Data.Matrix.Banded.IO( IOBanded )
import Data.Matrix.Banded.ST( STBanded )
import Data.Matrix.Dense( Matrix )
import Data.Matrix.Dense.Class hiding ( BaseMatrix )
import Data.Matrix.Dense.IO( IOMatrix )
import Data.Matrix.Dense.ST( STMatrix, runSTMatrix )
import Data.Vector.Dense.Class
import Data.Vector.Dense.ST( runSTVector )


data Herm a nn e = Herm UpLo (a nn e)

coerceHerm :: Herm a mn e -> Herm a mn' e
coerceHerm = unsafeCoerce

mapHerm :: (a (n,n) e -> b (n,n) e) -> Herm a (n,n) e -> Herm b (n,n) e
mapHerm f (Herm u a) = Herm u $ f a

fromBase :: UpLo -> a (n,n) e -> Herm a (n,n) e
fromBase = Herm
        
toBase :: Herm a (n,n) e -> (UpLo, a (n,n) e)
toBase (Herm u a) = (u,a)

hermL :: a (n,n) e -> Herm a (n,n) e
hermL = Herm Lower

hermU :: a (n,n) e -> Herm a (n,n) e
hermU = Herm Upper
      
instance BaseMatrix a e => BaseTensor (Herm a) (Int,Int) e where
    shape  (Herm _ a) = (n,n)             where n = min (numRows a) (numCols a)
    bounds (Herm _ a) = ((0,0),(n-1,n-1)) where n = min (numRows a) (numCols a)
      
instance BaseMatrix a e => BaseMatrix (Herm a) e where
    herm = coerceHerm
    
instance Show (a mn e) => Show (Herm a mn e) where
    show (Herm u a) = constructor ++ " (" ++ show a ++ ")"
      where
        constructor = case u of
            Lower -> "hermL"
            Upper -> "hermU"


------------------------- Dense Matrix instances ----------------------------

hemv :: (ReadMatrix a z e m, ReadVector x e m, WriteVector y e m, BLAS2 e) => 
    e -> Herm a (k,k) e -> x k e -> e -> y k e -> m ()
hemv alpha h x beta y
    | numRows h == 0 =
        return ()
    | isConj y = do
        doConj y
        hemv alpha h x beta (conj y)
        doConj y
    | isConj x = do
        x' <- newCopyVector x
        doConj x'
        hemv alpha h (conj x') beta y
    | otherwise =
        let order = colMajor
            (u,a) = toBase h
            n     = numCols a
            u'    = case isHermMatrix a of
                        True  -> flipUpLo u
                        False -> u
            uploA = cblasUpLo u'
            ldA   = ldaOfMatrix a
            incX  = stride x
            incY  = stride y
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withVectorPtr x $ \pX ->
               withVectorPtr y $ \pY ->
                   BLAS.hemv order uploA n alpha pA ldA pX incX beta pY incY

hemm :: (ReadMatrix a x e m, ReadMatrix b y e m, WriteMatrix c z e m, BLAS3 e) => 
    e -> Herm a (k,k) e -> b (k,l) e -> e -> c (k,l) e -> m ()
hemm alpha h b beta c
    | numRows b == 0 || numCols b == 0 || numCols c == 0 = return ()
    | (isHermMatrix a) /= (isHermMatrix c) || (isHermMatrix a) /= (isHermMatrix b) =
        zipWithM_ (\x y -> hemv alpha h x beta y) (colViews b) (colViews c)
    | otherwise =
        let order   = colMajor
            (m,n)   = shape c
            (side,u',m',n')
                    = if isHermMatrix a
                          then (rightSide, flipUpLo u, n, m)
                          else (leftSide,  u,          m, n)
            uploA   = cblasUpLo u'
            ldA     = ldaOfMatrix a
            ldB     = ldaOfMatrix b
            ldC     = ldaOfMatrix c
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withMatrixPtr b $ \pB ->
               withMatrixPtr c $ \pC ->
                   BLAS.hemm order side uploA m' n' alpha pA ldA pB ldB beta pC ldC
    where
      (u,a) = toBase h

hemv' :: (ReadMatrix a z e m, ReadVector x e m, WriteVector y e m, BLAS2 e) => 
    e -> Herm a (r,s) e -> x s e -> e -> y r e -> m ()
hemv' alpha a x beta y = 
    hemv alpha (coerceHerm a) x beta (coerceVector y)

hemm' :: (ReadMatrix a x e m, ReadMatrix b y e m, WriteMatrix c z e m, BLAS3 e) => 
    e -> Herm a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
hemm' alpha a b beta c = 
    hemm alpha (coerceHerm a) b beta (coerceMatrix c)

instance (BLAS3 e) => IMatrix (Herm Matrix) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    

instance (BLAS3 e) => MMatrix (Herm (STMatrix s)) e (ST s) where
    unsafeDoSApplyAdd    = hemv'
    unsafeDoSApplyAddMat = hemm'

instance (BLAS3 e) => MMatrix (Herm IOMatrix) e IO where
    unsafeDoSApplyAdd    = hemv'
    unsafeDoSApplyAddMat = hemm'

instance (BLAS3 e, UnsafeIOToM m) => MMatrix (Herm Matrix) e m where
    unsafeDoSApplyAdd    = hemv'
    unsafeDoSApplyAddMat = hemm'


------------------------- Banded Matrix instances ----------------------------

hbmv :: (ReadBanded a z e m, ReadVector x e m, WriteVector y e m, BLAS2 e) => 
    e -> Herm a (k,k) e -> x k e -> e -> y k e -> m ()
hbmv alpha h x beta y
    | numRows h == 0 =
        return ()
    | isConj y = do
        doConj y
        hbmv alpha h x beta (conj y)
        doConj y
    | isConj x = do
        x' <- newCopyVector x
        doConj x'
        hbmv alpha h (conj x') beta y
    | otherwise =
        let order = colMajor
            (u,a) = toBase h
            n     = numCols a
            k     = case u of 
                        Upper -> numUpper a
                        Lower -> numLower a      
            u'    = case (isHermBanded a) of
                        True  -> flipUpLo u
                        False -> u
            uploA = cblasUpLo u'
            ldA   = ldaOfBanded a
            incX  = stride x
            incY  = stride y
            withPtrA 
                  = case u' of Upper -> withBandedPtr a
                               Lower -> withBandedElemPtr a (0,0)
        in unsafeIOToM $
               withPtrA $ \pA ->
               withVectorPtr x $ \pX ->
               withVectorPtr y $ \pY -> do
                   BLAS.hbmv order uploA n k alpha pA ldA pX incX beta pY incY

hbmm :: (ReadBanded a x e m, ReadMatrix b y e m, WriteMatrix c z e m, BLAS2 e) => 
    e -> Herm a (k,k) e -> b (k,l) e -> e -> c (k,l) e -> m ()
hbmm alpha h b beta c =
    zipWithM_ (\x y -> hbmv alpha h x beta y) (colViews b) (colViews c)

hbmv' :: (ReadBanded a z e m, ReadVector x e m, WriteVector y e m, BLAS2 e) => 
    e -> Herm a (r,s) e -> x s e -> e -> y r e -> m ()
hbmv' alpha a x beta y = 
    hbmv alpha (coerceHerm a) x beta (coerceVector y)

hbmm' :: (ReadBanded a x e m, ReadMatrix b y e m, WriteMatrix c z e m, BLAS3 e) => 
    e -> Herm a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
hbmm' alpha a b beta c = 
    hbmm alpha (coerceHerm a) b beta (coerceMatrix c)

instance (BLAS3 e) => IMatrix (Herm Banded) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    

instance (BLAS3 e) => MMatrix (Herm (STBanded s)) e (ST s) where
    unsafeDoSApplyAdd    = hbmv'
    unsafeDoSApplyAddMat = hbmm'

instance (BLAS3 e) => MMatrix (Herm IOBanded) e IO where
    unsafeDoSApplyAdd    = hbmv'
    unsafeDoSApplyAddMat = hbmm'

instance (BLAS3 e, UnsafeIOToM m) => MMatrix (Herm Banded) e m where
    unsafeDoSApplyAdd    = hbmv'
    unsafeDoSApplyAddMat = hbmm'

