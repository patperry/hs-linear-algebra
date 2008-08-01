{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Herm.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Herm.Dense (
    module Data.Matrix.Herm,
    module BLAS.Matrix.Immutable,
    module BLAS.Matrix.ReadOnly,

    hemv,
    hemm,
    unsafeHemv,
    unsafeHemm,
    ) where

import Control.Monad ( zipWithM_ )
import System.IO.Unsafe
import Unsafe.Coerce

import BLAS.Access
import BLAS.Elem ( BLAS2, BLAS3 )
import BLAS.Internal ( checkMatVecMultAdd, checkMatMatMultAdd )
import BLAS.C ( colMajor, rightSide, leftSide, cblasUpLo )
import BLAS.Types ( flipUpLo )
import qualified BLAS.Elem as E
import qualified BLAS.C as BLAS

import Data.Matrix.Dense.Internal
import Data.Vector.Dense.Internal
import qualified Data.Matrix.Dense.Internal as M
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V

import Data.Matrix.Herm
import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly


instance (BLAS3 e) => IMatrix (Herm (DMatrix Imm)) e where
    unsafeSApply k h x = unsafePerformIO $ unsafeGetSApply k h x
    {-# NOINLINE unsafeSApply #-}

    unsafeSApplyMat k h a = unsafePerformIO $ unsafeGetSApplyMat k h a
    {-# NOINLINE unsafeSApplyMat #-}


instance (BLAS3 e) => RMatrix (Herm (DMatrix s)) e where
    unsafeGetSApply k h x = do
        y <- newZero (dim x)
        unsafeHemv k (unsafeCoerce h) x 1 y
        return (unsafeCoerce y)
    
    unsafeGetSApplyMat k h a = do
        b <- newZero (shape a)
        unsafeHemm k (unsafeCoerce h) a 1 b
        return (unsafeCoerce b)

hemv :: (BLAS2 e) => e -> Herm (DMatrix t) (n,n) e -> DVector s n e -> e -> IOVector n e -> IO ()
hemv alpha h x beta y =
    checkMatVecMultAdd (numRows h, numCols h) (dim x) (dim y) $
        unsafeHemv alpha h x beta y

unsafeHemv :: (BLAS2 e) => e -> Herm (DMatrix t) (n,n) e -> DVector s n e -> e -> IOVector n e -> IO ()
unsafeHemv alpha h x beta y
    | numRows h == 0 =
        return ()
    | isConj y = do
        V.doConj y
        hemv alpha h x beta (V.conj y)
        V.doConj y
    | isConj x = do
        x' <- newCopy x
        V.doConj (V.unsafeThaw x')
        hemv alpha h (conj x') beta y
    | otherwise =
        let order   = colMajor
            (u,e,a) = toBase h
            alpha'  = alpha * e
            n       = numCols a
            (u',alpha'') 
                    = case (isHerm a) of
                          True  -> (flipUpLo u, E.conj alpha')
                          False -> (u, alpha')
            uploA   = cblasUpLo u'
            ldA     = ldaOf a
            incX    = strideOf x
            incY    = strideOf y
        in M.unsafeWithElemPtr a (0,0) $ \pA ->
               V.unsafeWithElemPtr x 0 $ \pX ->
                    V.unsafeWithElemPtr y 0 $ \pY -> do
                        BLAS.hemv order uploA n alpha'' pA ldA pX incX beta pY incY

hemm :: (BLAS3 e) => e -> Herm (DMatrix t) (m,m) e -> DMatrix s (m,n) e -> e -> IOMatrix (m,n) e -> IO ()
hemm alpha h b beta c =
    checkMatMatMultAdd (numRows h, numCols h) (shape b) (shape c) $
        unsafeHemm alpha h b beta c

unsafeHemm :: (BLAS3 e) => e -> Herm (DMatrix t) (m,m) e -> DMatrix s (m,n) e -> e -> IOMatrix (m,n) e -> IO ()
unsafeHemm alpha h b beta c
    | numRows b == 0 || numCols b == 0 || numCols c == 0 = return ()
    | (isHerm a) /= (isHerm c) || (isHerm a) /= (isHerm b) =
        zipWithM_ (\x y -> hemv alpha h x beta y) (cols b) (cols c)
    | otherwise =
        let order   = colMajor
            (m,n)   = shape c
            alpha'  = alpha * e
            (side,u',m',n', alpha'')
                    = if isHerm a
                          then (rightSide, flipUpLo u, n, m, E.conj alpha')
                          else (leftSide,  u, m, n, alpha')
            uploA   = cblasUpLo u'
            ldA     = ldaOf a
            ldB     = ldaOf b
            ldC     = ldaOf c
        in M.unsafeWithElemPtr a (0,0) $ \pA ->
               M.unsafeWithElemPtr b (0,0) $ \pB ->
                   M.unsafeWithElemPtr c (0,0) $ \pC ->
                       BLAS.hemm order side uploA m' n' alpha'' pA ldA pB ldB beta pC ldC
    where
      (u,e,a) = toBase h
