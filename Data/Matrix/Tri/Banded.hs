{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Tri.Banded
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Tri.Banded (
    module Data.Matrix.Tri,
    module BLAS.Matrix.Immutable,
    module BLAS.Matrix.ReadOnly,
    module BLAS.Matrix.Solve,

    tbmv,
    tbsv,
    tbmm,    
    tbsm
    ) where

-- import Control.Monad ( when )
import System.IO.Unsafe
import Unsafe.Coerce

import Data.Matrix.Banded.Internal
import Data.Matrix.Dense.Internal ( IOMatrix )
import qualified Data.Matrix.Dense.Internal as M
import Data.Vector.Dense.Internal
import qualified Data.Vector.Dense.Internal as V
-- import qualified Data.Vector.Dense.Operations as V

import BLAS.Access
import BLAS.Elem ( BLAS2 )
-- import qualified BLAS.Elem as E
-- import BLAS.C.Types ( cblasDiag, cblasUpLo, cblasTrans, colMajor, 
--    noTrans, conjTrans, leftSide, rightSide )
-- import BLAS.Types ( Trans(..), flipTrans, flipUpLo )
                                   
-- import qualified BLAS.C as BLAS
-- import qualified BLAS.C.Types as BLAS

import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly
import BLAS.Matrix.Solve

import Data.Matrix.Tri

instance (BLAS2 e) => IMatrix (Tri (BMatrix Imm)) e where
    (<*>) t x = unsafePerformIO $ getApply t x
    {-# NOINLINE (<*>) #-}

    (<**>) t a = unsafePerformIO $ getApplyMat t a
    {-# NOINLINE (<**>) #-}

instance (BLAS2 e) => ISolve (Tri (BMatrix Imm)) e where
    (<\>) t x = unsafePerformIO $ getSolve t x
    {-# NOINLINE (<\>) #-}

    (<\\>) t a = unsafePerformIO $ getSolveMat t a
    {-# NOINLINE (<\\>) #-}

instance (BLAS2 e) => RMatrix (Tri (BMatrix s)) e where
    getApply t x = do
        x' <- newCopy x
        tbmv (unsafeCoerce t) (V.unsafeThaw x')
        return (unsafeCoerce x')
    
    getApplyMat t a = do
        a' <- newCopy a
        tbmm (unsafeCoerce t) (M.unsafeThaw a')
        return (unsafeCoerce a')

instance (BLAS2 e) => RSolve (Tri (BMatrix s)) e where
    getSolve t x = do
        x' <- newCopy x
        tbsv (unsafeCoerce t) (V.unsafeThaw x')
        return (unsafeCoerce x')
    
    getSolveMat t a = do
        a' <- newCopy a
        tbsm (unsafeCoerce t) (M.unsafeThaw a')
        return (unsafeCoerce a')


tbmv :: (BLAS2 e) => Tri (BMatrix t) (n,n) e -> IOVector n e -> IO ()
tbmv = undefined {- t x =
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        (transA,u') = if isHerm a then (conjTrans, flipUpLo u) else (noTrans, u)
        uploA     = cblasUpLo u'
        diagA     = cblasDiag d
        n         = numCols a
        ldA       = ldaOf a
        incX      = strideOf x
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           V.unsafeWithElemPtr x 0 $ \pX -> do
               BLAS.trmv order uploA transA diagA n pA ldA pX incX
               when (alpha /= 1) $ V.scaleBy alpha x
-}               
               
tbmm :: (BLAS2 e) => Tri (BMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
tbmm t b = mapM_ (\x -> tbmv t x) (M.cols b)

tbsv :: (BLAS2 e) =>Tri (BMatrix t) (n,n) e -> IOVector n e -> IO ()
tbsv = undefined {- t x =
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        (transA,u') = if isHerm a then (conjTrans, flipUpLo u) else (noTrans, u)
        uploA     = cblasUpLo u'
        diagA     = cblasDiag d
        n         = numCols a
        ldA       = ldaOf a
        incX      = strideOf x
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           V.unsafeWithElemPtr x 0 $ \pX -> do
               BLAS.trsv order uploA transA diagA n pA ldA pX incX
               when (alpha /= 1) $ V.invScaleBy alpha x
-}

tbsm :: (BLAS2 e) => Tri (BMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
tbsm t b = mapM_ (\x -> tbsv t x) (M.cols b)
