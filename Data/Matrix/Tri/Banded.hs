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

import Control.Monad ( when )
import System.IO.Unsafe
import Unsafe.Coerce

import Data.Matrix.Banded.Internal
import qualified Data.Matrix.Banded.Internal as B
import Data.Matrix.Dense.Internal ( IOMatrix )
import qualified Data.Matrix.Dense.Internal as M
import Data.Vector.Dense.Internal
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V

import BLAS.Access
import BLAS.Types ( flipUpLo )

import BLAS.C ( BLAS2, cblasDiag, cblasUpLo, colMajor, noTrans, conjTrans )                                   
import qualified BLAS.C as BLAS

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
tbmv t x | isConj x = do
    V.doConj x
    tbmv t (conj x)
    V.doConj x
    
tbmv t x =
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        (transA,u') = if isHerm a then (conjTrans, flipUpLo u) else (noTrans, u)
        uploA     = cblasUpLo u'
        diagA     = cblasDiag d
        n         = numCols a
        k         = case u of Upper -> numUpper a 
                              Lower -> numLower a
        ldA       = ldaOf a
        incX      = strideOf x
        withPtrA  = case u' of 
                        Upper -> B.unsafeWithBasePtr a
                        Lower -> B.unsafeWithElemPtr a (0,0)
    in withPtrA $ \pA ->
           V.unsafeWithElemPtr x 0 $ \pX -> do
               BLAS.tbmv order uploA transA diagA n k pA ldA pX incX
               when (alpha /= 1) $ V.scaleBy alpha x
              
               
tbmm :: (BLAS2 e) => Tri (BMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
tbmm t b = mapM_ (\x -> tbmv t x) (M.cols b)

tbsv :: (BLAS2 e) => Tri (BMatrix t) (n,n) e -> IOVector n e -> IO ()
tbsv t x | isConj x = do
    V.doConj x
    tbsv t (conj x)
    V.doConj x
    
tbsv t x = 
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        (transA,u') = if isHerm a then (conjTrans, flipUpLo u) else (noTrans, u)
        uploA     = cblasUpLo u'
        diagA     = cblasDiag d
        n         = numCols a
        k         = case u of Upper -> numUpper a 
                              Lower -> numLower a        
        ldA       = ldaOf a
        incX      = strideOf x
        withPtrA  = case u' of 
                        Upper -> B.unsafeWithBasePtr a
                        Lower -> B.unsafeWithElemPtr a (0,0)
    in withPtrA $ \pA ->
           V.unsafeWithElemPtr x 0 $ \pX -> do
               BLAS.tbsv order uploA transA diagA n k pA ldA pX incX
               when (alpha /= 1) $ V.invScaleBy alpha x


tbsm :: (BLAS2 e) => Tri (BMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
tbsm t b = mapM_ (\x -> tbsv t x) (M.cols b)
