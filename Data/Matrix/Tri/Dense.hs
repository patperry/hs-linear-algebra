{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Tri.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Tri.Dense (
    module Data.Matrix.Tri,
    module BLAS.Matrix.Immutable,
    module BLAS.Matrix.ReadOnly,
    module BLAS.Matrix.Solve,

    trmv,
    trsv,
    trmm,    
    trsm
    ) where

import Control.Monad ( when )
import Data.Maybe ( fromJust )
import System.IO.Unsafe
import Unsafe.Coerce

import Data.Matrix.Dense.Internal
import qualified Data.Matrix.Dense.Internal as M
import Data.Vector.Dense.Internal
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V

import BLAS.Access
import BLAS.Elem ( BLAS3 )
import qualified BLAS.Elem as E
import BLAS.C.Types ( cblasDiag, cblasUpLo, cblasTrans, colMajor, 
    noTrans, conjTrans, leftSide, rightSide )
import BLAS.Types ( Trans(..), flipTrans )
                                   
import qualified BLAS.C as BLAS
import qualified BLAS.C.Types as BLAS

import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly
import BLAS.Matrix.Solve

import Data.Matrix.Tri

instance (BLAS3 e) => IMatrix (Tri (DMatrix Imm)) e where
    (<*>) t x = unsafePerformIO $ getApply t x
    {-# NOINLINE (<*>) #-}

    (<**>) t a = unsafePerformIO $ getApplyMat t a
    {-# NOINLINE (<**>) #-}

instance (BLAS3 e) => ISolve (Tri (DMatrix Imm)) e where
    (<\>) t x = unsafePerformIO $ getSolve t x
    {-# NOINLINE (<\>) #-}

    (<\\>) t a = unsafePerformIO $ getSolveMat t a
    {-# NOINLINE (<\\>) #-}

instance (BLAS3 e) => RMatrix (Tri (DMatrix s)) e where
    getApply t x = do
        x' <- newCopy x
        trmv (unsafeCoerce t) (V.unsafeThaw x')
        return (unsafeCoerce x')
    
    getApplyMat t a = do
        a' <- newCopy a
        trmm (unsafeCoerce t) (M.unsafeThaw a')
        return (unsafeCoerce a')

instance (BLAS3 e) => RSolve (Tri (DMatrix s)) e where
    getSolve t x = do
        x' <- newCopy x
        trsv (unsafeCoerce t) (V.unsafeThaw x')
        return (unsafeCoerce x')
    
    getSolveMat t a = do
        a' <- newCopy a
        trsm (unsafeCoerce t) (M.unsafeThaw a')
        return (unsafeCoerce a')


trmv :: (BLAS3 e) => Tri (DMatrix t) (n,n) e -> IOVector n e -> IO ()
trmv _ x
    | dim x == 0 = return ()
trmv t x
    | isConj x =
        let b = fromJust $ maybeFromCol x
        in trmm t b
trmv t x =
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        uploA     = cblasUpLo u
        transA    = if isHerm a then conjTrans else noTrans
        diagA     = cblasDiag d
        n         = numCols a
        ldA       = ldaOf a
        incX      = strideOf x
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           V.unsafeWithElemPtr x 0 $ \pX -> do
               BLAS.trmv order uploA transA diagA n pA ldA pX incX
               when (alpha /= 1) $ V.scaleBy alpha x
               
               
trmm :: (BLAS3 e) => Tri (DMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
trmm _ b
    | M.numRows b == 0 || M.numCols b == 0 = return ()
trmm t b =
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        h         = if isHerm a then ConjTrans else NoTrans
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if M.isHerm b
                        then (rightSide, flipTrans h, n, m, E.conj alpha)
                        else (leftSide , h          , m, n, alpha       )
        uploA     = cblasUpLo u
        transA    = cblasTrans h'
        diagA     = cblasDiag d
        ldA       = ldaOf a
        ldB       = ldaOf b
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           M.unsafeWithElemPtr b (0,0) $ \pB ->
               BLAS.trmm order side uploA transA diagA m' n' alpha' pA ldA pB ldB

trsv :: (BLAS3 e) =>Tri (DMatrix t) (n,n) e -> IOVector n e -> IO ()
trsv _ x
    | dim x == 0 = return ()
trsv t x
    | isConj x =
        let b = fromJust $ maybeFromCol x
        in trsm t b
trsv t x =
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        uploA     = cblasUpLo u
        transA    = if isHerm a then conjTrans else noTrans
        diagA     = cblasDiag d
        n         = numCols a
        ldA       = ldaOf a
        incX      = strideOf x
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           V.unsafeWithElemPtr x 0 $ \pX -> do
               BLAS.trsv order uploA transA diagA n pA ldA pX incX
               when (alpha /= 1) $ V.invScaleBy alpha x

trsm :: (BLAS3 e) => Tri (DMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
trsm _ b
    | M.numRows b == 0 || M.numCols b == 0 = return ()
trsm t b =
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        h         = if isHerm a then ConjTrans else NoTrans
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if isHerm b
                        then (rightSide, flipTrans h, n, m, E.conj alpha)
                        else (leftSide , h          , m, n, alpha       )
        uploA     = cblasUpLo u
        transA    = cblasTrans h'
        diagA     = cblasDiag d
        ldA       = ldaOf a
        ldB       = ldaOf b
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           M.unsafeWithElemPtr b (0,0) $ \pB -> do
               BLAS.trsm order side uploA transA diagA m' n' (1/alpha') pA ldA pB ldB
               