{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
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
    trsm,
    
    unsafeTrmv,
    unsafeTrsv,
    unsafeTrmm,
    unsafeTrsm,
    
    ) where

import Control.Monad ( when )
import Data.Maybe ( fromJust )
import System.IO.Unsafe
import Unsafe.Coerce

import Data.Matrix.Dense.Internal
import qualified Data.Matrix.Dense.Internal as M
import qualified Data.Matrix.Dense.Operations as M
import Data.Vector.Dense.Internal
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V

import BLAS.Access
import BLAS.Elem ( BLAS3 )
import qualified BLAS.Elem as E
import BLAS.C.Types ( cblasDiag, cblasUpLo, cblasTrans, colMajor, 
    noTrans, conjTrans, leftSide, rightSide )
import BLAS.Internal ( checkMatVecMult, checkMatMatMult, checkMatVecSolv,
    checkMatMatSolv )
import BLAS.Types ( Trans(..), flipTrans, flipUpLo )
                                   
import qualified BLAS.C as BLAS
import qualified BLAS.C.Types as BLAS

import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly
import BLAS.Matrix.Solve

import Data.Matrix.Tri

instance (BLAS3 e) => IMatrix (Tri (DMatrix Imm)) e where
    unsafeSApply k t x = unsafePerformIO $ unsafeGetSApply k t x
    {-# NOINLINE unsafeSApply #-}

    unsafeSApplyMat k t a = unsafePerformIO $ unsafeGetSApplyMat k t a
    {-# NOINLINE unsafeSApplyMat #-}

instance (BLAS3 e) => ISolve (Tri (DMatrix Imm)) e where
    unsafeSolve t x = unsafePerformIO $ unsafeGetSolve t x
    {-# NOINLINE unsafeSolve #-}

    unsafeSolveMat t a = unsafePerformIO $ unsafeGetSolveMat t a
    {-# NOINLINE unsafeSolveMat #-}

instance (BLAS3 e) => RMatrix (Tri (DMatrix s)) e where
    unsafeGetSApply k t x = do
        x' <- newCopy x
        unsafeTrmv (unsafeCoerce t) (V.unsafeThaw x')
        V.scaleBy k (V.unsafeThaw x')
        return (unsafeCoerce x')
    
    unsafeGetSApplyMat k t a = do
        a' <- newCopy a
        unsafeTrmm (unsafeCoerce t) (M.unsafeThaw a')
        M.scaleBy k (M.unsafeThaw a')
        return (unsafeCoerce a')

instance (BLAS3 e) => RSolve (Tri (DMatrix s)) e where
    unsafeGetSolve t x = do
        x' <- newCopy x
        unsafeTrsv (unsafeCoerce t) (V.unsafeThaw x')
        return (unsafeCoerce x')
    
    unsafeGetSolveMat t a = do
        a' <- newCopy a
        unsafeTrsm (unsafeCoerce t) (M.unsafeThaw a')
        return (unsafeCoerce a')

trmv :: (BLAS3 e) => Tri (DMatrix t) (n,n) e -> IOVector n e -> IO ()
trmv t x =
    checkMatVecMult (numRows t, numCols t) (dim x) $
        unsafeTrmv t x

unsafeTrmv :: (BLAS3 e) => Tri (DMatrix t) (n,n) e -> IOVector n e -> IO ()
unsafeTrmv _ x
    | dim x == 0 = return ()
unsafeTrmv t x
    | isConj x =
        let b = fromJust $ maybeFromCol x
        in trmm t b
unsafeTrmv t x =
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
               

trmm :: (BLAS3 e) => Tri (DMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
trmm t b = 
    checkMatMatMult (numRows t, numCols t) (shape b) $
        unsafeTrmm t b
               
unsafeTrmm :: (BLAS3 e) => Tri (DMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
unsafeTrmm _ b
    | M.numRows b == 0 || M.numCols b == 0 = return ()
unsafeTrmm t b =
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        (h,u')    = if isHerm a then (ConjTrans, flipUpLo u) else (NoTrans, u)
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if M.isHerm b
                        then (rightSide, flipTrans h, n, m, E.conj alpha)
                        else (leftSide , h          , m, n, alpha       )
        uploA     = cblasUpLo u'
        transA    = cblasTrans h'
        diagA     = cblasDiag d
        ldA       = ldaOf a
        ldB       = ldaOf b
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           M.unsafeWithElemPtr b (0,0) $ \pB ->
               BLAS.trmm order side uploA transA diagA m' n' alpha' pA ldA pB ldB

trsv :: (BLAS3 e) =>Tri (DMatrix t) (n,n) e -> IOVector n e -> IO ()
trsv t x =
    checkMatVecSolv (numRows t, numCols t) (dim x) $
        unsafeTrsv t x

unsafeTrsv :: (BLAS3 e) =>Tri (DMatrix t) (n,n) e -> IOVector n e -> IO ()
unsafeTrsv _ x
    | dim x == 0 = return ()
unsafeTrsv t x
    | isConj x =
        let b = fromJust $ maybeFromCol x
        in trsm t b
unsafeTrsv t x =
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

trsm :: (BLAS3 e) => Tri (DMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
trsm t b =
    checkMatMatSolv (numRows t, numCols t) (shape b) $
        unsafeTrsm t b

unsafeTrsm :: (BLAS3 e) => Tri (DMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
unsafeTrsm _ b
    | M.numRows b == 0 || M.numCols b == 0 = return ()
unsafeTrsm t b =
    let (u,d,alpha,a) = toBase t
        order     = colMajor
        (h,u')    = if isHerm a then (ConjTrans, flipUpLo u) else (NoTrans, u)
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if isHerm b
                        then (rightSide, flipTrans h, n, m, E.conj alpha)
                        else (leftSide , h          , m, n, alpha       )
        uploA     = cblasUpLo u'
        transA    = cblasTrans h'
        diagA     = cblasDiag d
        ldA       = ldaOf a
        ldB       = ldaOf b
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           M.unsafeWithElemPtr b (0,0) $ \pB -> do
               BLAS.trsm order side uploA transA diagA m' n' (1/alpha') pA ldA pB ldB
               