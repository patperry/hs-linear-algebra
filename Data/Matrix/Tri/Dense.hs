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
    
    ) where

import Control.Monad ( when )
import Data.Maybe ( fromJust )

import Data.Matrix.Dense.Internal hiding ( diag )
import qualified Data.Matrix.Dense.Internal as M
import qualified Data.Matrix.Dense.Operations as M
import Data.Vector.Dense.Internal
import Data.Vector.Dense.Operations( unsafeCopyVector )
import Data.Matrix.Dense.Operations( unsafeCopyMatrix )
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V

import BLAS.Access

import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly
import BLAS.Matrix.Solve

import Data.Matrix.Tri

instance (BLAS3 e) => IMatrix (Tri (DMatrix Imm)) e where
instance (BLAS3 e) => ISolve (Tri (DMatrix Imm)) e where



instance (BLAS3 e) => RMatrix (Tri (DMatrix s)) e where
    unsafeDoSApply_    = trmv
    unsafeDoSApplyMat_ = trmm
    


        
    
instance (BLAS3 e) => RSolve (Tri (DMatrix s)) e where
    unsafeDoSSolve_    = trsv
    unsafeDoSSolveMat_ = trsm

    unsafeDoSSolve alpha t y x =
        case (u, toLower d a, toUpper d a) of
            (Lower,Left t',_) -> do
                unsafeCopyVector x (coerceVector y)
                trsv alpha t' (coerceVector x)
                
            (Lower,Right (t',_),_) -> do
                let y1 = unsafeSubvector y 0            (numRows t')
                unsafeCopyVector x y1
                trsv alpha t' x
                
            (Upper,_,Left t') -> do
                unsafeCopyVector x (coerceVector y)
                trsv alpha t' x

            (Upper,_,Right (t',r)) ->
                let x1 = unsafeSubvector x 0            (numCols t')
                    x2 = unsafeSubvector x (numCols t') (numCols r)
                in do
                    unsafeCopyVector x1 y
                    trsv alpha t' x1
                    setZero x2
      where
        (u,d,a) = toBase t


    unsafeDoSSolveMat alpha t c b =
        case (u, toLower d a, toUpper d a) of
            (Lower,Left t',_) -> do
                unsafeCopyMatrix b (coerceMatrix c)
                trsm alpha t' (coerceMatrix b)
                
            (Lower,Right (t',_),_) -> do
                let c1 = unsafeSubmatrix c (0,0)          (numRows t',numCols c)
                unsafeCopyMatrix b c1
                trsm alpha t' b
                
            (Upper,_,Left t') -> do
                unsafeCopyMatrix (coerceMatrix b) c
                trsm alpha t' (coerceMatrix b)

            (Upper,_,Right (t',r)) ->
                let b1 = unsafeSubmatrix b (0,0)          (numCols t',numCols b)
                    b2 = unsafeSubmatrix b (numCols t',0) (numCols r ,numCols b)
                in do
                    unsafeCopyMatrix b1 c
                    trsm alpha t' b1
                    setZero b2
      where
        (u,d,a) = toBase t


               

trsv :: (BLAS3 e) => e -> Tri (DMatrix t) (n,n) e -> IOVector n e -> IO ()
trsv _ _ x
    | dim x == 0 = return ()
trsv alpha t x
    | isConj x =
        let b = fromJust $ maybeFromCol x
        in trsm alpha t b
trsv alpha t x =
    let (u,d,a) = toBase t
        order     = colMajor
        (transA,u') = if isHerm a then (conjTrans, flipUpLo u) else (noTrans, u)
        uploA     = cblasUpLo u'
        diagA     = cblasDiag d
        n         = dim x
        ldA       = ldaOf a
        incX      = strideOf x
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           V.unsafeWithElemPtr x 0 $ \pX -> do
               when (alpha /= 1) $ V.scaleBy alpha x
               BLAS.trsv order uploA transA diagA n pA ldA pX incX


trsm :: (BLAS3 e) => e -> Tri (DMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
trsm _ _ b
    | M.numRows b == 0 || M.numCols b == 0 = return ()
trsm alpha t b =
    let (u,d,a)   = toBase t
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
               BLAS.trsm order side uploA transA diagA m' n' alpha' pA ldA pB ldB

               