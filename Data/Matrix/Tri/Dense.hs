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
import System.IO.Unsafe
import Unsafe.Coerce

import Data.Matrix.Dense.Internal hiding ( diag )
import qualified Data.Matrix.Dense.Internal as M
import qualified Data.Matrix.Dense.Operations as M
import Data.Vector.Dense.Internal
import Data.Vector.Dense.Operations( unsafeCopyVector )
import Data.Matrix.Dense.Operations( unsafeCopyMatrix )
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V

import BLAS.Access
import BLAS.Elem ( Elem, BLAS3 )
import qualified BLAS.Elem as E
import BLAS.C.Types ( cblasDiag, cblasUpLo, cblasTrans, colMajor, 
    noTrans, conjTrans, leftSide, rightSide )
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

toLower :: (Elem e) => Diag -> DMatrix s (m,n) e 
        -> Either (Tri (DMatrix s) (m,m) e) 
                  (Tri (DMatrix s) (n,n) e, DMatrix s (d,n) e)
toLower diag a =
    if m <= n
        then Left $  fromBase Lower diag (unsafeSubmatrix a (0,0) (m,m))
        else let t = fromBase Lower diag (unsafeSubmatrix a (0,0) (n,n))
                 r = unsafeSubmatrix a (n,0) (d,n)
             in Right (t,r)
  where
    (m,n) = shape a
    d     = m - n
    
toUpper :: (Elem e) => Diag -> DMatrix s (m,n) e
        -> Either (Tri (DMatrix s) (n,n) e)
                  (Tri (DMatrix s) (m,m) e, DMatrix s (m,d) e)
toUpper diag a =
    if n <= m
        then Left $  fromBase Upper diag (unsafeSubmatrix a (0,0) (n,n))
        else let t = fromBase Upper diag (unsafeSubmatrix a (0,0) (m,m))
                 r = unsafeSubmatrix a (0,m) (m,d)
             in Right (t,r)
  where
    (m,n) = shape a
    d     = n - m


instance (BLAS3 e) => RMatrix (Tri (DMatrix s)) e where
    unsafeDoSApply alpha t x y =
        case (u, toLower d a, toUpper d a) of
            (Lower,Left t',_) -> do
                unsafeCopyVector y (coerceVector x)
                trmv alpha t' y
                
            (Lower,Right (t',r),_) -> do
                let y1 = unsafeSubvector y 0            (numRows t')
                    y2 = unsafeSubvector y (numRows t') (numRows r)
                unsafeCopyVector y1 x
                trmv alpha t' y1
                unsafeDoSApply alpha r x y2
                
            (Upper,_,Left t') -> do
                unsafeCopyVector (coerceVector y) x
                trmv alpha t' (coerceVector y)

            (Upper,_,Right (t',r)) ->
                let x1 = unsafeSubvector x 0            (numCols t')
                    x2 = unsafeSubvector x (numCols t') (numCols r)
                in do
                    unsafeCopyVector y x1
                    trmv alpha t' y
                    unsafeDoSApplyAdd alpha r x2 1 y
      where
        (u,d,a) = toBase t


    unsafeDoSApplyMat alpha t b c =
        case (u, toLower d a, toUpper d a) of
            (Lower,Left t',_) -> do
                unsafeCopyMatrix c (coerceMatrix b)
                trmm alpha t' c
                
            (Lower,Right (t',r),_) -> do
                let c1 = unsafeSubmatrix c (0,0)          (numRows t',numCols c)
                    c2 = unsafeSubmatrix c (numRows t',0) (numRows r ,numCols c)
                unsafeCopyMatrix c1 b
                trmm alpha t' c1
                unsafeDoSApplyMat alpha r b c2
                
            (Upper,_,Left t') -> do
                unsafeCopyMatrix (coerceMatrix c) b
                trmm alpha t' (coerceMatrix c)

            (Upper,_,Right (t',r)) ->
                let b1 = unsafeSubmatrix b (0,0)          (numCols t',numCols b)
                    b2 = unsafeSubmatrix b (numCols t',0) (numCols r ,numCols b)
                in do
                    unsafeCopyMatrix c b1
                    trmm alpha t' c
                    unsafeDoSApplyAddMat alpha r b2 1 c
      where
        (u,d,a) = toBase t
        
    

instance (BLAS3 e) => RSolve (Tri (DMatrix s)) e where
    unsafeGetSolve t x = do
        x' <- newCopy x
        trsv 1 (unsafeCoerce t) (V.unsafeThaw x')
        return (unsafeCoerce x')
    
    unsafeGetSolveMat t a = do
        a' <- newCopy a
        trsm 1 (unsafeCoerce t) (M.unsafeThaw a')
        return (unsafeCoerce a')

trmv :: (BLAS3 e) => e -> Tri (DMatrix t) (n,n) e -> IOVector n e -> IO ()
trmv alpha t x 
    | dim x == 0 = 
        return ()
    | isConj x =
        let b = fromJust $ maybeFromCol x
        in trmm alpha t b
    | otherwise =
        let (u,d,a)   = toBase t
            order     = colMajor
            (transA,u') = if isHerm a then (conjTrans, flipUpLo u) else (noTrans, u)
            uploA     = cblasUpLo u'
            diagA     = cblasDiag d
            n         = dim x
            ldA       = ldaOf a
            incX      = strideOf x
        in M.unsafeWithElemPtr a (0,0) $ \pA ->
               V.unsafeWithElemPtr x 0 $ \pX -> do
                   BLAS.trmv order uploA transA diagA n pA ldA pX incX
                   when (alpha /= 1) $ V.scaleBy alpha x

               
trmm :: (BLAS3 e) => e -> Tri (DMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
trmm _ _ b
    | M.numRows b == 0 || M.numCols b == 0 = return ()
trmm alpha t b =
    let (u,d,a)   = toBase t
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
        n         = numCols a
        ldA       = ldaOf a
        incX      = strideOf x
    in M.unsafeWithElemPtr a (0,0) $ \pA ->
           V.unsafeWithElemPtr x 0 $ \pX -> do
               BLAS.trsv order uploA transA diagA n pA ldA pX incX
               when (alpha /= 1) $ V.scaleBy alpha x



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
               