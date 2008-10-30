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
    ) where


import Data.Matrix.Banded.Internal
import qualified Data.Matrix.Banded.Internal as B
import Data.Matrix.Dense.IO ( IOMatrix, coerceMatrix, unsafeCopyMatrix )
import qualified Data.Matrix.Dense.IO as M
import Data.Vector.Dense.IO ( IOVector, isConj, strideOf, conj, coerceVector, 
    unsafeCopyVector )
import qualified Data.Vector.Dense.IO as V

import BLAS.Access
import BLAS.Types ( flipUpLo )

import BLAS.C ( BLAS2, cblasDiag, cblasUpLo, colMajor, noTrans, conjTrans )                                   
import qualified BLAS.C as BLAS

import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly
import BLAS.Matrix.Solve

import Data.Matrix.Tri

instance (BLAS2 e) => IMatrix (Tri (BMatrix Imm)) e where

instance (BLAS2 e) => RMatrix (Tri (BMatrix s)) e where

    unsafeDoSApply alpha a x y = do
        unsafeCopyVector (coerceVector y) x
        unsafeDoSApply_ alpha (coerceTri a) y
    
    unsafeDoSApplyMat alpha a b c = do
        unsafeCopyMatrix (coerceMatrix c) b
        unsafeDoSApplyMat_ alpha (coerceTri a) c
        
    unsafeDoSApply_    = tbmv
    unsafeDoSApplyMat_ = tbmm

instance (BLAS2 e) => ISolve (Tri (BMatrix Imm)) e where

instance (BLAS2 e) => RSolve (Tri (BMatrix s)) e where
    unsafeDoSSolve alpha a y x = do
        unsafeCopyVector (coerceVector x) y
        unsafeDoSSolve_ alpha (coerceTri a) x
    
    unsafeDoSSolveMat alpha a c b = do
        unsafeCopyMatrix (coerceMatrix b) c
        unsafeDoSSolveMat_ alpha (coerceTri a) b

    unsafeDoSSolve_    = tbsv
    unsafeDoSSolveMat_ = tbsm


tbmv :: (BLAS2 e) => e -> Tri (BMatrix t) (n,n) e -> IOVector n e -> IO ()
              
               
tbmm :: (BLAS2 e) => e -> Tri (BMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()


tbsv :: (BLAS2 e) => e -> Tri (BMatrix t) (n,n) e -> IOVector n e -> IO ()
tbsv alpha t x | isConj x = do
    V.doConj x
    tbsv alpha t (conj x)
    V.doConj x
    
tbsv alpha t x = 
    let (u,d,a) = toBase t
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
               V.scaleBy alpha x
               BLAS.tbsv order uploA transA diagA n k pA ldA pX incX


tbsm :: (BLAS2 e) => e -> Tri (BMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
tbsm 1     t b = mapM_ (\x -> tbsv 1 t x) (M.cols b)
tbsm alpha t b = M.scaleBy alpha b >> tbsm 1 t b
