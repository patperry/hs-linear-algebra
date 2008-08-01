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
    tbsm,
    
    unsafeTbmv,
    unsafeTbsv,
    unsafeTbmm,
    unsafeTbsm,
    
    ) where

import Control.Monad ( when )
import System.IO.Unsafe
import Unsafe.Coerce

import Data.Matrix.Banded.Internal
import qualified Data.Matrix.Banded.Internal as B
import Data.Matrix.Dense.Internal ( IOMatrix )
import qualified Data.Matrix.Dense.Internal as M
import qualified Data.Matrix.Dense.Operations as M
import Data.Vector.Dense.Internal
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V

import BLAS.Access
import BLAS.Types ( flipUpLo )

import BLAS.C ( BLAS2, cblasDiag, cblasUpLo, colMajor, noTrans, conjTrans )                                   
import qualified BLAS.C as BLAS
import BLAS.Internal ( checkMatVecMult, checkMatMatMult, checkMatVecSolv,
    checkMatMatSolv )

import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly
import BLAS.Matrix.Solve

import Data.Matrix.Tri

instance (BLAS2 e) => IMatrix (Tri (BMatrix Imm)) e where
    unsafeSApply k t x = unsafePerformIO $ unsafeGetSApply k t x
    {-# NOINLINE unsafeSApply #-}

    unsafeSApplyMat k t a = unsafePerformIO $ unsafeGetSApplyMat k t a
    {-# NOINLINE unsafeSApplyMat #-}

instance (BLAS2 e) => ISolve (Tri (BMatrix Imm)) e where
    unsafeSolve t x = unsafePerformIO $ unsafeGetSolve t x
    {-# NOINLINE unsafeSolve #-}

    unsafeSolveMat t a = unsafePerformIO $ unsafeGetSolveMat t a
    {-# NOINLINE unsafeSolveMat #-}

instance (BLAS2 e) => RMatrix (Tri (BMatrix s)) e where
    unsafeGetSApply k t x = do
        x' <- newCopy x
        unsafeTbmv (unsafeCoerce t) (V.unsafeThaw x')
        V.scaleBy k (V.unsafeThaw x')
        return (unsafeCoerce x')
    
    unsafeGetSApplyMat k t a = do
        a' <- newCopy a
        unsafeTbmm (unsafeCoerce t) (M.unsafeThaw a')
        M.scaleBy k (M.unsafeThaw a')
        return (unsafeCoerce a')

instance (BLAS2 e) => RSolve (Tri (BMatrix s)) e where
    unsafeGetSolve t x = do
        x' <- newCopy x
        unsafeTbsv (unsafeCoerce t) (V.unsafeThaw x')
        return (unsafeCoerce x')
    
    unsafeGetSolveMat t a = do
        a' <- newCopy a
        unsafeTbsm (unsafeCoerce t) (M.unsafeThaw a')
        return (unsafeCoerce a')




tbmv :: (BLAS2 e) => Tri (BMatrix t) (n,n) e -> IOVector n e -> IO ()
tbmv t x =
    checkMatVecMult (numRows t, numCols t) (dim x) $
        unsafeTbmv t x

unsafeTbmv :: (BLAS2 e) => Tri (BMatrix t) (n,n) e -> IOVector n e -> IO ()
unsafeTbmv t x | isConj x = do
    V.doConj x
    tbmv t (conj x)
    V.doConj x

unsafeTbmv t x =
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
tbmm t b =
    checkMatMatMult (numRows t, numCols t) (shape b) $
        unsafeTbmm t b

unsafeTbmm :: (BLAS2 e) => Tri (BMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
unsafeTbmm t b = mapM_ (\x -> unsafeTbmv t x) (M.cols b)

tbsv :: (BLAS2 e) => Tri (BMatrix t) (n,n) e -> IOVector n e -> IO ()
tbsv t x =
    checkMatVecSolv (numRows t, numCols t) (dim x) $
        unsafeTbsv t x

unsafeTbsv :: (BLAS2 e) => Tri (BMatrix t) (n,n) e -> IOVector n e -> IO ()
unsafeTbsv t x | isConj x = do
    V.doConj x
    tbsv t (conj x)
    V.doConj x
    
unsafeTbsv t x = 
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
tbsm t b =
    checkMatMatSolv (numRows t, numCols t) (shape b) $
        unsafeTbsm t b

unsafeTbsm :: (BLAS2 e) => Tri (BMatrix t) (m,m) e -> IOMatrix (m,n) e -> IO ()
unsafeTbsm t b = mapM_ (\x -> unsafeTbsv t x) (M.cols b)
