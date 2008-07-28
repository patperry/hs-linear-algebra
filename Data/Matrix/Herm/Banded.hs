{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Herm.Banded
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Herm.Banded (
    module Data.Matrix.Herm,
    module BLAS.Matrix.Immutable,
    module BLAS.Matrix.ReadOnly,

    hbmv,
    hbmm,
    ) where

import Control.Monad ( zipWithM_ )
import System.IO.Unsafe
import Unsafe.Coerce

import BLAS.Access
import BLAS.Elem ( BLAS2 )
import BLAS.C ( colMajor, cblasUpLo )
import BLAS.Types ( flipUpLo )
import qualified BLAS.Elem as E
import qualified BLAS.C as BLAS

import Data.Matrix.Banded.Internal
import Data.Matrix.Dense.Internal ( DMatrix, IOMatrix )
import Data.Vector.Dense.Internal
import qualified Data.Matrix.Banded.Internal as B
import qualified Data.Matrix.Dense.Internal as M
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V

import Data.Matrix.Herm
import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly


instance (BLAS2 e) => IMatrix (Herm (BMatrix Imm)) e where
    (<*>) h x = unsafePerformIO $ getApply h x
    {-# NOINLINE (<*>) #-}

    (<**>) h a = unsafePerformIO $ getApplyMat h a
    {-# NOINLINE (<**>) #-}


instance (BLAS2 e) => RMatrix (Herm (BMatrix s)) e where
    getApply h x = do
        y <- newZero (dim x)
        hbmv 1 (unsafeCoerce h) x 1 y
        return (unsafeCoerce y)
    
    getApplyMat h a = do
        b <- newZero (shape a)
        hbmm 1 (unsafeCoerce h) a 1 b
        return (unsafeCoerce b)


hbmv :: (BLAS2 e) => e -> Herm (BMatrix t) (n,n) e -> DVector s n e -> e -> IOVector n e -> IO ()
hbmv alpha h x beta y
    | numRows h == 0 =
        return ()
    | isConj y = do
        V.doConj y
        hbmv alpha h x beta (V.conj y)
        V.doConj y
    | isConj x = do
        x' <- newCopy x
        V.doConj (V.unsafeThaw x')
        hbmv alpha h (conj x') beta y
    | otherwise =
        let order   = colMajor
            (u,e,a) = toBase h
            alpha'  = alpha * e
            n       = numCols a
            k       = case u of 
                          Upper -> numUpper a
                          Lower -> numLower a      
            (u',alpha'') 
                    = case (isHerm a) of
                          True  -> (flipUpLo u, E.conj alpha')
                          False -> (u, alpha')
            uploA   = cblasUpLo u'
            ldA     = ldaOf a
            incX    = strideOf x
            incY    = strideOf y
            withPtrA 
                    = case u' of Upper -> B.unsafeWithBasePtr a
                                 Lower -> B.unsafeWithElemPtr a (0,0)
                    
        in withPtrA  $ \pA ->
               V.unsafeWithElemPtr x 0 $ \pX ->
                    V.unsafeWithElemPtr y 0 $ \pY -> do
                        BLAS.hbmv order uploA n k alpha'' pA ldA pX incX beta pY incY

hbmm :: (BLAS2 e) => e -> Herm (BMatrix t) (m,m) e -> DMatrix s (m,n) e -> e -> IOMatrix (m,n) e -> IO ()
hbmm alpha h b beta c =
    zipWithM_ (\x y -> hbmv alpha h x beta y) (M.cols b) (M.cols c)
