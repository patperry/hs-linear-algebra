{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Apply.Immutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Apply.Immutable (
    (<*>),
    (<**>),
    sapply,
    sapplyMat,
    
    IApply(..),
    unsafeApply,
    unsafeApplyMat,
    ) where

import BLAS.Elem( Elem, BLAS3 )
import BLAS.Internal ( checkMatVecMult, checkMatMatMult )
import BLAS.Matrix.Base
import BLAS.Matrix.Apply.Read

import Data.Vector.Dense
import Data.Vector.Dense.ST( runSTVector )
import Data.Matrix.Dense.Internal
import Data.Matrix.Dense.ST( runSTMatrix )

infixr 7 <*>, <**>

class (Elem e, BaseMatrix a e) => IApply a e where
    unsafeSApply :: e -> a (m,n) e -> Vector n e -> Vector m e
    unsafeSApplyMat :: e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e


-- | Apply to a vector
(<*>) :: (IApply a e) => a (m,n) e -> Vector n e -> Vector m e
(<*>) a x = checkMatVecMult (shape a) (dim x) $ unsafeApply a x
    
-- | Apply to a matrix
(<**>) :: (IApply a e) => a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
(<**>) a b = checkMatMatMult (shape a) (shape b) $ unsafeApplyMat a b

sapply :: (IApply a e) => e -> a (m,n) e -> Vector n e -> Vector m e
sapply k a x = checkMatVecMult (shape a) (dim x) $ unsafeSApply k a x
    
sapplyMat :: (IApply a e) => e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e    
sapplyMat k a b = checkMatMatMult (shape a) (shape b) $ unsafeSApplyMat k a b

unsafeApply :: (IApply a e) => a (m,n) e -> Vector n e -> Vector m e
unsafeApply = unsafeSApply 1

unsafeApplyMat :: (IApply a e) => a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
unsafeApplyMat = unsafeSApplyMat 1

instance (BLAS3 e) => IApply Matrix e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b

{-# RULES
"scale.apply/sapply"       forall k a x. (<*>) (k *> a) x = sapply k a x
"scale.applyMat/sapplyMat" forall k a b. (<**>) (k *> a) b = sapplyMat k a b
  #-}
