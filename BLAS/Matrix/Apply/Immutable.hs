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
    IApply(..)
    ) where

import BLAS.Elem ( BLAS3 )
import BLAS.Internal ( checkMatVecMult, checkMatMatMult )
import BLAS.Matrix.Base

import Data.Vector.Dense
import Data.Matrix.Dense

import System.IO.Unsafe ( unsafePerformIO )

infixr 7 <*>, <**>

class (BaseMatrix a e) => IApply a e where
    -- | Apply to a vector
    (<*>) :: a (m,n) e -> Vector n e -> Vector m e
    (<*>) a x = 
        checkMatVecMult (shape a) (dim x) $ unsafeApply a x
        
    -- | Apply to a matrix
    (<**>) :: a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
    (<**>) a b = 
        checkMatMatMult (numRows a, numCols a) 
                        (numRows b, numCols b) $ unsafeApplyMat a b

    sapply :: e -> a (m,n) e -> Vector n e -> Vector m e
    sapply k a x =
        checkMatVecMult (numRows a, numCols a) (dim x) $ unsafeSApply k a x
        
    sapplyMat :: e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e    
    sapplyMat k a b = 
        checkMatMatMult (numRows a, numCols a) 
                        (numRows b, numCols b) $ unsafeSApplyMat k a b


    unsafeApply :: a (m,n) e -> Vector n e -> Vector m e
    unsafeApply = unsafeSApply 1
    
    unsafeApplyMat :: a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
    unsafeApplyMat = unsafeSApplyMat 1
    
    unsafeSApply :: e -> a (m,n) e -> Vector n e -> Vector m e
    unsafeSApply alpha a x = unsafePerformIO $ unsafeGetSApply alpha a x
    {-# NOINLINE unsafeSApply #-}

    unsafeSApplyMat :: e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
    unsafeSApplyMat alpha a b = unsafePerformIO $ unsafeGetSApplyMat alpha a b
    {-# NOINLINE unsafeSApplyMat #-}


-- instance (BLAS3 e) => IApply Matrix e where

{-# RULES
"scale.apply/sapply"       forall k a x. (<*>) (k *> a) x = sapply k a x
"scale.applyMat/sapplyMat" forall k a b. (<**>) (k *> a) b = sapplyMat k a b
  #-}
