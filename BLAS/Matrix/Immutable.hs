{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Immutable
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Immutable (
    IMatrix(..)
    ) where

import BLAS.Access
import BLAS.Elem ( Elem, BLAS3 )
import BLAS.Internal ( checkMatVecMult, checkMatMatMult )
import qualified BLAS.Matrix.Base as Base
import Data.Vector.Dense
import Data.Matrix.Dense.Internal
import qualified Data.Matrix.Dense.Operations as M

infixr 7 <*>, <**>

class (Base.Matrix a, Elem e) => IMatrix a e where
    -- | Apply to a vector
    (<*>) :: a (m,n) e -> Vector n e -> Vector m e
    (<*>) a x = 
        checkMatVecMult (numRows a, numCols a) (dim x) $ unsafeApply a x
        
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
    
    unsafeSApplyMat :: e -> a (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e    


instance (BLAS3 e) => IMatrix (DMatrix Imm) e where
    unsafeSApply    = M.unsafeSApply
    unsafeSApplyMat = M.unsafeSApplyMat

{-# RULES
"scale.apply/sapply"       forall k a x. (<*>) (k *> a) x = sapply k a x
"scale.applyMat/sapplyMat" forall k a b. (<**>) (k *> a) b = sapplyMat k a b
  #-}
