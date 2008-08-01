{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.ReadOnly
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.ReadOnly (
    RMatrix(..)
    ) where

import BLAS.Elem ( Elem, BLAS3 )
import BLAS.Internal ( checkMatVecMult, checkMatMatMult )
import qualified BLAS.Matrix.Base as Base
import Data.Vector.Dense.Internal
import Data.Matrix.Dense.Internal
import qualified Data.Matrix.Dense.Operations as M

class (Base.Matrix a, Elem e) => RMatrix a e where
    -- | Apply to a vector
    getApply :: a (m,n) e -> DVector t n e -> IO (DVector r m e)
    getApply a x =
        checkMatVecMult (numRows a, numCols a) (dim x) $ unsafeGetApply a x
    
    -- | Apply to a matrix
    getApplyMat :: a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
    getApplyMat a b =
        checkMatMatMult (numRows a, numCols a) 
                        (numRows b, numCols b) $ unsafeGetApplyMat a b

    -- | Scale and apply to a vector
    getSApply :: e -> a (m,n) e -> DVector t n e -> IO (DVector r m e)
    getSApply k a x =
        checkMatVecMult (numRows a, numCols a) (dim x) $ unsafeGetSApply k a x
    
    -- | Scale and apply to a matrix
    getSApplyMat :: e -> a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
    getSApplyMat k a b =
        checkMatMatMult (numRows a, numCols a) 
                        (numRows b, numCols b) $ unsafeGetSApplyMat k a b
    
    unsafeGetApply :: a (m,n) e -> DVector t n e -> IO (DVector r m e)
    unsafeGetApply = unsafeGetSApply 1
    
    unsafeGetApplyMat :: a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
    unsafeGetApplyMat = unsafeGetSApplyMat 1
    
    unsafeGetSApply :: e -> a (m,n) e -> DVector t n e -> IO (DVector r m e)
    
    unsafeGetSApplyMat :: e -> a (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
    

instance (BLAS3 e) => RMatrix (DMatrix t) e where
    unsafeGetSApply    = M.unsafeGetSApply
    unsafeGetSApplyMat = M.unsafeGetSApplyMat
