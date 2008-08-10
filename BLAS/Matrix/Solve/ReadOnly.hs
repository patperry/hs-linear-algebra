{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.Solve.ReadOnly
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Matrix.Solve.ReadOnly (
    RSolve(..),
    
    -- * Matrix and vector solving
    getSolve,
    getSolveMat,
    getSSolve,
    getSSolveMat,
    
    -- * In-place solving
    doSolve,
    doSolveMat,
    doSSolve,
    doSSolveMat,
    doSolve_,
    doSolveMat_,
    doSSolve_,
    doSSolveMat_,
    
    -- * unsafe operations
    unsafeGetSolve,
    unsafeGetSolveMat,
    unsafeGetSSolve,
    unsafeGetSSolveMat,

    ) where

import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv, checkMatVecSolvTo,
    checkMatMatSolvTo, checkSquare )
import BLAS.Matrix.Base ( numRows, numCols )
import BLAS.Matrix.ReadOnly

import Data.Vector.Dense.IO ( DVector, IOVector, dim, newVector_ )
import qualified Data.Vector.Dense.IO as V
import Data.Matrix.Dense.IO ( DMatrix, IOMatrix, shape, newMatrix_ )
import qualified Data.Matrix.Dense.IO as M

import Unsafe.Coerce


class RMatrix a e => RSolve a e where
    unsafeDoSolve :: a (m,n) e -> DVector t m e -> IOVector n e -> IO ()
    unsafeDoSolve = unsafeDoSSolve 1
    
    unsafeDoSolveMat :: a (m,n) e -> DMatrix t (m,k) e -> IOMatrix (n,k) e -> IO ()
    unsafeDoSolveMat = unsafeDoSSolveMat 1
    
    unsafeDoSSolve :: e -> a (m,n) e -> DVector t m e -> IOVector n e -> IO ()
    unsafeDoSSolve alpha a y x = do
        V.scaleBy alpha x
        unsafeDoSolve a y x
    
    unsafeDoSSolveMat :: e -> a (m,n) e -> DMatrix t (m,k) e -> IOMatrix (n,k) e -> IO ()
    unsafeDoSSolveMat alpha a c b = do
        M.scaleBy alpha b
        unsafeDoSolveMat a c b

    unsafeDoSolve_ :: a (n,n) e -> IOVector n e -> IO ()
    unsafeDoSolve_ = unsafeDoSSolve_ 1

    unsafeDoSSolve_ :: e -> a (n,n) e -> IOVector n e -> IO ()
    unsafeDoSSolve_ alpha a x = do
        V.scaleBy alpha x
        unsafeDoSolve_ a x
        
    unsafeDoSolveMat_ :: a (m,m) e -> IOMatrix (m,n) e -> IO ()
    unsafeDoSolveMat_ = unsafeDoSSolveMat_ 1
        
    unsafeDoSSolveMat_ :: e -> a (m,m) e -> IOMatrix (m,n) e -> IO ()
    unsafeDoSSolveMat_ alpha a b = do
        M.scaleBy alpha b
        unsafeDoSolveMat_ a b


unsafeGetSolve :: (RSolve a e) => 
    a (m,n) e -> DVector t m e -> IO (DVector r n e)
unsafeGetSolve a y = do
    x <- newVector_ (numCols a)
    unsafeDoSolve a y x
    return (unsafeCoerce x)
    
unsafeGetSSolve :: (RSolve a e) => 
    e -> a (m,n) e -> DVector t m e -> IO (DVector r n e)
unsafeGetSSolve alpha a y = do
    x <- newVector_ (numCols a)
    unsafeDoSSolve alpha a y x
    return (unsafeCoerce x)
    
unsafeGetSolveMat :: (RSolve a e) => 
    a (m,n) e -> DMatrix t (m,k) e -> IO (DMatrix r (n,k) e)
unsafeGetSolveMat a c = do
    b <- newMatrix_ (numCols a, numCols c)
    unsafeDoSolveMat a c b
    return (unsafeCoerce b)

unsafeGetSSolveMat :: (RSolve a e) => 
    e -> a (m,n) e -> DMatrix t (m,k) e -> IO (DMatrix r (n,k) e)
unsafeGetSSolveMat alpha a c = do
    b <- newMatrix_ (numCols a, numCols c)
    unsafeDoSSolveMat alpha a c b
    return (unsafeCoerce b)

-- | Solve for a vector
getSolve :: (RSolve a e) => 
    a (m,n) e -> DVector t m e -> IO (DVector r n e)
getSolve a y = 
    checkMatVecSolv (shape a) (dim y) $
        unsafeGetSolve a y

-- | Solve for a vector and scale
getSSolve :: (RSolve a e) => 
    e -> a (m,n) e -> DVector t m e -> IO (DVector r n e)
getSSolve alpha a y = 
    checkMatVecSolv (shape a) (dim y) $
        unsafeGetSSolve alpha a y

-- | Solve for a matrix
getSolveMat :: (RSolve a e) => 
    a (m,n) e -> DMatrix t (m,k) e -> IO (DMatrix r (n,k) e)
getSolveMat a b =
    checkMatMatSolv (shape a) (shape b) $
            unsafeGetSolveMat a b
            
-- | Solve for a matrix and scale
getSSolveMat :: (RSolve a e) => 
    e -> a (m,n) e -> DMatrix t (m,k) e -> IO (DMatrix r (n,k) e)
getSSolveMat alpha a b =
    checkMatMatSolv (shape a) (shape b) $
            unsafeGetSSolveMat alpha a b


doSolve :: (RSolve a e) => 
    a (m,n) e -> DVector t m e -> IOVector n e -> IO ()
doSolve a y x =
    checkMatVecSolvTo (shape a) (dim y) (dim x) $
        unsafeDoSolve a y x
        
doSolveMat :: (RSolve a e) => 
    a (m,n) e -> DMatrix t (m,k) e -> IOMatrix (n,k) e -> IO ()
doSolveMat a c b =
    checkMatMatSolvTo (shape a) (shape c) (shape b) $
        unsafeDoSolveMat a c b
    
doSSolve :: (RSolve a e) => 
    e -> a (m,n) e -> DVector t m e -> IOVector n e -> IO ()
doSSolve alpha a y x =
    checkMatVecSolvTo (shape a) (dim y) (dim x) $
        unsafeDoSSolve alpha a y x

doSSolveMat :: (RSolve a e) => 
    e -> a (m,n) e -> DMatrix t (m,k) e -> IOMatrix (n,k) e -> IO ()
doSSolveMat alpha a c b =
    checkMatMatSolvTo (shape a) (shape c) (shape b) $
        unsafeDoSSolveMat alpha a c b

doSolve_ :: (RSolve a e) => 
    a (n,n) e -> IOVector n e -> IO ()
doSolve_ a x =
    checkSquare (shape a) $
        checkMatVecSolv (shape a) (dim x) $
            unsafeDoSolve_ a x

doSSolve_ :: (RSolve a e) => 
    e -> a (n,n) e -> IOVector n e -> IO ()
doSSolve_ alpha a x =
    checkSquare (shape a) $
        checkMatVecSolv (shape a) (dim x) $
            unsafeDoSSolve_ alpha a x

doSolveMat_ :: (RSolve a e) => 
    a (m,m) e -> IOMatrix (m,n) e -> IO ()
doSolveMat_ a b =
    checkSquare (shape a) $
        checkMatMatSolv (shape a) (shape b) $
            unsafeDoSolveMat_ a b

doSSolveMat_ :: (RSolve a e) => 
    e -> a (m,m) e -> IOMatrix (m,n) e -> IO ()
doSSolveMat_ alpha a b =
    checkSquare (shape a) $
        checkMatMatSolv (shape a) (shape b) $
            unsafeDoSSolveMat_ alpha a b

