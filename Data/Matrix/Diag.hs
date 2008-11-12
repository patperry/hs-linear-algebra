{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances, ExistentialQuantification #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Diag
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Diag (
    -- * The diagonal matrix types
    Diag,
    
    -- * Converting to and from @Vector@s
    diagFromVector,
    vectorFromDiag,
    
    -- * Coercing shapes
    coerceDiag,
    
    module BLAS.Matrix,
    
    ) where

import Control.Monad( zipWithM_ )
import Control.Monad.ST( ST )

import BLAS.Elem( BLAS1 )
import BLAS.Matrix hiding ( BaseMatrix )
import qualified BLAS.Matrix as BLAS
import BLAS.Tensor
import BLAS.UnsafeIOToM
import Unsafe.Coerce

import Data.AEq

import Data.Matrix.Dense.Class( ReadMatrix, WriteMatrix, unsafeCopyMatrix, 
    rowViews, coerceMatrix )
import Data.Matrix.Dense.ST( runSTMatrix )
import Data.Vector.Dense( Vector )
import Data.Vector.Dense.IO( IOVector )
import Data.Vector.Dense.ST( STVector, runSTVector )
import Data.Vector.Dense.Class( BaseVector, ReadVector, WriteVector,
    conj, dim, coerceVector, scaleBy, unsafeCopyVector, unsafeAxpyVector,
    unsafeMulVector, unsafeDivVector, newCopyVector )

data Diag x nn e = forall n. Diag (x n e)

coerceDiag :: Diag x mn e -> Diag x mn' e
coerceDiag = unsafeCoerce

diagFromVector :: (BaseVector x) => x n e -> Diag x (n,n) e
diagFromVector = Diag . coerceVector

vectorFromDiag :: (BaseVector x) => Diag x (n,n) e -> x n e
vectorFromDiag (Diag x) = coerceVector x

instance (BaseVector x) => BaseTensor (Diag x) (Int,Int) where
    shape  (Diag x) = (n,n) where n = dim x
    bounds (Diag x) = ((0,0),(n-1,n-1)) where n = dim x

instance (BaseVector x) => BLAS.BaseMatrix (Diag x) where
    herm (Diag x) = Diag (conj x)
    
        
instance (BLAS1 e) => IMatrix (Diag Vector) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b

instance (BLAS1 e) => MMatrix (Diag IOVector) e IO where
    unsafeDoSApplyAdd    = unsafeDoSApplyAddDiagVector
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatDiagVector
    unsafeDoSApply_      = unsafeDoSApplyDiagVector_
    unsafeDoSApplyMat_   = unsafeDoSApplyMatDiagVector_

instance (BLAS1 e) => MMatrix (Diag (STVector s)) e (ST s) where
    unsafeDoSApplyAdd    = unsafeDoSApplyAddDiagVector
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatDiagVector
    unsafeDoSApply_      = unsafeDoSApplyDiagVector_
    unsafeDoSApplyMat_   = unsafeDoSApplyMatDiagVector_

instance (BLAS1 e, UnsafeIOToM m) => MMatrix (Diag Vector) e m where
    unsafeDoSApplyAdd    = unsafeDoSApplyAddDiagVector
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatDiagVector
    unsafeDoSApply_      = unsafeDoSApplyDiagVector_
    unsafeDoSApplyMat_   = unsafeDoSApplyMatDiagVector_


unsafeDoSApplyAddDiagVector :: (ReadVector z m, ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> Diag z (r,s) e -> x s e -> e -> y r e -> m ()
unsafeDoSApplyAddDiagVector alpha a x beta y = do
    x' <- newCopyVector x
    unsafeDoSApplyDiagVector_ 1 (coerceDiag a) x'
    scaleBy beta y
    unsafeAxpyVector alpha x' (coerceVector y)

unsafeDoSApplyAddMatDiagVector :: (ReadVector x m, ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    e -> Diag x (r,s) e -> b (s,t) e ->  e -> c (r,t) e -> m ()
unsafeDoSApplyAddMatDiagVector alpha a b beta c = do
    scaleBy beta c
    ks <- getElems (vectorFromDiag $ coerceDiag a)
    let (kxs) = zip ks (rowViews b)
        ys    = rowViews c
    zipWithM_ (\(k,x) y -> unsafeAxpyVector (alpha*k) x y) kxs ys

unsafeDoSApplyDiagVector_ :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> Diag x (s,s) e -> y s e -> m ()
unsafeDoSApplyDiagVector_ alpha a x = do
    unsafeMulVector x (vectorFromDiag a)
    scaleBy alpha x

unsafeDoSApplyMatDiagVector_ :: (ReadVector x m, WriteMatrix b m, BLAS1 e) =>
    e -> Diag x (s,s) e -> b (s,t) e ->  m ()
unsafeDoSApplyMatDiagVector_ alpha a b = do
    ks <- getElems (vectorFromDiag a)
    zipWithM_ (\k r -> scaleBy (alpha*k) r) ks (rowViews b)


instance (BLAS1 e) => ISolve (Diag Vector) e where
    unsafeSSolve alpha a y    = runSTVector $ unsafeGetSSolve    alpha a y
    unsafeSSolveMat alpha a c = runSTMatrix $ unsafeGetSSolveMat alpha a c

instance (BLAS1 e) => MSolve (Diag IOVector) e IO where
    unsafeDoSSolve     = unsafeDoSSolveDiagVector
    unsafeDoSSolveMat  = unsafeDoSSolveMatDiagVector
    unsafeDoSSolve_    = unsafeDoSSolveDiagVector_
    unsafeDoSSolveMat_ = unsafeDoSSolveMatDiagVector_

instance (BLAS1 e) => MSolve (Diag (STVector s)) e (ST s) where
    unsafeDoSSolve     = unsafeDoSSolveDiagVector
    unsafeDoSSolveMat  = unsafeDoSSolveMatDiagVector
    unsafeDoSSolve_    = unsafeDoSSolveDiagVector_
    unsafeDoSSolveMat_ = unsafeDoSSolveMatDiagVector_

instance (BLAS1 e, UnsafeIOToM m) => MSolve (Diag Vector) e m where
    unsafeDoSSolve     = unsafeDoSSolveDiagVector
    unsafeDoSSolveMat  = unsafeDoSSolveMatDiagVector
    unsafeDoSSolve_    = unsafeDoSSolveDiagVector_
    unsafeDoSSolveMat_ = unsafeDoSSolveMatDiagVector_



unsafeDoSSolveDiagVector :: (ReadVector z m, ReadVector y m, WriteVector x m, BLAS1 e) =>
    e -> Diag z (r,s) e -> y r e -> x s e -> m ()
unsafeDoSSolveDiagVector alpha a y x = do
    unsafeCopyVector x (coerceVector y)
    unsafeDoSSolveDiagVector_ alpha (coerceDiag a) x

unsafeDoSSolveMatDiagVector :: (ReadVector x m, ReadMatrix c m, WriteMatrix b m, BLAS1 e) =>
    e -> Diag x (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
unsafeDoSSolveMatDiagVector alpha a c b = do
    unsafeCopyMatrix b (coerceMatrix c)
    unsafeDoSSolveMatDiagVector_ alpha (coerceDiag a) b

unsafeDoSSolveDiagVector_ :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> Diag x (k,k) e -> y k e -> m ()
unsafeDoSSolveDiagVector_ alpha a x = do
    scaleBy alpha x
    unsafeDivVector x (vectorFromDiag a)

unsafeDoSSolveMatDiagVector_ :: (ReadVector x m, WriteMatrix a m, BLAS1 e) =>
    e -> Diag x (k,k) e -> a (k,l) e -> m ()
unsafeDoSSolveMatDiagVector_ alpha a b = do
    scaleBy alpha b
    ks <- unsafeInterleaveM $ getElems (vectorFromDiag a)
    zipWithM_ (\k r -> scaleBy (1/k) r) ks (rowViews b)


instance (Show e, BLAS1 e) => Show (Diag Vector (n,n) e) where
    show x = "diagFromVector (" ++ show (vectorFromDiag x) ++ ")"

instance (Eq e, BLAS1 e) => Eq (Diag Vector (n,n) e) where
    (==) x y = (==) (vectorFromDiag x) (vectorFromDiag y)

instance (AEq e, BLAS1 e) => AEq (Diag Vector (n,n) e) where
    (===) x y = (===) (vectorFromDiag x) (vectorFromDiag y)
    (~==) x y = (~==) (vectorFromDiag x) (vectorFromDiag y)

{-
instance (BLAS1 e) => ITensor (Diag Vector) (Int,Int) e where
    zero (m,n) | m /= n = error "tried to make a non-square diagonal matrix"
               | otherwise = coerceDiag $ diagFromVector $ zero n

    constant (m,n) e | m /= n = error "tried to make a non-square diagonal matrix"
                     | otherwise = coerceDiag $ diagFromVector $ constant n e
               
    size = size . vectorFromDiag . coerceDiag
    
    assocs a =
        let ies = assocs $ vectorFromDiag $ coerceDiag a
        in map (\(i,e) -> ((i,i),e)) ies
    
    (//) = replaceHelp (//)
    
    tmap f a = (coerceDiag . diagFromVector) 
                   (tmap f $ vectorFromDiag $ coerceDiag a)
    
    unsafeAt a (i,j) | i /= j = 0
                     | otherwise = unsafeAt (vectorFromDiag $ coerceDiag a) i
                     
    unsafeReplace = replaceHelp unsafeReplace
    
replaceHelp :: (BLAS1 e) => 
       (Vector n e -> [(Int,e)] -> Vector n e) 
    -> Diag Vector nn e -> [((Int,Int),e)] -> Diag Vector nn e
replaceHelp f a ijes =
    let iies = filter (\((i,j),_) -> i == j) ijes
        ies  = map (\((i,_),e) -> (i,e)) iies
        x'   = f (vectorFromDiag $ coerceDiag a) ies
    in coerceDiag $ diagFromVector x'
    
    
instance (BLAS1 e) => ReadTensor (Diag Vector) (Int,Int) e IO where
    getSize = getSize . toVector
    
    newCopy a = do
        x' <- newCopy $ toVector a
        return $ fromVector x'
    
    unsafeReadElem a (i,j)
        | i /= j    = return 0
        | otherwise = unsafeReadElem (toVector a) i
        

        
instance (BLAS1 e) => MTensor (DiagMatrix Mut (n,n)) (Int,Int) e IO where
    setZero = setZero . toVector
    
    setConstant k = setConstant k . toVector
    
    canModifyElem a (i,j) = return (i == j && i >= 0 && i < numRows a)
    
    unsafeWriteElem a (i,_) = unsafeWriteElem (toVector a) i
    
    modifyWith f a = modifyWith f (toVector a)
    
-}
