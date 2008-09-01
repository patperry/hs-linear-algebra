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
    
    ) where

import Control.Monad( zipWithM_ )
import Control.Monad.ST( ST )

import BLAS.Elem( BLAS1 )
import BLAS.Matrix hiding ( BaseMatrix )
import qualified BLAS.Matrix as BLAS
import BLAS.Tensor
import BLAS.UnsafeInterleaveM
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
    conj, dim, coerceVector, scaleBy, unsafeCopyVector, unsafeDivVector )

data Diag x nn e = forall n. Diag (x n e)

coerceDiag :: Diag x mn e -> Diag x mn' e
coerceDiag = unsafeCoerce

diagFromVector :: (BaseVector x e) => x n e -> Diag x (n,n) e
diagFromVector = Diag . coerceVector

vectorFromDiag :: (BaseVector x e) => Diag x (n,n) e -> x n e
vectorFromDiag (Diag x) = coerceVector x

instance (BaseVector x e) => BaseTensor (Diag x) (Int,Int) e where
    shape  (Diag x) = (n,n) where n = dim x
    bounds (Diag x) = ((0,0),(n-1,n-1)) where n = dim x

instance (BaseVector x e) => BLAS.BaseMatrix (Diag x) e where
    herm (Diag x) = Diag (conj x)
    
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
    
        
instance (BLAS2 e) => IMatrix (DiagMatrix Imm) e

instance (BLAS2 e) => RMatrix (DiagMatrix t) e where
    unsafeDoSApplyAdd alpha a x beta y = do
        x' <- newCopy x
        unsafeDoSApply_ 1 (coerceDiag a) (V.unsafeThaw x')
        scaleBy beta y
        V.unsafeAxpy alpha x' (V.coerceVector y)

    unsafeDoSApplyAddMat alpha a b beta c = do
        M.scaleBy beta c
        ks <- unsafeInterleaveIO $ getElems (toVector $ coerceDiag a)
        let (kxs) = zip ks (rows b)
            ys    = rows c
        zipWithM_ (\(k,x) y -> V.unsafeAxpy (alpha*k) x y) kxs ys

    unsafeDoSApply_ alpha a x = do
        V.unsafeTimesEquals x (toVector a)
        V.scaleBy alpha x

    unsafeDoSApplyMat_ alpha a b = do
        ks <- unsafeInterleaveIO $ getElems (toVector a)
        zipWithM_ (\k r -> scaleBy (alpha*k) r) ks (rows b)

-}

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

instance (BLAS1 e, UnsafeIOToM m, UnsafeInterleaveM m) => MSolve (Diag Vector) e m where
    unsafeDoSSolve     = unsafeDoSSolveDiagVector
    unsafeDoSSolveMat  = unsafeDoSSolveMatDiagVector
    unsafeDoSSolve_    = unsafeDoSSolveDiagVector_
    unsafeDoSSolveMat_ = unsafeDoSSolveMatDiagVector_



unsafeDoSSolveDiagVector :: (ReadVector z e m, ReadVector y e m, WriteVector x e m, BLAS1 e) =>
    e -> Diag z (r,s) e -> y r e -> x s e -> m ()
unsafeDoSSolveDiagVector alpha a y x = do
    unsafeCopyVector x (coerceVector y)
    unsafeDoSSolveDiagVector_ alpha (coerceDiag a) x

unsafeDoSSolveMatDiagVector :: (ReadVector x e m, ReadMatrix c z e m, WriteMatrix b y e m, BLAS1 e) =>
    e -> Diag x (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
unsafeDoSSolveMatDiagVector alpha a c b = do
    unsafeCopyMatrix b (coerceMatrix c)
    unsafeDoSSolveMatDiagVector_ alpha (coerceDiag a) b

unsafeDoSSolveDiagVector_ :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    e -> Diag x (k,k) e -> y k e -> m ()
unsafeDoSSolveDiagVector_ alpha a x = do
    scaleBy alpha x
    unsafeDivVector x (vectorFromDiag a)

unsafeDoSSolveMatDiagVector_ :: (ReadVector x e m, WriteMatrix a y e m, BLAS1 e) =>
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
