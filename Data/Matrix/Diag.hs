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
    
    ) where

import BLAS.Elem( BLAS1 )
import BLAS.Matrix hiding ( BaseMatrix )
import qualified BLAS.Matrix as BLAS
import BLAS.Tensor
import Unsafe.Coerce

import Data.Vector.Dense( Vector )
import Data.Vector.Dense.Class( BaseVector, conj, dim, coerceVector )

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
    
{-
    
instance (BLAS1 e) => RTensor (DiagMatrix t (n,n)) (Int,Int) e IO where
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


instance (BLAS2 e) => ISolve (DiagMatrix Imm) e


instance (BLAS2 e) => RSolve (DiagMatrix Imm) e where
    unsafeDoSSolve alpha a y x = do
        unsafeCopyVector x (unsafeCoerce y)
        unsafeDoSSolve_ alpha (coerceDiag a) x
        
    unsafeDoSSolveMat alpha a c b = do
        unsafeCopyMatrix b (unsafeCoerce c)
        unsafeDoSSolveMat_ alpha (coerceDiag a) b
    
    unsafeDoSSolve_ alpha a x = do
        V.scaleBy alpha x
        x //= (toVector a)

    unsafeDoSSolveMat_ alpha a b = do
        M.scaleBy alpha b
        ks <- unsafeInterleaveIO $ getElems (toVector a)
        zipWithM_ (\k r -> invScaleBy k r) ks (rows b)


instance (Show e, BLAS1 e) => Show (DiagMatrix Imm (n,n) e) where
    show (Diag x) = "fromVector (" ++ show x ++ ")"


instance (Eq e, BLAS1 e) => Eq (DiagMatrix Imm (n,n) e) where
    (==) (Diag x) (Diag y) = (==) x y


instance (AEq e, BLAS1 e) => AEq (DiagMatrix Imm (n,n) e) where
    (===) (Diag x) (Diag y) = (===) x y
    (~==) (Diag x) (Diag y) = (~==) x y
-}
