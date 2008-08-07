{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Diag
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Diag (
    module BLAS.Matrix,
    module BLAS.Tensor,

    -- * The diagonal matrix types
    DiagMatrix,
    Diag,
    IODiag,
    
    -- * Converting to and from @Vector@s
    fromVector,
    toVector,

    -- * In-place operations
    doSolve,
    doSolveMat,
    
    -- * Unsafe operations
    unsafeDoSolve,
    unsafeDoSolveMat,
    
    ) where

import BLAS.Access
import BLAS.Elem ( BLAS1, BLAS2 )
import BLAS.Internal ( checkMatVecSolv, checkMatMatSolv )
import BLAS.Matrix hiding ( Matrix )
import qualified BLAS.Matrix as Base
import BLAS.Tensor

import Control.Monad ( zipWithM_ )

import Data.AEq
import Data.Matrix.Dense.Internal ( IOMatrix, rows )
import qualified Data.Matrix.Dense.Operations as M
import qualified Data.Matrix.Dense.Internal as M
import Data.Vector.Dense.Internal
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V
import Data.Vector.Dense.Operations ( (//=), scaleBy, invScaleBy )

import System.IO.Unsafe
import Unsafe.Coerce

newtype DiagMatrix t mn e = Diag (DVector t mn e)
type Diag   = DiagMatrix Imm
type IODiag = DiagMatrix Mut

coerceDiag :: DiagMatrix t mn e -> DiagMatrix t mn' e
coerceDiag = unsafeCoerce

fromVector :: DVector t n e -> DiagMatrix t (n,n) e
fromVector = Diag . unsafeCoerce

toVector :: DiagMatrix t (n,n) e -> DVector t n e
toVector (Diag x) = unsafeCoerce x

doSolve :: (BLAS2 e) => DiagMatrix t (n,n) e -> IOVector n e -> IO ()
doSolve a x =
    checkMatVecSolv (shape a) (dim x) $ unsafeDoSolve a x
    
unsafeDoSolve :: (BLAS2 e) => DiagMatrix t (n,n) e -> IOVector n e -> IO ()
unsafeDoSolve a x = x //= (toVector a)

doSolveMat :: (BLAS1 e) => DiagMatrix t (m,m) e -> IOMatrix (m,n) e -> IO ()
doSolveMat a b =
    checkMatMatSolv (shape a) (shape b) $ unsafeDoSolveMat a b

unsafeDoSolveMat :: (BLAS1 e) => DiagMatrix t (m,m) e -> IOMatrix (m,n) e -> IO ()
unsafeDoSolveMat a b = do
    ks <- unsafeInterleaveIO $ getElems (toVector a)
    zipWithM_ (\k r -> invScaleBy k r) ks (rows b)


instance (BLAS1 e) => Scalable (DiagMatrix Imm (n,n)) e where
    (*>) k (Diag x) = Diag $ k *> x


replaceHelp :: (BLAS1 e) => 
       (Vector n e -> [(Int,e)] -> Vector n e) 
    -> Diag (n,n) e -> [((Int,Int),e)] -> Diag (n,n) e
replaceHelp f a ijes =
    let iies = filter (\((i,j),_) -> i == j) ijes
        ies  = map (\((i,_),e) -> (i,e)) iies
        x'   = f (toVector a) ies
    in fromVector x'
    
    
instance (BLAS1 e) => ITensor (DiagMatrix Imm (n,n)) (Int,Int) e where
    size = size . toVector
    
    assocs a =
        let ies = assocs $ toVector a
        in map (\(i,e) -> ((i,i),e)) ies
    
    (//) = replaceHelp (//)
    
    amap f a = fromVector (amap f $ toVector a)
    
    unsafeAt a (i,j) | i /= j = 0
                     | otherwise = unsafeAt (toVector a) i
                     
    unsafeReplace = replaceHelp unsafeReplace
    
    
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
    
        
instance Base.Matrix (DiagMatrix t) where
    numRows (Diag x) = dim x
    numCols = numCols

    herm (Diag x) = unsafeCoerce $ Diag (conj x)


instance (BLAS2 e) => IMatrix (DiagMatrix Imm) e

instance (BLAS2 e) => RMatrix (DiagMatrix t) e where
    unsafeDoSApplyAdd alpha a x beta y = do
        x' <- newCopy x
        unsafeDoApply_ (coerceDiag a) (V.unsafeThaw x')
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



instance (BLAS2 e) => ISolve (DiagMatrix Imm) e where
    unsafeSolve a x = unsafePerformIO $ unsafeGetSolve a x
    {-# NOINLINE unsafeSolve #-}
    
    unsafeSolveMat a b = unsafePerformIO $ unsafeGetSolveMat a b
    {-# NOINLINE unsafeSolveMat #-}


instance (BLAS2 e) => RSolve (DiagMatrix Imm) e where
    unsafeGetSolve a x = do
        x' <- newCopy x
        unsafeDoSolve (unsafeCoerce a) (V.unsafeThaw x')
        return (unsafeCoerce x')

    unsafeGetSolveMat a b = do
        b' <- newCopy b
        unsafeDoSolveMat (unsafeCoerce a) (M.unsafeThaw b')
        return (unsafeCoerce b')


instance (Show e, BLAS1 e) => Show (DiagMatrix Imm (n,n) e) where
    show (Diag x) = "fromVector (" ++ show x ++ ")"


instance (Eq e, BLAS1 e) => Eq (DiagMatrix Imm (n,n) e) where
    (==) (Diag x) (Diag y) = (==) x y


instance (AEq e, BLAS1 e) => AEq (DiagMatrix Imm (n,n) e) where
    (===) (Diag x) (Diag y) = (===) x y
    (~==) (Diag x) (Diag y) = (~==) x y
