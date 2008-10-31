{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Perm
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Perm (
    module BLAS.Matrix,

    Perm(..),
    
    -- * The identity permutation
    identity,
    
    -- * Converting to/from @Permutation@s
    fromPermutation,
    toPermutation,
    
    -- * Coercing
    coercePerm

    ) where

import Control.Monad ( forM_ )
import Data.AEq

import BLAS.Elem
import BLAS.Matrix
import BLAS.Tensor
import BLAS.UnsafeIOToM

import Data.Matrix.Dense.Class( ReadMatrix, WriteMatrix, unsafeRowView,
    unsafeAxpyMatrix, coerceMatrix )
import Data.Vector.Dense.Class( ReadVector, WriteVector, dim, isConj, 
    scaleBy, unsafeAxpyVector, unsafeSwapVector, coerceVector )
import Data.Vector.Dense.ST( runSTVector )
import Data.Matrix.Dense.ST( runSTMatrix )

{-


import Data.Vector.Dense.IO ( dim, isConj, conj, unsafeSwapVectors, 
    unsafeCopyVector, unsafeWithElemPtr, unsafeReadElem, unsafeWriteElem, 
    coerceVector )
import qualified Data.Vector.Dense.IO as V
    
import Data.Matrix.Dense.IO ( unsafeRow, unsafeCopyMatrix, coerceMatrix )
import qualified Data.Matrix.Dense.IO as M
-}

import Data.Permutation ( Permutation )
import qualified Data.Permutation as P

import Unsafe.Coerce


data Perm mn e = 
      P { baseOf :: !Permutation
        , isHerm :: !Bool
        }
    | I !Int

identity :: Int -> Perm (n,n) e
identity = I

fromPermutation :: Permutation -> Perm (n,n) e
fromPermutation = flip P False

toPermutation :: Perm (n,n) e -> Permutation
toPermutation (I n)       = P.identity n
toPermutation (P sigma h) = if h then P.inverse sigma else sigma

coercePerm :: Perm mn e -> Perm mn' e
coercePerm = unsafeCoerce

instance BaseTensor Perm (Int,Int) e where
    shape (P sigma _) = (n,n) where n = P.size sigma
    shape (I n)       = (n,n)
    
    bounds p = ((0,0), (m-1,n-1)) where (m,n) = shape p
    
instance BaseMatrix Perm e where
    herm a@(P _ _) = a{ isHerm=(not . isHerm) a }
    herm a@(I _)   = coercePerm a


{-          
instance (BLAS1 e) => IMatrix Perm e where
instance (BLAS1 e) => RMatrix Perm e where

-}

instance (BLAS1 e) => IMatrix Perm e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b

instance (BLAS1 e, UnsafeIOToM m) => MMatrix Perm e m where
    unsafeDoSApplyAdd    = unsafeDoSApplyAddPerm
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatPerm
    unsafeDoSApply_      = unsafeDoSApplyPerm_
    unsafeDoSApplyMat_   = unsafeDoSApplyMatPerm_


unsafeDoSApplyPerm_ :: (WriteVector y e m, BLAS1 e) => 
    e -> Perm (k,k) e -> y k e -> m ()
unsafeDoSApplyPerm_ alpha   (I _)   x = scaleBy alpha x
unsafeDoSApplyPerm_ alpha p@(P _ _) x
    | isHerm p  = P.invertWith swap sigma >> scaleBy alpha x
    | otherwise = P.applyWith swap sigma  >> scaleBy alpha x
  where
    sigma = baseOf p
    swap  = unsafeSwapElem x

unsafeDoSApplyMatPerm_ :: (WriteMatrix c z e m, BLAS1 e) => 
    e -> Perm (k,k) e -> c (k,l) e -> m ()
unsafeDoSApplyMatPerm_ alpha   (I _)   a = scaleBy alpha a
unsafeDoSApplyMatPerm_ alpha p@(P _ _) a
    | isHerm p  = P.invertWith swap sigma >> scaleBy alpha a
    | otherwise = P.applyWith  swap sigma >> scaleBy alpha a
  where
    sigma    = baseOf p
    swap i j = unsafeSwapVector (unsafeRowView a i) (unsafeRowView a j)


unsafeDoSApplyAddPerm :: (ReadVector x e m, WriteVector y e m, BLAS1 e) =>
    e -> Perm (k,l) e -> x l e -> e -> y k e -> m ()
unsafeDoSApplyAddPerm alpha (I _) x beta y = do
    scaleBy beta y
    unsafeAxpyVector alpha (coerceVector x) y
unsafeDoSApplyAddPerm alpha p x beta y
    | isConj x =
        unsafeDoSApplyAddPerm (conj alpha) p (conj x) (conj beta) (conj y)
    | otherwise =
        let n     = dim x
            sigma = baseOf p
        in do
            scaleBy beta y
            forM_ [0..(n-1)] $ \i ->
                let i' = P.unsafeApply sigma i
                in case (isHerm p) of
                       False -> do
                           e <- unsafeReadElem x i  
                           f <- unsafeReadElem y i'
                           unsafeWriteElem y i' (alpha*e + f)
                           
                       True  -> do
                           e <- unsafeReadElem x i'
                           f <- unsafeReadElem y i
                           unsafeWriteElem y i (alpha*e + f)

unsafeDoSApplyAddMatPerm :: (ReadMatrix b x e m, WriteMatrix c y e m, BLAS1 e) =>
    e -> Perm (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
unsafeDoSApplyAddMatPerm alpha (I _) b beta c = do
    scaleBy beta c
    unsafeAxpyMatrix alpha b (coerceMatrix c)
unsafeDoSApplyAddMatPerm alpha p b beta c =
    let m     = numCols p
        sigma = baseOf p
    in do
        scaleBy beta c
        forM_ [0..(m-1)] $ \i ->
            let i' = P.unsafeApply sigma i
            in case (isHerm p) of
                   False -> unsafeAxpyVector alpha (unsafeRowView b i)
                                (unsafeRowView c i') 
                   True  -> unsafeAxpyVector alpha (unsafeRowView b i') 
                                (unsafeRowView c i)

instance (BLAS1 e) => ISolve Perm e where
    unsafeSSolve alpha a y    = runSTVector $ unsafeGetSSolve    alpha a y
    unsafeSSolveMat alpha a c = runSTMatrix $ unsafeGetSSolveMat alpha a c
    
instance (BLAS1 e, UnsafeIOToM m) => MSolve Perm e m where    
    unsafeDoSSolve_ alpha p       = unsafeDoSApplyPerm_ alpha (herm p)
    unsafeDoSSolveMat_ alpha p    = unsafeDoSApplyMatPerm_ alpha (herm p)
    unsafeDoSSolve alpha p x y    = unsafeDoSApplyAddPerm alpha (coercePerm $ herm p) x 0 y
    unsafeDoSSolveMat alpha p a b = unsafeDoSApplyAddMatPerm alpha (coercePerm $ herm p) a 0 b


instance (Elem e) => Show (Perm (n,n) e) where
    show (I n)         = "identity " ++ show n
    show p | isHerm p  = "herm (" ++ show (herm p) ++ ")"
           | otherwise = "fromPermutation (" ++ show (baseOf p) ++ ")"
    
    
instance Eq (Perm (n,n) e) where
    (==) (I n) (I n') = n == n'
    (==) (I n) p
        | isHerm p   = (==) (I n) (herm p)
        | otherwise  = (==) (fromPermutation $ P.identity n) p
    (==) p     (I n) = (==) (I n) p

    (==) (P sigma h) (P sigma' h') 
        | h == h'   = sigma == sigma'
        | otherwise = P.size sigma == P.size sigma'
                      && sigma == (P.inverse sigma')


instance AEq (Perm (n,n) e) where
    (===) = (==)
    (~==) = (==)
    