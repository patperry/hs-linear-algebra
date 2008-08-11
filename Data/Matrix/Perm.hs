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
    module BLAS.Matrix.Base,
    module BLAS.Matrix.Immutable,
    module BLAS.Matrix.ReadOnly,
    module BLAS.Matrix.Solve.Immutable,
    module BLAS.Matrix.Solve.ReadOnly,

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
import Foreign ( peek, poke )

import BLAS.Elem ( Elem, BLAS1 )
import qualified BLAS.Elem as E

import BLAS.Matrix.Base hiding ( Matrix )
import qualified BLAS.Matrix.Base as Base
import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly
import BLAS.Matrix.Solve.Immutable
import BLAS.Matrix.Solve.ReadOnly

import Data.AEq

import Data.Vector.Dense.IO ( dim, isConj, conj, unsafeSwapVectors, 
    unsafeCopyVector, unsafeWithElemPtr, unsafeReadElem, unsafeWriteElem, 
    coerceVector )
import qualified Data.Vector.Dense.IO as V
    
import Data.Matrix.Dense.IO ( unsafeRow, unsafeCopyMatrix, coerceMatrix )
import qualified Data.Matrix.Dense.IO as M

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

          
instance Base.Matrix Perm where
    numRows (P sigma _) = P.size sigma
    numRows (I n)       = n
    
    numCols a = numRows a
    
    herm a@(P _ _) = a{ isHerm=(not . isHerm) a }
    herm a@(I _)   = (unsafeCoerce a)


instance (BLAS1 e) => IMatrix Perm e where
          
         
instance (BLAS1 e) => RMatrix Perm e where
    unsafeDoApply_   (I _)   _ = return ()
    unsafeDoApply_ p@(P _ _) x
        | isHerm p  = P.invertWith swap sigma
        | otherwise = P.applyWith swap sigma
      where
        sigma = baseOf p
        swap i j = 
            unsafeWithElemPtr x i $ \pI ->
                unsafeWithElemPtr x j $ \pJ -> do
                    xi <- peek pI
                    xj <- peek pJ
                    poke pI xj
                    poke pJ xi


    unsafeDoApplyMat_   (I _)   _ = return ()
    unsafeDoApplyMat_ p@(P _ _) a
        | isHerm p  = P.invertWith swap sigma
        | otherwise = P.applyWith swap sigma
      where
        sigma = baseOf p
        swap i j = unsafeSwapVectors (unsafeRow a i) (unsafeRow a j)


    unsafeDoSApply k (I _) x y = do
        unsafeCopyVector y (coerceVector x)
        V.scaleBy k y
    unsafeDoSApply k p x y
        | isConj x =
            unsafeDoSApply (E.conj k) p (conj x) (conj y)
        | otherwise =
            let n     = dim x
                sigma = baseOf p
            in do
                forM_ [0..(n-1)] $ \i ->
                    let i' = P.unsafeApply sigma i
                    in case (isHerm p) of
                           False -> unsafeReadElem x i  
                                    >>= unsafeWriteElem y i'
                           True  -> unsafeReadElem x i' 
                                    >>= unsafeWriteElem y i
                V.scaleBy k y


    unsafeDoSApplyMat alpha (I _) b c = do
        unsafeCopyMatrix c (coerceMatrix b)
        M.scaleBy alpha c
    unsafeDoSApplyMat alpha p b c =
        let m     = numCols p
            sigma = baseOf p
        in do
            forM_ [0..(m-1)] $ \i ->
                let i' = P.unsafeApply sigma i
                in case (isHerm p) of
                       False -> unsafeCopyVector (unsafeRow c i') 
                                                 (unsafeRow b i)
                       True  -> unsafeCopyVector (unsafeRow c i)
                                                 (unsafeRow b i')
            M.scaleBy alpha c


instance (BLAS1 e) => ISolve Perm e where
    
instance (BLAS1 e) => RSolve Perm e where    
    unsafeDoSolve_ p          = unsafeDoApply_ (herm p)
    unsafeDoSolveMat_ p       = unsafeDoApplyMat_ (herm p)
    unsafeDoSSolve alpha p    = unsafeDoSApply alpha (coercePerm $ herm p)
    unsafeDoSSolveMat alpha p = unsafeDoSApplyMat alpha (coercePerm $ herm p)


instance (Elem e) => Show (Perm (n,n) e) where
    show (I n)         = "identity " ++ show n
    show p | isHerm p  = "herm (" ++ show (herm p) ++ ")"
           | otherwise = "fromPermutation (" ++ show (baseOf p) ++ ")"
    
    
instance (Elem e) => Eq (Perm (n,n) e) where
    (==) (I n) (I n') = n == n'
    (==) (I n) p
        | isHerm p   = (==) (I n) (herm p)
        | otherwise  = (==) (fromPermutation $ P.identity n) p
    (==) p     (I n) = (==) (I n) p

    (==) (P sigma h) (P sigma' h') 
        | h == h'   = sigma == sigma'
        | otherwise = P.size sigma == P.size sigma'
                      && sigma == (P.inverse sigma')


instance (Elem e) => AEq (Perm (n,n) e) where
    (===) = (==)
    (~==) = (==)
    