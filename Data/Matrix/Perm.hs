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

    ) where

import Control.Monad ( forM_ )
import Foreign ( peek, poke )

import BLAS.Elem ( Elem, BLAS1 )
import qualified BLAS.Elem as E

import BLAS.Internal ( checkMatVecMult, checkMatMatMult )

import BLAS.Matrix.Base hiding ( Matrix )
import qualified BLAS.Matrix.Base as Base
import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly
import BLAS.Matrix.Solve.Immutable
import BLAS.Matrix.Solve.ReadOnly


import Data.AEq

import Data.Vector.Dense.IO ( IOVector, newVector_, dim, isConj, conj,
    unsafeSwapVectors, unsafeCopyVector, unsafeWithElemPtr )
import qualified Data.Vector.Dense.IO as V
    
import Data.Matrix.Dense.IO ( IOMatrix, newMatrix_, shape, unsafeRow )
import qualified Data.Matrix.Dense.IO as M

import Data.Permutation ( Permutation )
import qualified Data.Permutation as P

import System.IO.Unsafe ( unsafePerformIO )
import Unsafe.Coerce ( unsafeCoerce )

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

doApply :: (Elem e) => Perm (n,n) e -> IOVector n e -> IO ()
doApply p x = 
    checkMatVecMult (numRows p, numCols p) (dim x) $
        unsafeDoApply p x

unsafeDoApply :: (Elem e) => Perm (n,n) e -> IOVector n e -> IO ()
unsafeDoApply   (I _)   _ = return ()
unsafeDoApply p@(P _ _) x
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

doApplyMat :: (BLAS1 e) => Perm (m,m) e -> IOMatrix (m,n) e -> IO ()
doApplyMat p a = 
    checkMatMatMult (numRows p, numCols p) (shape a) $
        unsafeDoApplyMat p a

unsafeDoApplyMat :: (BLAS1 e) => Perm (m,m) e -> IOMatrix (m,n) e -> IO ()
unsafeDoApplyMat   (I _)   _ = return ()
unsafeDoApplyMat p@(P _ _) a
    | isHerm p  = P.invertWith swap sigma
    | otherwise = P.applyWith swap sigma
  where
    sigma = baseOf p
    swap i j = unsafeSwapVectors (unsafeRow a i) (unsafeRow a j)

          
instance Base.Matrix Perm where
    numRows (P sigma _) = P.size sigma
    numRows (I n)       = n
    
    numCols a = numRows a
    
    herm a@(P _ _) = a{ isHerm=(not . isHerm) a }
    herm a@(I _)   = (unsafeCoerce a)


instance (BLAS1 e) => IMatrix Perm e where
          
         
instance (BLAS1 e) => RMatrix Perm e where
    unsafeGetSApply k (I _) x = do
        y <- V.getScaled k x
        return (unsafeCoerce y)
    
    unsafeGetSApply k p x 
        | isConj x = do
            y <- unsafeGetSApply (E.conj k) p (conj x)
            return (conj y)
        | otherwise =
            let n     = dim x
                sigma = baseOf p
            in do
                y <- newVector_ n
                forM_ [0..(n-1)] $ \i ->
                    let i' = P.unsafeApply sigma i
                    in case (isHerm p) of
                           False -> V.unsafeReadElem x i  
                                    >>= V.unsafeWriteElem (V.unsafeThaw y) i'
                           True  -> V.unsafeReadElem x i' 
                                    >>= V.unsafeWriteElem (V.unsafeThaw y) i
                V.scaleBy k (V.unsafeThaw y)
                return y

    unsafeGetSApplyMat k (I _) a = do
        b <- M.getScaled k a
        return (unsafeCoerce b)
    
    unsafeGetSApplyMat k p a =
        let (m,n) = shape a
            sigma = baseOf p
        in do
            b <- newMatrix_ (m,n)
            forM_ [0..(m-1)] $ \i ->
                let i' = P.unsafeApply sigma i
                in case (isHerm p) of
                       False -> unsafeCopyVector (V.unsafeThaw $ unsafeRow b i') 
                                                 (unsafeRow a i)
                       True  -> unsafeCopyVector (V.unsafeThaw $ unsafeRow b i)  
                                                 (unsafeRow a i')
            M.scaleBy k (M.unsafeThaw b)
            return b


instance (BLAS1 e) => ISolve Perm e where
    unsafeSolve p    = unsafeApply (herm p)
    unsafeSolveMat p = unsafeApplyMat (herm p)
    
    
instance (BLAS1 e) => RSolve Perm e where    
    unsafeGetSolve p    = unsafeGetApply (herm p)
    unsafeGetSolveMat p = unsafeGetApplyMat (herm p)


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
    