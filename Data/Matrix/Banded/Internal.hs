{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Internal (
    -- * Banded matrix type
    Banded(..),

    -- * Banded shape
    module BLAS.Tensor.Base,
    module BLAS.Matrix.Base,
    bandwidth,
    numLower,
    numUpper,
    coerceBanded,

    -- * Creating banded matrices
    banded,
    listsBanded,
    unsafeBanded,

    -- * Reading banded matrix elements
    module BLAS.Tensor.Immutable,
    
    -- * Special banded matrices
    zeroBanded,
    constantBanded,

    -- * Vector views
    diagBanded,
    unsafeDiagBanded,

    -- * Converting to lists
    listsFromBanded,

    -- * Low-level properties
    ldaOfBanded,
    isHermBanded,
    
    -- * Matrix and vector multiplication
    module BLAS.Matrix.Immutable,

    ) where

import Data.AEq
import System.IO.Unsafe


import BLAS.Internal ( diagLen, checkedDiag, inlinePerformIO )
import BLAS.Elem( Elem, BLAS1, BLAS2 )
import BLAS.Tensor.Base
import BLAS.Tensor.Immutable
import BLAS.Tensor.Read
import BLAS.UnsafeIOToM

import BLAS.Matrix.Base hiding ( BaseMatrix )
import BLAS.Matrix.Immutable
import BLAS.Matrix.Mutable
import qualified BLAS.Matrix.Base as BLAS

import Data.Ix( inRange, range )
import Data.Matrix.Banded.Class.Internal( BaseBanded(..), ReadBanded,
    IOBanded, coerceBanded, numLower, numUpper, bandwidth, isHermBanded,
    shapeBanded, boundsBanded, ldaOfBanded, gbmv, gbmm, unsafeGetRowBanded,
    unsafeGetColBanded )
import Data.Matrix.Banded.Class.Creating( newListsBanded, unsafeNewBanded, 
    newBanded )
import Data.Matrix.Banded.Class.Elements( writeElem, unsafeWriteElem )
import Data.Matrix.Banded.Class.Special( newZeroBanded, newConstantBanded )
import Data.Matrix.Banded.Class.Views( unsafeDiagViewBanded )
import Data.Matrix.Banded.Class.Copying( newCopyBanded )

import Data.Vector.Dense( Vector, zeroVector )
import Data.Vector.Dense.ST( runSTVector )
import Data.Matrix.Dense.ST( runSTMatrix )

newtype Banded mn e = B (IOBanded mn e)

unsafeFreezeIOBanded :: IOBanded mn e -> Banded mn e
unsafeFreezeIOBanded = B

unsafeThawIOBanded :: Banded mn e -> IOBanded mn e
unsafeThawIOBanded (B a) = a


liftBanded :: (IOBanded mn e -> a) -> Banded mn e -> a
liftBanded f (B x) = f x
{-# INLINE liftBanded #-}


-- liftBanded2 :: 
--     (IOBanded mn e -> IOBanded mn e -> a) -> 
--         Banded mn e -> Banded mn e -> a
-- liftBanded2 f x = liftBanded (liftBanded f x)
-- {-# INLINE liftBanded2 #-}
-- 
-- unsafeLiftBanded :: (IOBanded mn e -> IO a) -> Banded mn e -> a
-- unsafeLiftBanded f = unsafePerformIO . liftBanded f
-- {-# NOINLINE unsafeLiftBanded #-}
-- 
-- unsafeLiftBanded2 :: 
--     (IOBanded mn e -> IOBanded mn e -> IO a) -> 
--         Banded mn e -> Banded mn e -> a
-- unsafeLiftBanded2 f x y = unsafePerformIO $ liftBanded2 f x y
-- {-# NOINLINE unsafeLiftBanded2 #-}


inlineLiftBanded :: (IOBanded n e -> IO a) -> Banded n e -> a
inlineLiftBanded f = inlinePerformIO . liftBanded f
{-# INLINE inlineLiftBanded #-}


banded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (m,n) e
banded mn kl ijes = 
    unsafeFreezeIOBanded $ unsafePerformIO $ newBanded mn kl ijes
{-# NOINLINE banded #-}

unsafeBanded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (m,n) e
unsafeBanded mn kl ijes = 
    unsafeFreezeIOBanded $ unsafePerformIO $ unsafeNewBanded mn kl ijes
{-# NOINLINE unsafeBanded #-}

listsBanded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> [[e]] -> Banded (m,n) e
listsBanded mn kl xs = 
    unsafeFreezeIOBanded $ unsafePerformIO $ newListsBanded mn kl xs
{-# NOINLINE listsBanded #-}

zeroBanded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> Banded (m,n) e
zeroBanded mn kl =
    unsafeFreezeIOBanded $ unsafePerformIO $ newZeroBanded mn kl
{-# NOINLINE zeroBanded #-}

constantBanded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> e -> Banded (m,n) e
constantBanded mn kl e =
    unsafeFreezeIOBanded $ unsafePerformIO $ newConstantBanded mn kl e
{-# NOINLINE constantBanded #-}

-- | Get a the given diagonal in a banded matrix.  Negative indices correspond 
-- to sub-diagonals.
diagBanded :: (BLAS1 e) => Banded mn e -> Int -> Vector k e
diagBanded a = checkedDiag (shape a) (unsafeDiagBanded a)

-- | Same as 'diagBanded' but index is not range-checked.
unsafeDiagBanded :: (BLAS1 e) => Banded mn e -> Int -> Vector k e
unsafeDiagBanded a i 
    | inRange (bandwidth a) i = unsafeDiagViewBanded a i
    | otherwise               = zeroVector $ diagLen (shape a) i


instance (Elem e) => BaseTensor Banded (Int,Int) e where
    shape  = shapeBanded . unsafeThawIOBanded
    bounds = boundsBanded . unsafeThawIOBanded

instance (BLAS1 e) => ITensor Banded (Int,Int) e where
    (//)          = replaceHelp writeElem
    unsafeReplace = replaceHelp unsafeWriteElem
    
    unsafeAt x i  = inlineLiftBanded (flip unsafeReadElem i) x
    {-# INLINE unsafeAt #-}
    
    size          = inlineLiftBanded getSize
    elems         = inlineLiftBanded getElems
    indices       = inlineLiftBanded getIndices
    assocs        = inlineLiftBanded getAssocs

    tmap f a      = coerceBanded $ listsBanded mn bw (map (map f) es)
      where (mn,bw,es) = listsFromBanded a

listsFromBanded :: (BLAS1 e) => Banded mn e -> ((Int,Int), (Int,Int),[[e]])
listsFromBanded a = ( (m,n)
            , (kl,ku)
            , map paddedDiag [(-kl)..ku]
            )
  where
    (m,n)   = shape a
    (kl,ku) = (numLower a, numUpper a)
    
    padBegin i   = replicate (max (-i) 0)    0
    padEnd   i   = replicate (max (m-n+i) 0) 0
    paddedDiag i = (  padBegin i
                   ++ elems (unsafeDiagViewBanded a i) 
                   ++ padEnd i 
                   )

replaceHelp :: (BLAS1 e) => 
       (IOBanded mn e -> (Int,Int) -> e -> IO ())
    -> Banded mn e -> [((Int,Int), e)] -> Banded mn e
replaceHelp set x ies =
    unsafeFreezeIOBanded $ unsafePerformIO $ do
        y  <- newCopyBanded (unsafeThawIOBanded x)
        mapM_ (uncurry $ set y) ies
        return y
{-# NOINLINE replaceHelp #-}


instance (BLAS1 e, Monad m) => ReadTensor Banded (Int,Int) e m where
    getSize        = return . size
    getAssocs      = return . assocs
    getIndices     = return . indices
    getElems       = return . elems
    getAssocs'     = getAssocs
    getIndices'    = getIndices
    getElems'      = getElems
    unsafeReadElem x i = return (unsafeAt x i)

instance (Elem e) => BLAS.BaseMatrix Banded e where
    herm (B a) = B (herm a)
    
instance (Elem e) => BaseBanded Banded Vector e where
    bandedViewArray f p m n kl ku l h = B $ bandedViewArray f p m n kl ku l h
    arrayFromBanded (B a )            = arrayFromBanded a

instance (BLAS1 e, UnsafeIOToM m) => 
    ReadBanded Banded Vector e m where

instance (BLAS2 e) => IMatrix Banded e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    
    unsafeRow a i             = runSTVector $ unsafeGetRow a i
    unsafeCol a j             = runSTVector $ unsafeGetCol a j

instance (BLAS2 e, UnsafeIOToM m) => MMatrix Banded e m where
    unsafeDoSApplyAdd    = gbmv
    unsafeDoSApplyAddMat = gbmm
    unsafeGetRow         = unsafeGetRowBanded
    unsafeGetCol         = unsafeGetColBanded


instance (BLAS1 e) => Show (Banded mn e) where
    show a 
        | isHermBanded a = 
           "herm (" ++ show (herm $ coerceBanded a) ++ ")"
        | otherwise = 
             let (mn,kl,es) = listsFromBanded a 
             in "listsBanded " ++ show mn ++ " " ++ show kl ++ " " ++ show es

compareHelp :: (BLAS1 e) => 
    (e -> e -> Bool) -> Banded mn e -> Banded mn e -> Bool
compareHelp cmp a b
    | shape a /= shape b =
        False
    | isHermBanded a == isHermBanded b && bandwidth a == bandwidth b =
        let elems' = if isHermBanded a then elems . herm .coerceBanded
                                       else elems
        in
            and $ zipWith cmp (elems' a) (elems' b)
    | otherwise =
        let l = max (numLower a) (numLower b)
            u = max (numUpper a) (numUpper b)
        in
            and $ zipWith cmp (diagElems (-l,u) a) (diagElems (-l,u) b)
  where
    diagElems bw c = concatMap elems [ diagBanded c i | i <- range bw ]

instance (BLAS1 e, Eq e) => Eq (Banded mn e) where
    (==) = compareHelp (==)

instance (BLAS1 e, AEq e) => AEq (Banded mn e) where
    (===) = compareHelp (===)
    (~==) = compareHelp (~==)
