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
    module BLAS.Matrix.Shaped,
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
import Foreign( Storable )
import System.IO.Unsafe


import BLAS.Internal ( diagLen, checkedDiag, inlinePerformIO )
import BLAS.Elem( BLAS1, BLAS2, BLAS3 )
import BLAS.Tensor.Base
import BLAS.Tensor.Immutable
import BLAS.Tensor.Read
import BLAS.UnsafeIOToM

import BLAS.Matrix.Shaped
import BLAS.Matrix.Immutable
import BLAS.Matrix.Mutable
import BLAS.Matrix.Solve.Immutable
import BLAS.Matrix.Solve.Mutable

import Data.Ix( inRange, range )
import Data.Matrix.Herm
import Data.Matrix.Tri.Internal
import Data.Matrix.Banded.Class.Internal( BaseBanded_(..), BaseBanded,
    ReadBanded,
    IOBanded, coerceBanded, numLower, numUpper, bandwidth, isHermBanded,
    shapeBanded, boundsBanded, ldaOfBanded, gbmv, gbmm, hbmv', hbmm',
    tbmv, tbmm, tbmv', tbmm', unsafeDoSSolveTriBanded, 
    unsafeDoSSolveMatTriBanded, tbsv, tbsm, unsafeGetRowBanded, 
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


banded :: (BLAS2 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (m,n) e
banded mn kl ijes = 
    unsafeFreezeIOBanded $ unsafePerformIO $ newBanded mn kl ijes
{-# NOINLINE banded #-}

unsafeBanded :: (BLAS2 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (m,n) e
unsafeBanded mn kl ijes = 
    unsafeFreezeIOBanded $ unsafePerformIO $ unsafeNewBanded mn kl ijes
{-# NOINLINE unsafeBanded #-}

listsBanded :: (BLAS2 e) => (Int,Int) -> (Int,Int) -> [[e]] -> Banded (m,n) e
listsBanded mn kl xs = 
    unsafeFreezeIOBanded $ unsafePerformIO $ newListsBanded mn kl xs
{-# NOINLINE listsBanded #-}

zeroBanded :: (BLAS2 e) => (Int,Int) -> (Int,Int) -> Banded (m,n) e
zeroBanded mn kl =
    unsafeFreezeIOBanded $ unsafePerformIO $ newZeroBanded mn kl
{-# NOINLINE zeroBanded #-}

constantBanded :: (BLAS2 e) => (Int,Int) -> (Int,Int) -> e -> Banded (m,n) e
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


instance (Storable e) => BaseTensor Banded (Int,Int) e where
    shape  = shapeBanded . unsafeThawIOBanded
    bounds = boundsBanded . unsafeThawIOBanded

instance (BLAS2 e) => ITensor Banded (Int,Int) e where
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

replaceHelp :: (BLAS2 e) => 
       (IOBanded mn e -> (Int,Int) -> e -> IO ())
    -> Banded mn e -> [((Int,Int), e)] -> Banded mn e
replaceHelp set x ies =
    unsafeFreezeIOBanded $ unsafePerformIO $ do
        y  <- newCopyBanded (unsafeThawIOBanded x)
        mapM_ (uncurry $ set y) ies
        return y
{-# NOINLINE replaceHelp #-}


instance (BLAS2 e, Monad m) => ReadTensor Banded (Int,Int) e m where
    getSize        = return . size
    getAssocs      = return . assocs
    getIndices     = return . indices
    getElems       = return . elems
    getAssocs'     = getAssocs
    getIndices'    = getIndices
    getElems'      = getElems
    unsafeReadElem x i = return (unsafeAt x i)

instance (Storable e) => MatrixShaped Banded e where
    herm (B a) = B (herm a)
    
instance (Storable e) => BaseBanded_ Banded e where
    type VectorViewB Banded = Vector
    bandedViewArray f p m n kl ku l h = B $ bandedViewArray f p m n kl ku l h
    arrayFromBanded (B a )            = arrayFromBanded a

instance (Storable e) => BaseBanded Banded e

instance (BLAS2 e, UnsafeIOToM m) => ReadBanded Banded e m where

instance (BLAS3 e) => IMatrix Banded e where
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

compareHelp :: (BLAS2 e) => 
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

instance (BLAS2 e, Eq e) => Eq (Banded mn e) where
    (==) = compareHelp (==)

instance (BLAS2 e, AEq e) => AEq (Banded mn e) where
    (===) = compareHelp (===)
    (~==) = compareHelp (~==)

instance (BLAS3 e) => IMatrix (Herm Banded) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    

instance (BLAS2 e, UnsafeIOToM m) => MMatrix (Herm Banded) e m where
    unsafeDoSApplyAdd    = hbmv'
    unsafeDoSApplyAddMat = hbmm'

instance (BLAS3 e) => IMatrix (Tri Banded) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    

instance (BLAS2 e, UnsafeIOToM m) => MMatrix (Tri Banded) e m where
    unsafeDoSApply_      = tbmv
    unsafeDoSApplyMat_   = tbmm
    unsafeDoSApplyAdd    = tbmv'
    unsafeDoSApplyAddMat = tbmm'

instance (BLAS3 e) => ISolve (Tri Banded) e where
    unsafeSSolve    alpha a y = runSTVector $ unsafeGetSSolve    alpha a y
    unsafeSSolveMat alpha a c = runSTMatrix $ unsafeGetSSolveMat alpha a c

instance (BLAS2 e, UnsafeIOToM m) => MSolve (Tri Banded) e m where
    unsafeDoSSolve     = unsafeDoSSolveTriBanded
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriBanded
    unsafeDoSSolve_    = tbsv
    unsafeDoSSolveMat_ = tbsm
