{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses, UndecidableInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Class.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Class.Internal (
    -- * Banded type classes
    BaseBanded(..),
    ReadBanded,
    WriteBanded,

    -- * Low-level Banded properties
    lda,
    isHerm,
    bandedViewMatrix,
    matrixFromBanded,

    -- * Bandwidth properties
    bandwidth,
    numLower,
    numUpper,

    -- * Coercing the Banded shape
    coerceBanded,
    
    ) where

import Control.Monad
import Control.Monad.ST
import Data.Ix
import Data.List( foldl' )
import Foreign
import Unsafe.Coerce

import BLAS.Elem
import BLAS.C.Types
import qualified BLAS.C.Level2 as BLAS
import qualified BLAS.C.Level3 as BLAS
import BLAS.Internal( diagStart, diagLen )
import BLAS.UnsafeIOToM
import BLAS.UnsafeInterleaveM


import BLAS.Numeric
import BLAS.Tensor

import Data.Vector.Dense.Class.Internal( IOVector, STVector,
    BaseVector(..), ReadVector, WriteVector, 
    newCopyVector, unsafeCopyVector, unsafeSwapVector, 
    doConjVector, scaleByVector, shiftByVector, unsafeAxpyVector, 
    unsafeMulVector, unsafeDivVector, withVectorPtr, dim, stride, isConj )

import qualified Data.Matrix.Dense.Class as M
import Data.Matrix.Dense.Class( BaseMatrix, arrayFromMatrix, matrixViewArray,
    coerceMatrix )

import BLAS.Matrix.Base hiding ( BaseMatrix )
import qualified BLAS.Matrix.Base as BLAS


class (BLAS.BaseMatrix a e, BaseVector x e) => 
    BaseBanded a x e | a -> x where
        bandedViewArray :: ForeignPtr e -> Int -> (Int,Int) -> (Int,Int) -> Int -> Bool -> a mn e
        arrayFromBanded :: a mn e -> (ForeignPtr e, Int, (Int,Int), (Int,Int), Int, Bool)

class (Elem e, UnsafeInterleaveM m, ReadTensor a (Int,Int) e m, 
           ReadNumeric a (Int,Int) e m, BaseBanded a x e, 
           ReadVector x e m) => 
    ReadBanded a x e m | a -> x where

class (WriteTensor a (Int,Int) e m, WriteNumeric a (Int,Int) e m,
           WriteVector x e m, ReadBanded a x e m) => 
    WriteBanded a x e m | a -> m, m -> a, a -> x where


------------------------- Basic Banded Properties ---------------------------

fptrOfBanded :: (BaseBanded a x e) => a mn e -> ForeignPtr e
fptrOfBanded a = let (f,_,_,_,_,_) = arrayFromBanded a in f
{-# INLINE fptrOfBanded #-}

offsetOfBanded :: (BaseBanded a x e) => a mn e -> Int
offsetOfBanded a = let (_,o,_,_,_,_) = arrayFromBanded a in o
{-# INLINE offsetOfBanded #-}

size1 :: (BaseBanded a x e) => a mn e -> Int
size1 a = let (_,_,(m,_),_,_,_) = arrayFromBanded a in m
{-# INLINE size1  #-}

size2 :: (BaseBanded a x e) => a mn e -> Int
size2 a = let (_,_,(_,n),_,_,_) = arrayFromBanded a in n
{-# INLINE size2 #-}

lowBW :: (BaseBanded a x e) => a mn e -> Int
lowBW a = let (_,_,_,(kl,_),_,_) = arrayFromBanded a in kl
{-# INLINE lowBW  #-}

upBW :: (BaseBanded a x e) => a mn e -> Int
upBW a = let (_,_,_,(_,ku),_,_) = arrayFromBanded a in ku
{-# INLINE upBW #-}

lda :: (BaseBanded a x e) => a mn e -> Int
lda a = let (_,_,_,_,l,_) = arrayFromBanded a in l
{-# INLINE lda #-}

isHerm :: (BaseBanded a x e) => a mn e -> Bool
isHerm a = let (_,_,_,_,_,h) = arrayFromBanded a in h
{-# INLINE isHerm #-}

matrixFromBanded :: (BaseBanded b x e, BaseMatrix a x e) => 
    b mn e -> ((Int,Int), (Int,Int), a mn' e, Bool)
matrixFromBanded b =
    let (f,o,(m,n),(kl,ku),ld,h) = arrayFromBanded b
        a = matrixViewArray f o (kl+1+ku,n) ld False
    in ((m,n), (kl,ku), a, h)

bandedViewMatrix :: (BaseMatrix a x e, BaseBanded b x e) => 
    (Int,Int) -> (Int,Int) -> a mn e -> Bool -> Maybe (b mn' e)
bandedViewMatrix (m,n) (kl,ku) a h = 
    if M.isHerm a 
        then Nothing
        else let (f,o,(m',n'),ld,_) = arrayFromMatrix a
             in case undefined of
                 _ | m' /= kl+1+ku -> 
                     error $ "bandedViewMatrix:"
                        ++ " number of rows must be equal to number of diagonals"
                 _ | n' /= n ->
                     error $ "bandedViewMatrix:"
                        ++ " numbers of columns must be equal"
                 _ ->
                     Just $ bandedViewArray f o (m,n) (kl,ku) ld h

bandwidth :: (BaseBanded a x e) => a mn e -> (Int,Int)
bandwidth a =
    let (kl,ku) = (numLower a, numUpper a)
    in (negate kl, ku)
{-# INLINE bandwidth #-}

numLower :: (BaseBanded a x e) =>  a mn e -> Int
numLower a | isHerm a  = upBW a
           | otherwise = lowBW a
{-# INLINE numLower #-}

numUpper :: (BaseBanded a x e) =>  a mn e -> Int
numUpper a | isHerm a  = lowBW a
           | otherwise = upBW a
{-# INLINE numUpper #-}


-- | Cast the shape type of the matrix.
coerceBanded :: (BaseBanded a x e) => a mn e -> a mn' e
coerceBanded = unsafeCoerce
{-# INLINE coerceBanded #-}


-------------------------- BaseTensor functions -----------------------------

shapeBanded :: (BaseBanded a x e) => a mn e -> (Int,Int)
shapeBanded a | isHerm a  = (size2 a, size1 a)
              | otherwise = (size1 a, size2 a)
{-# INLINE shapeBanded #-}

boundsBanded :: (BaseBanded a x e) => a mn e -> ((Int,Int), (Int,Int))
boundsBanded a = ((0,0), (m-1,n-1)) where (m,n) = shapeBanded a
{-# INLINE boundsBanded #-}


-------------------------- BaseMatrix functions -----------------------------

hermBanded :: (BaseBanded a x e) => a (m,n) e -> a (n,m) e
hermBanded a = let (f,o,mn,kl,l,h) = arrayFromBanded a
               in bandedViewArray f o mn kl l (not h)
{-# INLINE hermBanded #-}


-------------------------- ReadTensor functions -----------------------------

getSizeBanded :: (ReadBanded a x e m) => a mn e -> m Int
getSizeBanded = return . sizeBanded
{-# INLINE getSizeBanded #-}

getIndicesBanded :: (ReadBanded a x e m) => a mn e -> m [(Int,Int)]
getIndicesBanded = return . indicesBanded
{-# INLINE getIndicesBanded #-}

getElemsBanded :: (ReadBanded a x e m) => a mn e -> m [e]
getElemsBanded a = getAssocsBanded a >>= return . (map snd)

getAssocsBanded :: (ReadBanded a x e m) => a mn e -> m [((Int,Int),e)]
getAssocsBanded a = do
    is <- getIndicesBanded a
    unsafeInterleaveM $ mapM (\i -> unsafeReadElem a i >>= \e -> return (i,e)) is
    
getIndicesBanded' :: (ReadBanded a x e m) => a mn e -> m [(Int,Int)]
getIndicesBanded' = getIndicesBanded
{-# INLINE getIndicesBanded' #-}

getElemsBanded' :: (ReadBanded a x e m) => a mn e -> m [e]
getElemsBanded' a = getAssocsBanded' a >>= return . (map snd)

getAssocsBanded' :: (ReadBanded a x e m) => a mn e -> m [((Int,Int),e)]
getAssocsBanded' a = do
    is <- getIndicesBanded a
    mapM (\i -> unsafeReadElem a i >>= \e -> return (i,e)) is

unsafeReadElemBanded :: (ReadBanded a x e m) => a mn e -> (Int,Int) -> m e
unsafeReadElemBanded a (i,j)
    | isHerm a = 
        unsafeReadElemBanded (hermBanded $ coerceBanded a) (j,i) 
        >>= return . conj
    | hasStorageBanded a (i,j) =
        unsafeIOToM $
            withForeignPtr (fptrOfBanded a) $ \ptr ->
                peekElemOff ptr (indexOfBanded a (i,j))
    | otherwise =
        return 0
{-# INLINE unsafeReadElemBanded #-}


------------------------- WriteTensor functions -----------------------------

-- | Create a new banded matrix of given shape and (lower,upper), bandwidths,
-- but do not initialize the elements.
newBanded_ :: (WriteBanded a x e m) => (Int,Int) -> (Int,Int) -> m (a mn e)
newBanded_ (m,n) (kl,ku)
    | m < 0 || n < 0 =
        err "dimensions must be non-negative."
    | kl < 0 =
        err "lower bandwdth must be non-negative."
    | m /= 0 && kl >= m =
        err "lower bandwidth must be less than m."
    | ku < 0 =
        err "upper bandwidth must be non-negative."
    | n /= 0 && ku >= n =
        err "upper bandwidth must be less than n."
    | otherwise =
        let off = 0
            m'  = kl + 1 + ku
            l   = m'
            h   = False
        in unsafeIOToM $ do    
            ptr <- mallocForeignPtrArray (m' * n)
            return $ bandedViewArray ptr off (m,n) (kl,ku) l h
    where
      err s = fail $ "newBanded_ " ++ show (m,n) ++ " " ++ show (kl,ku) ++ ": " ++ s
                  
-- | Create a zero banded matrix of the specified shape and bandwidths.
newZeroBanded :: (WriteBanded a x e m) => (Int,Int) -> (Int,Int) -> m (a mn e)
newZeroBanded mn bw = do
    a <- newBanded_ mn bw
    setZeroBanded a
    return a

-- | Create a constant banded matrix of the specified shape and bandwidths.
newConstantBanded :: (WriteBanded a x e m) => (Int,Int) -> (Int,Int) -> e -> m (a mn e)
newConstantBanded mn bw e = do
    a <- newBanded_ mn bw
    setConstantBanded e a
    return a

setZeroBanded :: (WriteBanded a x e m) => a mn e -> m ()    
setZeroBanded = setConstantBanded 0

setConstantBanded :: (WriteBanded a x e m) => e -> a mn e -> m ()
setConstantBanded e a
    | isHerm a = setConstantBanded (conj e) a'
    | otherwise = do
        is <- getIndicesBanded a
        mapM_ (\i -> unsafeWriteElemBanded a i e) is
  where
    a' = (hermBanded . coerceBanded) a

unsafeWriteElemBanded :: (WriteBanded a x e m) => 
    a mn e -> (Int,Int) -> e -> m ()
unsafeWriteElemBanded a (i,j) e
    | isHerm a  = unsafeWriteElemBanded a' (j,i) $ conj e
    | otherwise = unsafeIOToM $
                      withForeignPtr (fptrOfBanded a) $ \ptr ->
                          pokeElemOff ptr (indexOfBanded a (i,j)) e
  where
    a' = (hermBanded . coerceBanded) a

modifyWithBanded :: (WriteBanded a x e m) => (e -> e) -> a mn e -> m ()
modifyWithBanded f a = do
    ies <- getAssocsBanded a
    mapM_ (\(ij,e) -> unsafeWriteElemBanded a ij (f e)) ies

canModifyElemBanded :: (WriteBanded a x e m) => a mn e -> (Int,Int) -> m Bool
canModifyElemBanded a ij = return $ hasStorageBanded a ij
{-# INLINE canModifyElemBanded #-}


--------------------------- Utility functions -------------------------------

indexOfBanded :: (BaseBanded a x e) => a mn e -> (Int,Int) -> Int
indexOfBanded a (i,j) =
    let (_,off,_,(_,ku),ld,h) = arrayFromBanded a
        (i',j')               = if h then (j,i) else (i,j)
    in off + ku + (i' - j') + j' * ld

hasStorageBanded :: (BaseBanded a x e) => a mn e -> (Int,Int) -> Bool
hasStorageBanded a (i,j) =
    let (_,_,(m,_),(kl,ku),_,h) = arrayFromBanded a
        (i',j')                 = if h then (j,i) else (i,j)
    in inRange (max 0 (j'-ku), min (m-1) (j'+kl)) i'

sizeBanded :: (BaseBanded a x e) => a mn e -> Int
sizeBanded a =
    let (_,_,(m,n),(kl,ku),_,_) = arrayFromBanded a
    in foldl' (+) 0 $ map (diagLen (m,n)) [(-kl)..ku]

indicesBanded :: (BaseBanded a x e) => a mn e -> [(Int,Int)]
indicesBanded a =
    let is = if isHerm a 
                 then [ (i,j) | i <- range (0,m-1), j <- range (0,n-1) ]
                 else [ (i,j) | j <- range (0,n-1), i <- range (0,m-1) ]
    in filter (hasStorageBanded a) is
  where (m,n) = shapeBanded a
