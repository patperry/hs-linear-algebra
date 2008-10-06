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
    -- * Banded matrix types
    IOBanded,
    STBanded,
    unsafeIOBandedToSTBanded,
    unsafeSTBandedToIOBanded,

    -- * Banded type classes
    BaseBanded(..),
    ReadBanded,
    WriteBanded,

    -- * Low-level Banded properties
    bandedViewMatrix,
    matrixFromBanded,
    lda,
    isHerm,
    hermBanded,

    -- * Bandwidth properties
    bandwidth,
    numLower,
    numUpper,

    -- * Coercing the Banded shape
    coerceBanded,

    -- * WriteTensor functions
    newBanded_,
    newZeroBanded,
    setZeroBanded,
    newConstantBanded,
    setConstantBanded,
    modifyWithBanded,
    canModifyElemBanded,
    unsafeWriteElemBanded,

    -- * Vector views
    unsafeBandRowView,
    unsafeBandColView,
    unsafeGetRowBanded,
    unsafeGetColBanded,
    
    -- * Utility functions
    withBandedPtr,
    withBandedElemPtr,
    fptrOfBanded,
    offsetOfBanded,
    indexOfBanded,
    indicesBanded,
    
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
import BLAS.Internal( diagLen )
import BLAS.UnsafeIOToM
import BLAS.UnsafeInterleaveM

import BLAS.Matrix.Base hiding ( BaseMatrix )
import qualified BLAS.Matrix.Base as BLAS
import BLAS.Matrix.Mutable

import BLAS.Tensor

import Data.Vector.Dense.Class.Internal( IOVector, STVector,
    BaseVector(..), ReadVector, WriteVector, doConjVector,
    withVectorPtr, stride, isConj )
import Data.Vector.Dense.Class.Creating( newListVector )
import Data.Vector.Dense.Class.Operations( getConjVector )

import qualified Data.Matrix.Dense.Class.Internal as M
import Data.Matrix.Dense.Class( BaseMatrix, ReadMatrix, WriteMatrix,
    arrayFromMatrix, matrixViewArray, colViews )


class (BLAS.BaseMatrix a e, BaseVector x e) => 
    BaseBanded a x e | a -> x where
        bandedViewArray :: ForeignPtr e -> Int -> (Int,Int) -> (Int,Int) -> Int -> Bool -> a mn e
        arrayFromBanded :: a mn e -> (ForeignPtr e, Int, (Int,Int), (Int,Int), Int, Bool)

class (Elem e, UnsafeInterleaveM m, ReadTensor a (Int,Int) e m, 
           BaseBanded a x e, ReadVector x e m) => 
    ReadBanded a x e m | a -> x where

class (WriteTensor a (Int,Int) e m,
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


------------------------------ Vector views ---------------------------------

unsafeBandRowView :: (BaseBanded a x e) => a mn e -> Int -> (Int, x k e, Int)
unsafeBandRowView a i =
    if h then
        case unsafeBandColView a' i of (nb, v, na) -> (nb, conj v, na)        
    else
        let nb   = max (i - kl)         0
            na   = max (n - 1 - i - ku) 0
            r    = min (ku + i)         (kl + ku)
            c    = max (i - kl)         0 
            off' = off + r + c * ld
            inc  = ld - 1
            len  = n - (nb + na)
        in if len >= 0 
            then (nb, vectorViewArray f off' len inc False, na)
            else (n , vectorViewArray f off' 0   inc False,  0)
  where
    (f,off,(_,n),(kl,ku),ld,h) = arrayFromBanded a
    a' = (hermBanded . coerceBanded) a

unsafeBandColView :: (BaseBanded a x e) => a mn e -> Int -> (Int, x k e, Int)
unsafeBandColView a j =
    if h then
        case unsafeBandRowView a' j of (nb, v, na) -> (nb, conj v, na)
    else
        let nb    = max (j - ku)         0
            na    = max (m - 1 - j - kl) 0
            r     = max (ku - j) 0 
            c     = j 
            off'  = off + r + c * ld
            inc   = 1
            len   = m - (nb + na)
        in if len >= 0
            then (nb, vectorViewArray f off' len inc False, na)
            else (m , vectorViewArray f off' 0   inc False,  0)
  where
    (f,off,(m,_),(kl,ku),ld,h) = arrayFromBanded a
    a' = (hermBanded . coerceBanded) a

unsafeGetRowBanded :: (ReadBanded a x e m, WriteVector y e m) => 
    a (k,l) e -> Int -> m (y l e)
unsafeGetRowBanded a i = 
    let (nb,x,na) = unsafeBandRowView a i
        n = numCols a
    in do
        es <- getElems x
        newListVector n $ (replicate nb 0) ++ es ++ (replicate na 0)

unsafeGetColBanded :: (ReadBanded a x e m, WriteVector y e m) => 
    a (k,l) e -> Int -> m (y k e)
unsafeGetColBanded a j = unsafeGetRowBanded (hermBanded a) j >>= return . conj


-------------------------- Matrix multiplication ----------------------------

-- | @gbmv alpha a x beta y@ replaces @y := alpha a * x + beta y@
gbmv :: (ReadBanded a z e m, ReadVector x e m, WriteVector y e m, BLAS2 e) => 
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
gbmv alpha a x beta y
    | isConj x = do
        x' <- getConjVector (conj x)
        gbmv alpha a x' beta y
    | isConj y = do
        doConjVector y
        gbmv alpha a x beta (conj y)
        doConjVector y
    | otherwise =
        let order  = colMajor
            transA = blasTransOf a
            (m,n)  = case (isHerm a) of
                         False -> shape a
                         True  -> (flipShape . shape) a
            (kl,ku) = case (isHerm a) of
                          False -> (numLower a, numUpper a)
                          True  -> (numUpper a, numLower a)
            ldA    = lda a
            incX   = stride x
            incY   = stride y
        in unsafeIOToM $
               withBandedPtr a $ \pA ->
               withVectorPtr x $ \pX ->
               withVectorPtr y $ \pY -> do
                   BLAS.gbmv order transA m n kl ku alpha pA ldA pX incX beta pY incY

-- | @gbmm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gbmm :: (ReadBanded a x e m, ReadMatrix b y e m, WriteMatrix c z e m, BLAS2 e) => 
    e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
gbmm alpha a b beta c =
    sequence_ $
        zipWith (\x y -> gbmv alpha a x beta y) (colViews b) (colViews c)

--------------------------- Utility functions -------------------------------

withBandedPtr :: (BaseBanded a x e, Storable e) => 
    a mn e -> (Ptr e -> IO b) -> IO b
withBandedPtr a f =
    withForeignPtr (fptrOfBanded a) $ \ptr ->
        f $ ptr `advancePtr` (offsetOfBanded a)

withBandedElemPtr :: (BaseBanded a x e, Storable e) => 
    a mn e -> (Int,Int) -> (Ptr e -> IO b) -> IO b
withBandedElemPtr a (i,j) f
    | isHerm a  = withBandedElemPtr (hermBanded $ coerceBanded a) (j,i) f
    | otherwise = withForeignPtr (fptrOfBanded a) $ \ptr ->
                      f $ ptr `advancePtr` (indexOfBanded a (i,j))

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

blasTransOf :: (BaseBanded a x e) => a mn e -> CBLASTrans
blasTransOf a = 
    case (isHerm a) of
          False -> noTrans
          True  -> conjTrans

flipShape :: (Int,Int) -> (Int,Int)
flipShape (m,n) = (n,m)


------------------------------------ Instances ------------------------------

-- | The Banded matrix data type.
data IOBanded mn e =
    BM {-# UNPACK #-} !(ForeignPtr e) -- ^ base storage
       {-# UNPACK #-} !Int            -- ^ offset into storage
       {-# UNPACK #-} !Int            -- ^ numer of rows
       {-# UNPACK #-} !Int            -- ^ number of columns
       {-# UNPACK #-} !Int            -- ^ lower bandwidth
       {-# UNPACK #-} !Int            -- ^ upper bandwidth
       {-# UNPACK #-} !Int            -- ^ lda of storage
       {-# UNPACK #-} !Bool           -- ^ isHerm flag

newtype STBanded s mn e = ST (IOBanded mn e)

unsafeIOBandedToSTBanded :: IOBanded mn e -> STBanded s mn e
unsafeIOBandedToSTBanded = ST

unsafeSTBandedToIOBanded :: STBanded s mn e -> IOBanded mn e
unsafeSTBandedToIOBanded (ST x) = x

instance (Elem e) => BaseBanded IOBanded IOVector e where
    bandedViewArray f o (m,n) (kl,ku) ld h = BM f o m n kl ku ld h
    arrayFromBanded (BM f o m n kl ku ld h) = (f,o,(m,n),(kl,ku),ld,h)

instance (Elem e) => BaseBanded (STBanded s) (STVector s) e where
    bandedViewArray f o (m,n) (kl,ku) ld h       = ST (BM f o m n kl ku ld h)
    arrayFromBanded (ST (BM f o m n kl ku ld h)) = (f,o,(m,n),(kl,ku),ld,h)

instance (Elem e) => BaseTensor IOBanded (Int,Int) e where
    shape  = shapeBanded
    bounds = boundsBanded
    
instance (Elem e) => BaseTensor (STBanded s) (Int,Int) e where
    shape  = shapeBanded
    bounds = boundsBanded

instance (Elem e) => BLAS.BaseMatrix IOBanded e where
    herm = hermBanded
    
instance (Elem e) => BLAS.BaseMatrix (STBanded s) e where
    herm = hermBanded

instance (BLAS1 e) => ReadBanded IOBanded     IOVector     e IO
instance (BLAS1 e) => ReadBanded (STBanded s) (STVector s) e (ST s)

instance (BLAS1 e) => ReadTensor IOBanded (Int,Int) e IO where
    getSize        = getSizeBanded
    getAssocs      = getAssocsBanded
    getIndices     = getIndicesBanded
    getElems       = getElemsBanded
    getAssocs'     = getAssocsBanded'
    getIndices'    = getIndicesBanded'
    getElems'      = getElemsBanded'
    unsafeReadElem = unsafeReadElemBanded
    
instance (BLAS1 e) => ReadTensor (STBanded s) (Int,Int) e (ST s) where
    getSize        = getSizeBanded
    getAssocs      = getAssocsBanded
    getIndices     = getIndicesBanded
    getElems       = getElemsBanded
    getAssocs'     = getAssocsBanded'
    getIndices'    = getIndicesBanded'
    getElems'      = getElemsBanded'
    unsafeReadElem = unsafeReadElemBanded

instance (BLAS1 e) => WriteBanded IOBanded IOVector e IO where
instance (BLAS1 e) => WriteBanded (STBanded s) (STVector s) e (ST s) where

instance (BLAS1 e) => WriteTensor IOBanded (Int,Int) e IO where
    setConstant     = setConstantBanded
    setZero         = setZeroBanded
    modifyWith      = modifyWithBanded
    unsafeWriteElem = unsafeWriteElemBanded
    canModifyElem   = canModifyElemBanded

instance (BLAS1 e) => WriteTensor (STBanded s) (Int,Int) e (ST s) where
    setConstant     = setConstantBanded
    setZero         = setZeroBanded
    modifyWith      = modifyWithBanded
    unsafeWriteElem = unsafeWriteElemBanded
    canModifyElem   = canModifyElemBanded

instance (BLAS2 e) => MMatrix IOBanded e IO where
    unsafeDoSApplyAdd    = gbmv
    unsafeDoSApplyAddMat = gbmm
    unsafeGetRow         = unsafeGetRowBanded
    unsafeGetCol         = unsafeGetColBanded

instance (BLAS2 e) => MMatrix (STBanded s) e (ST s) where
    unsafeDoSApplyAdd    = gbmv
    unsafeDoSApplyAddMat = gbmm
    unsafeGetRow         = unsafeGetRowBanded
    unsafeGetCol         = unsafeGetColBanded
