{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances,
        TypeFamilies, ScopedTypeVariables #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.IOBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.IOBase
    where

import Data.Ix
import Data.List
import Control.Monad
import Foreign
import System.IO.Unsafe
import Unsafe.Coerce

import BLAS.Internal( diagLen )
import Data.Elem.BLAS( Elem, BLAS1, BLAS2, BLAS3, Trans(..), conjugate )
import qualified Data.Elem.BLAS as BLAS

import Data.Matrix.Class
import Data.Matrix.Class.MMatrixBase
import Data.Matrix.Class.MSolveBase

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Matrix.Dense.IOBase( IOMatrix(..) )
import Data.Matrix.Dense.Class
import Data.Vector.Dense.IOBase( IOVector(..), withIOVector )
import Data.Vector.Dense.Class

-- | The Banded matrix data type.
data IOBanded np e =
    IOBanded { fptrIOBanded     :: {-# UNPACK #-} !(ForeignPtr e)
             , ptrIOBanded      :: {-# UNPACK #-} !(Ptr e)
             , numRowsIOBanded  :: {-# UNPACK #-} !Int
             , numColsIOBanded  :: {-# UNPACK #-} !Int
             , numLowerIOBanded :: {-# UNPACK #-} !Int
             , numUpperIOBanded :: {-# UNPACK #-} !Int
             , ldaIOBanded      :: {-# UNPACK #-} !Int
             , isHermIOBanded   :: {-# UNPACK #-} !Bool
             }

hermIOBanded :: IOBanded np e -> IOBanded pn e
hermIOBanded a = a{ numRowsIOBanded  = numColsIOBanded a
                  , numColsIOBanded  = numRowsIOBanded a
                  , numLowerIOBanded = numUpperIOBanded a
                  , numUpperIOBanded = numLowerIOBanded a
                  , isHermIOBanded   = not (isHermIOBanded a)
                  }
{-# INLINE hermIOBanded #-}

indexIOBanded :: IOBanded np e -> (Int,Int) -> Int
indexIOBanded (IOBanded _ _ _ _ kl ku ld h) (i,j) =
    let (ku',i',j') = if h then (kl,j,i) else (ku,i,j)
    in ku + (i' - j') + j' * ld
{-# INLINE indexIOBanded #-}

hasStorageIOBanded :: IOBanded np e -> (Int,Int) -> Bool
hasStorageIOBanded (IOBanded _ _ m n kl ku _ h) (i,j) =
    let (m',kl',ku',i',j') = if h then (n,ku,kl,j,i) else (m,kl,ku,i,j)
    in inRange (max 0 (j'-ku), min (m-1) (j'+kl)) i'
{-# INLINE hasStorageIOBanded #-}

withIOBanded :: IOBanded (n,p) e -> (Ptr e -> IO a) -> IO a
withIOBanded b f = do
    a <- f (ptrIOBanded b)
    touchForeignPtr (fptrIOBanded b)
    return a
{-# INLINE withIOBanded #-}

withIOBandedElem :: (Elem e)
                 => IOBanded np e -> (Int,Int) -> (Ptr e -> IO a) -> IO a
withIOBandedElem a (i,j) f =
    withIOBanded (coerceIOBanded a) $ \p ->
        f $ p `advancePtr` (indexIOBanded a (i,j))
{-# INLINE withIOBandedElem #-}

matrixIOBanded :: IOBanded np e -> IOMatrix np' e
matrixIOBanded (IOBanded f p m n kl ku l h) =
    let n' = if h then m else n
    in IOMatrix f p (kl+1+ku) n' l False
{-# INLINE matrixIOBanded #-}

bandedFromIOMatrix :: (Int,Int)
                   -> (Int,Int)
                   -> IOMatrix np e
                   -> Maybe (IOBanded np' e)
bandedFromIOMatrix (m,n) (kl,ku) (IOMatrix f p m' n' l h)
    | h         = Nothing
    | otherwise =
         case undefined of
             _ | m' /= kl+1+ku ->
                  error $ "bandedFromMatrix:"
                        ++ " number of rows must be equal to number of diagonals"
             _ | n' /= n ->
                  error $ "bandedFromMatrix:"
                        ++ " numbers of columns must be equal"
             _ ->
                  Just $ IOBanded f p m n kl ku l False

bandwidthIOBanded :: IOBanded np e -> (Int,Int)
bandwidthIOBanded a =
    let (kl,ku) = (numLowerIOBanded a, numUpperIOBanded a)
    in (negate kl, ku)
{-# INLINE bandwidthIOBanded #-}

coerceIOBanded :: IOBanded np e -> IOBanded np' e
coerceIOBanded = unsafeCoerce
{-# INLINE coerceIOBanded #-}

shapeIOBanded :: IOBanded np e -> (Int,Int)
shapeIOBanded a = (numRowsIOBanded a, numColsIOBanded a)
{-# INLINE shapeIOBanded #-}

boundsIOBanded :: IOBanded np e -> ((Int,Int), (Int,Int))
boundsIOBanded a = ((0,0), (m-1,n-1)) where (m,n) = shapeIOBanded a
{-# INLINE boundsIOBanded #-}

sizeIOBanded :: IOBanded np e -> Int
sizeIOBanded (IOBanded _ _ m n kl ku _ _) =
    foldl' (+) 0 $ map (diagLen (m,n)) [(-kl)..ku]

indicesIOBanded :: IOBanded np e -> [(Int,Int)]
indicesIOBanded a =
    let is = if isHermIOBanded a
                 then [ (i,j) | i <- range (0,m-1), j <- range (0,n-1) ]
                 else [ (i,j) | j <- range (0,n-1), i <- range (0,m-1) ]
    in filter (hasStorageIOBanded a) is
  where (m,n) = shapeIOBanded a

getSizeIOBanded :: IOBanded np e -> IO Int
getSizeIOBanded = return . sizeIOBanded
{-# INLINE getSizeIOBanded #-}

getIndicesIOBanded :: IOBanded np e -> IO [(Int,Int)]
getIndicesIOBanded = return . indicesIOBanded
{-# INLINE getIndicesIOBanded #-}

getIndicesIOBanded' :: IOBanded np e -> IO [(Int,Int)]
getIndicesIOBanded' = getIndicesIOBanded
{-# INLINE getIndicesIOBanded' #-}

getAssocsIOBanded :: (Elem e) => IOBanded np e -> IO [((Int,Int),e)]
getAssocsIOBanded a =
    let go []     = return []
        go (i:is) = unsafeInterleaveIO $ do
                        e   <- unsafeReadElemIOBanded a i
                        ies <- go is
                        return $ (i,e):ies
    in go (indicesIOBanded a)
{-# INLINE getAssocsIOBanded #-}

getAssocsIOBanded' :: (Elem e) => IOBanded np e -> IO [((Int,Int),e)]
getAssocsIOBanded' a =
    sequence [ liftM ((,) i) (unsafeReadElemIOBanded a i)
             | i <- indicesIOBanded a ]
{-# INLINE getAssocsIOBanded' #-}

getElemsIOBanded :: (Elem e) => IOBanded np e -> IO [e]
getElemsIOBanded a = liftM (snd . unzip) $ getAssocsIOBanded a
{-# INLINE getElemsIOBanded #-}

getElemsIOBanded' :: (Elem e) => IOBanded np e -> IO [e]
getElemsIOBanded' a = liftM (snd . unzip) $ getAssocsIOBanded a
{-# INLINE getElemsIOBanded' #-}

unsafeReadElemIOBanded :: (Elem e) => IOBanded np e -> (Int,Int) -> IO e
unsafeReadElemIOBanded a (i,j)
    | hasStorageIOBanded a (i,j) =
        let f = if isHermIOBanded a then conjugate else id
        in liftM f $ withIOBandedElem a (i,j) peek
    | otherwise =
        return 0
{-# INLINE unsafeReadElemIOBanded #-}

unsafeWriteElemIOBanded :: (Elem e)
                      => IOBanded np e -> (Int,Int) -> e -> IO ()
unsafeWriteElemIOBanded a ij e =
    let e' = if isHermIOBanded a then conjugate e else e
    in withIOBandedElem a ij (`poke` e')
{-# INLINE unsafeWriteElemIOBanded #-}

modifyWithIOBanded :: (Elem e)
                   => (e -> e) -> IOBanded np e -> IO ()
modifyWithIOBanded f a = do
    ies <- getAssocsIOBanded a
    mapM_ (\(i,e) -> unsafeWriteElemIOBanded a i (f e)) ies

canModifyElemIOBanded :: IOBanded np e -> (Int,Int) -> IO Bool
canModifyElemIOBanded a ij = return $ hasStorageIOBanded a ij
{-# INLINE canModifyElemIOBanded #-}

-- | Create a new banded matrix of given shape and (lower,upper), bandwidths,
-- but do not initialize the elements.
newIOBanded_ :: (Elem e) => (Int,Int) -> (Int,Int) -> IO (IOBanded np e)
newIOBanded_ (m,n) (kl,ku)
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
        let m'  = kl + 1 + ku
            l   = m'
            h   = False
        in do
            fp <- mallocForeignPtrArray (m' * n)
            let p = unsafeForeignPtrToPtr fp
            return $ IOBanded fp p m n kl ku l h
    where
      err s = fail $ "newBanded_ " ++ show (m,n) ++ " " ++ show (kl,ku) ++ ": " ++ s


-- | Create a zero banded matrix of the specified shape and bandwidths.
newZeroIOBanded :: (Elem e) => (Int,Int) -> (Int,Int) -> IO (IOBanded np e)
newZeroIOBanded mn bw = do
    a <- newIOBanded_ mn bw
    setZeroIOBanded a
    return a

setZeroIOBanded :: (Elem e) => IOBanded np e -> IO ()
setZeroIOBanded = setConstantIOBanded 0
{-# INLINE setZeroIOBanded #-}

-- | Create a constant banded matrix of the specified shape and bandwidths.
newConstantIOBanded :: (Elem e)
                    => (Int,Int) -> (Int,Int) -> e -> IO (IOBanded np e)
newConstantIOBanded mn bw e = do
    a <- newIOBanded_ mn bw
    setConstantIOBanded e a
    return a

setConstantIOBanded :: (Elem e) => e -> IOBanded np e -> IO ()
setConstantIOBanded e a
    | isHermIOBanded a = setConstantIOBanded (conjugate e) (hermIOBanded a)
    | otherwise = do
        is <- getIndicesIOBanded a
        mapM_ (\i -> unsafeWriteElemIOBanded a i e) is
{-# INLINE setConstantIOBanded #-}

unsafeRowViewIOBanded :: (Elem e)
                      => IOBanded np e
                      -> Int
                      -> (Int, IOVector k e, Int)
unsafeRowViewIOBanded a@(IOBanded f p _ n kl ku ld h) i =
    if h then
        case unsafeColViewIOBanded (hermIOBanded a) i of
            (nb, v, na) -> (nb, conj v, na)
    else
        let nb  = min n $ max (i - kl)         0
            na  = min n $ max (n - 1 - i - ku) 0
            r   = min (ku + i)         (kl + ku)
            c   = max (i - kl)         0
            p'  = p `advancePtr` (r + c * ld)
            inc = max (ld - 1) 1
            len = n - (nb + na)
        in (nb, IOVector f p' len inc False, na)

unsafeColViewIOBanded :: (Elem e)
                      => IOBanded np e
                      -> Int
                      -> (Int, IOVector k e, Int)
unsafeColViewIOBanded a@(IOBanded f p m _ kl ku ld h) j =
    if h then
        case unsafeRowViewIOBanded (hermIOBanded a) j of
            (nb, v, na) -> (nb, conj v, na)
    else
        let nb  = min m $ max (j - ku)         0
            na  = min m $ max (m - 1 - j - kl) 0
            r   = max (ku - j) 0
            c   = j
            p'  = p `advancePtr` (r + c * ld)
            inc = 1
            len = m - (nb + na)
        in (nb, IOVector f p' len inc False, na)

unsafeGetRowIOBanded :: (BLAS1 e) => IOBanded np e -> Int -> IO (IOVector p e)
unsafeGetRowIOBanded a i =
    let (nb,x,na) = unsafeRowViewIOBanded a i
        n = numCols a
    in do
        es <- getElems x
        newListVector n $ (replicate nb 0) ++ es ++ (replicate na 0)

unsafeGetColIOBanded :: (BLAS1 e) => IOBanded np e -> Int -> IO (IOVector n e)
unsafeGetColIOBanded a j =
    liftM conj $ unsafeGetRowIOBanded (hermIOBanded a) j

-- | @gbmv alpha a x beta y@ replaces @y := alpha a * x + beta y@
gbmv :: (ReadVector x e IO, ReadVector y e IO, BLAS2 e)
     => e -> IOBanded (n,p) e -> x p e -> e -> y n e -> IO ()
gbmv alpha a x beta y =
    gbmvIO alpha a (unsafeVectorToIOVector x) beta (unsafeVectorToIOVector y)
gbmv alpha a x beta y
    | numRows a == 0 || numCols a == 0 =
        scaleByVector beta y
    | isConj x = do
        newCopyVector' x
        gbmv alpha a x' beta y
    | isConj y = do
        doConjVector y
        gbmv alpha a x beta (conj y)
        doConjVector y
    | otherwise =
        let transA = transIOBanded a
            (m,n)  = case (isHermIOBanded a) of
                         False -> shape a
                         True  -> (flipShape . shape) a
            (kl,ku) = case (isHermIOBanded a) of
                          False -> (numLowerIOBanded a, numUpperIOBanded a)
                          True  -> (numUpperIOBanded a, numLowerIOBanded a)
            ldA    = ldaIOBanded a
            incX   = stride x
            incY   = stride y
        in
            withIOBanded a   $ \pA ->
            withIOVector (unsafeVectorToIOVector x) $ \pX ->
            withIOVector (unsafeVectorToIOVector y) $ \pY -> do
                BLAS.gbmv transA m n kl ku alpha pA ldA pX incX beta pY incY

-- | @gbmm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gbmm :: (BLAS2 e)
     => e -> IOBanded (n,p) e -> IOMatrix (p,q) e -> e -> IOMatrix (n,q) e -> IO ()
gbmm alpha a b beta c =
    sequence_ $
        zipWith (\x y -> gbmv alpha a x beta y) (colViews b) (colViews c)



instance HasVectorView IOBanded where
    type VectorView IOBanded = IOVector

instance (Elem e) => Shaped IOBanded (Int,Int) e where
    shape  = shapeIOBanded
    bounds = boundsIOBanded

instance (Elem e) => MatrixShaped IOBanded e where
    herm = hermIOBanded

instance (BLAS3 e) => ReadTensor IOBanded (Int,Int) e IO where
    getSize        = getSizeIOBanded
    getAssocs      = getAssocsIOBanded
    getIndices     = getIndicesIOBanded
    getElems       = getElemsIOBanded
    getAssocs'     = getAssocsIOBanded'
    getIndices'    = getIndicesIOBanded'
    getElems'      = getElemsIOBanded'
    unsafeReadElem = unsafeReadElemIOBanded

instance (BLAS3 e) => WriteTensor IOBanded (Int,Int) e IO where
    setConstant     = setConstantIOBanded
    setZero         = setZeroIOBanded
    modifyWith      = modifyWithIOBanded
    unsafeWriteElem = unsafeWriteElemIOBanded
    canModifyElem   = canModifyElemIOBanded

instance (BLAS3 e) => MMatrix IOBanded e IO where
    unsafeDoSApplyAdd    = gbmv
    unsafeDoSApplyAddMat = gbmm
    unsafeGetRow         = unsafeGetRowIOBanded
    unsafeGetCol         = unsafeGetColIOBanded

transIOBanded :: IOBanded np e -> Trans
transIOBanded a =
    case (isHermIOBanded a) of
          False -> NoTrans
          True  -> ConjTrans
{-# INLINE transIOBanded #-}

