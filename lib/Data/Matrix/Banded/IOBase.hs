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

import BLAS.Internal( diagLen, diagStart )
import Data.Elem.BLAS( Elem, BLAS1, BLAS2, BLAS3, Trans(..), flipUpLo, 
    conjugate )
import qualified Data.Elem.BLAS as BLAS

import Data.Matrix.Class
import Data.Matrix.Class.MMatrixBase
import Data.Matrix.Class.MSolveBase

import Data.Matrix.Herm
import Data.Matrix.Tri

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Matrix.Dense.IOBase( IOMatrix(..) )
import Data.Matrix.Dense.Base( unsafeAxpyMatrix, unsafeCopyMatrix )
import Data.Matrix.Dense.Class
import Data.Vector.Dense.IOBase( IOVector(..), withIOVector )
import Data.Vector.Dense.Base( unsafeAxpyVector, unsafeCopyVector )
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

maybeBandedFromIOMatrix :: (Int,Int)
                   -> (Int,Int)
                   -> IOMatrix np e
                   -> Maybe (IOBanded np' e)
maybeBandedFromIOMatrix (m,n) (kl,ku) (IOMatrix f p m' n' l h)
    | h         = Nothing
    | otherwise =
         case undefined of
             _ | kl < 0 ->
                  error $ "maybeBandedFromMatrix:"
                        ++ " lower bandwidth must be non-negative"
             _ | ku < 0 ->
                  error $ "maybeBandedFromMatrix:"
                        ++ " upper bandwidth must be non-negative"
             _ | m' /= kl+1+ku ->
                  error $ "maybeBandedFromMatrix:"
                        ++ " number of rows must be equal to number of diagonals"
             _ | n' /= n ->
                  error $ "maybeBandedFromMatrix:"
                        ++ " numbers of columns must be equal"
             _ ->
                  Just $ IOBanded f p m n kl ku l False

bandwidthsIOBanded :: IOBanded np e -> (Int,Int)
bandwidthsIOBanded a =
    (numLowerIOBanded a, numUpperIOBanded a)
{-# INLINE bandwidthsIOBanded #-}

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

newCopyIOBanded :: (BLAS1 e) => IOBanded np e -> IO (IOBanded np e)
newCopyIOBanded a 
    | isHermIOBanded a =
        liftM hermIOBanded $ newCopyIOBanded (hermIOBanded a)
    | otherwise = do
        a' <- newIOBanded_ (shapeIOBanded a) (numLowerIOBanded a, numUpperIOBanded a)
        unsafeCopyIOBanded a' a
        return a'

unsafeCopyIOBanded :: (BLAS1 e)
                   => IOBanded mn e -> IOBanded mn e -> IO ()
unsafeCopyIOBanded dst src
    | isHermIOBanded dst = 
        unsafeCopyIOBanded (hermIOBanded dst) 
                           (hermIOBanded src)
                         
    | (not . isHermIOBanded) src =
        withIOBanded (coerceIOBanded dst) $ \pDst ->
        withIOBanded (coerceIOBanded src) $ \pSrc ->
            if ldDst == m && ldSrc == m
                then copyBlock pDst pSrc
                else copyCols  pDst pSrc n
                
    | otherwise =
        zipWithM_ unsafeCopyVector (diagViews dst) (diagViews src)
        
  where
    m     = numLowerIOBanded dst + numUpperIOBanded dst + 1 -- we can be sure dst is not herm
    n     = numColsIOBanded dst
    ldDst = ldaIOBanded dst
    ldSrc = ldaIOBanded src

    copyBlock pDst pSrc =
        BLAS.copy (m*n) pSrc 1 pDst 1

    copyCols pDst pSrc nleft
        | nleft == 0 = return ()
        | otherwise = do
            BLAS.copy m pSrc 1 pDst 1
            copyCols (pDst `advancePtr` ldDst) (pSrc `advancePtr` ldSrc) 
                     (nleft-1)

    diagViews a = 
        case bandwidthsIOBanded a of
            (kl,ku) -> map (unsafeDiagViewIOBanded a) $ range (-kl,ku)


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

unsafeDiagViewIOBanded :: (Elem e) => 
    IOBanded np e -> Int -> IOVector k e
unsafeDiagViewIOBanded a@(IOBanded fp p m n _ _ ld h) d
    | h         = conj $ unsafeDiagViewIOBanded (hermIOBanded a) (negate d)
    | otherwise =
        let off = indexIOBanded a (diagStart d)
            p'  = p `advancePtr` off
            len = diagLen (m,n) d
            inc = ld
            c   = False
        in (IOVector fp p' len inc c)


unsafeGetRowIOBanded :: (WriteVector y e IO, BLAS1 e) 
                     => IOBanded np e -> Int -> IO (y p e)
unsafeGetRowIOBanded a i =
    let (nb,x,na) = unsafeRowViewIOBanded a i
        n = numCols a
    in do
        es <- getElems x
        newListVector n $ (replicate nb 0) ++ es ++ (replicate na 0)

unsafeGetColIOBanded :: (WriteVector y e IO, BLAS1 e) 
                     => IOBanded np e -> Int -> IO (y n e)
unsafeGetColIOBanded a j =
    liftM conj $ unsafeGetRowIOBanded (hermIOBanded a) j

-- | @gbmv alpha a x beta y@ replaces @y := alpha a * x + beta y@
gbmv :: (ReadVector x e IO, WriteVector y e IO, BLAS2 e)
     => e -> IOBanded (n,p) e -> x p e -> e -> y n e -> IO ()
gbmv alpha a (x :: x p e) beta y
    | numRows a == 0 || numCols a == 0 =
        scaleByVector beta y
    | isConj x = do
        x' <- newCopyVector' x :: IO (IOVector p e) 
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
gbmm :: (ReadMatrix b e IO, WriteMatrix c e IO, BLAS2 e)
     => e -> IOBanded (n,p) e -> b (p,q) e -> e -> c (n,q) e -> IO ()
gbmm alpha a b beta c =
    sequence_ $
        zipWith (\x y -> gbmv alpha a x beta y) (colViews b) (colViews c)

hbmv :: (ReadVector x e IO, WriteVector y e IO, BLAS2 e) => 
    e -> Herm IOBanded (n,p) e -> x p e -> e -> y n e -> IO ()
hbmv alpha h (x :: x p e) beta (y :: y n e)
    | numRows h == 0 =
        return ()
    | isConj y = do
        doConjVector y
        hbmv alpha h x beta (conj y)
        doConjVector y
    | isConj x = do
        (x' :: IOVector p e) <- newCopyVector' x
        hbmv alpha h x' beta y
    | otherwise =
        let (u,a) = hermToBase (coerceHerm h)
            n     = numCols a
            k     = case u of 
                        Upper -> numUpperIOBanded a
                        Lower -> numLowerIOBanded a      
            u'    = case (isHermIOBanded a) of
                        True  -> flipUpLo u
                        False -> u
            uploA = u'
            ldA   = ldaIOBanded a
            incX  = stride x
            incY  = stride y
            withPtrA 
                  = case u' of Upper -> withIOBanded a
                               Lower -> withIOBandedElem a (0,0)
        in
            withPtrA $ \pA ->
            withIOVector (unsafeVectorToIOVector x) $ \pX ->
            withIOVector (unsafeVectorToIOVector y) $ \pY -> do
                BLAS.hbmv uploA n k alpha pA ldA pX incX beta pY incY

hbmm :: (ReadMatrix b e IO, WriteMatrix c e IO) => 
    e -> Herm IOBanded (n,p) e -> b (p,q) e -> e -> c (n,q) e -> IO ()
hbmm alpha h b beta c =
    zipWithM_ (\x y -> hbmv alpha h x beta y) (colViews b) (colViews c)

tbmv :: (WriteVector y e IO, BLAS2 e) => 
    e -> Tri IOBanded (n,n) e -> y n e -> IO ()
tbmv alpha t x | isConj x = do
    doConjVector x
    tbmv alpha t (conj x)
    doConjVector x

tbmv alpha t x =
    let (u,d,a) = triToBase t
        (transA,u') 
                  = if isHermIOBanded a then (ConjTrans, flipUpLo u) 
                                        else (NoTrans  , u)
        uploA     = u'
        diagA     = d
        n         = numCols a
        k         = case u of Upper -> numUpperIOBanded a 
                              Lower -> numLowerIOBanded a
        ldA       = ldaIOBanded a
        incX      = stride x
        withPtrA  = case u' of 
                        Upper -> withIOBanded a
                        Lower -> withIOBandedElem a (0,0)
    in do
        scaleByVector alpha x
        withPtrA $ \pA ->
            withVectorPtrIO x $ \pX -> do
                BLAS.tbmv uploA transA diagA n k pA ldA pX incX
  where withVectorPtrIO = withIOVector . unsafeVectorToIOVector

tbmm :: (WriteMatrix b e IO, BLAS2 e) =>
    e -> Tri IOBanded (n,n) e -> b (n,p) e -> IO ()
tbmm 1     t b = mapM_ (\x -> tbmv 1 t x) (colViews b)
tbmm alpha t b = scaleByMatrix alpha b >> tbmm 1 t b

tbmv' :: (ReadVector x e IO, WriteVector y e IO, BLAS2 e) => 
    e -> Tri IOBanded (n,p) e -> x p e -> e -> y n e -> IO ()
tbmv' alpha a (x :: x p e) beta (y  :: y n e)
    | beta /= 0 = do
        (x' :: y p e) <- newCopyVector x
        tbmv alpha (coerceTri a) x'
        scaleByVector beta y
        unsafeAxpyVector 1 x' (coerceVector y)
    | otherwise = do
        unsafeCopyVector (coerceVector y) x
        tbmv alpha (coerceTri a) (coerceVector y)

tbmm' :: (ReadMatrix b e IO, WriteMatrix c e IO, BLAS2 e) => 
    e -> Tri IOBanded (r,s) e -> b (s,t) e -> e -> c (r,t) e -> IO ()
tbmm' alpha a (b :: b (s,t) e) beta (c :: c (r,t) e)
    | beta /= 0 = do
        (b' :: c (s,t) e) <- newCopyMatrix b
        tbmm alpha (coerceTri a) b'
        scaleByMatrix beta c
        unsafeAxpyMatrix 1 b' (coerceMatrix c)
    | otherwise = do
        unsafeCopyMatrix (coerceMatrix c) b
        tbmm alpha (coerceTri a) (coerceMatrix c)

tbsv :: (WriteVector y e IO, BLAS2 e) => 
    e -> Tri IOBanded (k,k) e -> y n e -> IO ()
tbsv alpha t x | isConj x = do
    doConjVector x
    tbsv alpha t (conj x)
    doConjVector x
tbsv alpha t x = 
    let (u,d,a) = triToBase t
        (transA,u') = if isHermIOBanded a then (ConjTrans, flipUpLo u)
                                          else (NoTrans  , u)
        uploA     = u'
        diagA     = d
        n         = numCols a
        k         = case u of Upper -> numUpperIOBanded a 
                              Lower -> numLowerIOBanded a        
        ldA       = ldaIOBanded a
        incX      = stride x
        withPtrA  = case u' of 
                        Upper -> withIOBanded a
                        Lower -> withIOBandedElem a (0,0)
    in do
        scaleByVector alpha x
        withPtrA $ \pA ->
            withVectorPtrIO x $ \pX -> do
                BLAS.tbsv uploA transA diagA n k pA ldA pX incX
  where withVectorPtrIO = withIOVector . unsafeVectorToIOVector

tbsm :: (WriteMatrix b e IO, BLAS2 e) => 
    e -> Tri IOBanded (k,k) e -> b (k,l) e -> IO ()
tbsm 1     t b = mapM_ (\x -> tbsv 1 t x) (colViews b)
tbsm alpha t b = scaleByMatrix alpha b >> tbsm 1 t b

tbsv' :: (ReadVector y e IO, WriteVector x e IO, BLAS2 e)
      => e -> Tri IOBanded (k,l) e -> y k e -> x l e -> IO ()
tbsv' alpha a y x = do
    unsafeCopyVector (coerceVector x) y
    tbsv alpha (coerceTri a) (coerceVector x)

tbsm' :: (ReadMatrix c e IO, WriteMatrix b e IO, BLAS2 e) 
      => e -> Tri IOBanded (r,s) e -> c (r,t) e -> b (s,t) e -> IO ()
tbsm' alpha a c b = do
    unsafeCopyMatrix (coerceMatrix b) c
    tbsm alpha (coerceTri a) b

instance HasVectorView IOBanded where
    type VectorView IOBanded = IOVector

instance HasMatrixStorage IOBanded where
    type MatrixStorage IOBanded = IOMatrix

instance (Elem e) => Shaped IOBanded (Int,Int) e where
    shape  = shapeIOBanded
    {-# INLINE shape #-}
    bounds = boundsIOBanded
    {-# INLINE bounds #-}

instance (Elem e) => MatrixShaped IOBanded e where
    herm = hermIOBanded
    {-# INLINE herm #-}

instance (BLAS3 e) => ReadTensor IOBanded (Int,Int) e IO where
    getSize        = getSizeIOBanded
    {-# INLINE getSize #-}
    getAssocs      = getAssocsIOBanded
    {-# INLINE getAssocs #-}
    getIndices     = getIndicesIOBanded
    {-# INLINE getIndices #-}
    getElems       = getElemsIOBanded
    {-# INLINE getElems #-}
    getAssocs'     = getAssocsIOBanded'
    {-# INLINE getAssocs' #-}
    getIndices'    = getIndicesIOBanded'
    {-# INLINE getIndices' #-}
    getElems'      = getElemsIOBanded'
    {-# INLINE getElems' #-}
    unsafeReadElem = unsafeReadElemIOBanded
    {-# INLINE unsafeReadElem #-}

instance (BLAS3 e) => WriteTensor IOBanded (Int,Int) e IO where
    setConstant     = setConstantIOBanded
    {-# INLINE setConstant #-}
    setZero         = setZeroIOBanded
    {-# INLINE setZero #-}
    modifyWith      = modifyWithIOBanded
    {-# INLINE modifyWith #-}
    unsafeWriteElem = unsafeWriteElemIOBanded
    {-# INLINE unsafeWriteElem #-}
    canModifyElem   = canModifyElemIOBanded
    {-# INLINE canModifyElem #-}

instance (BLAS3 e) => MMatrix IOBanded e IO where
    unsafeDoSApplyAdd    = gbmv
    unsafeDoSApplyAddMat = gbmm
    unsafeGetRow         = unsafeGetRowIOBanded
    unsafeGetCol         = unsafeGetColIOBanded

instance (BLAS3 e) => MMatrix (Herm IOBanded) e IO where
    unsafeDoSApplyAdd    = hbmv
    unsafeDoSApplyAddMat = hbmm

instance (BLAS3 e) => MMatrix (Tri IOBanded) e IO where
    unsafeDoSApply_      = tbmv
    unsafeDoSApplyMat_   = tbmm
    unsafeDoSApplyAdd    = tbmv'
    unsafeDoSApplyAddMat = tbmm'

instance (BLAS3 e) => MSolve (Tri IOBanded) e IO where
    unsafeDoSSolve     = tbsv'
    unsafeDoSSolveMat  = tbsm'
    unsafeDoSSolve_    = tbsv
    unsafeDoSSolveMat_ = tbsm




transIOBanded :: IOBanded np e -> Trans
transIOBanded a =
    case (isHermIOBanded a) of
          False -> NoTrans
          True  -> ConjTrans
{-# INLINE transIOBanded #-}

