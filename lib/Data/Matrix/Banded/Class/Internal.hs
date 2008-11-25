{-# LANGUAGE FlexibleInstances, FlexibleContexts, MultiParamTypeClasses
    , FunctionalDependencies, TypeFamilies #-}
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
    BaseBanded_(..),
    BaseBanded,
    ReadBanded,
    WriteBanded,

    -- * Low-level Banded properties
    bandedViewMatrix,
    matrixFromBanded,
    ldaOfBanded,
    isHermBanded,
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
    unsafeRowViewBanded,
    unsafeColViewBanded,
    unsafeGetRowBanded,
    unsafeGetColBanded,
    
    -- * Utility functions
    shapeBanded,
    boundsBanded,
    withBandedPtr,
    withBandedElemPtr,
    indexOfBanded,
    indicesBanded,
    gbmv,
    gbmm,
    hbmv',
    hbmm',
    tbmv,
    tbmm,
    tbmv',
    tbmm',
    unsafeDoSSolveTriBanded,
    unsafeDoSSolveMatTriBanded,
    tbsv,
    tbsm,
    
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

import BLAS.Matrix.Shaped
import BLAS.Matrix.Mutable
import BLAS.Matrix.Solve.Mutable

import BLAS.Tensor
import BLAS.Types( flipUpLo )

import Data.Vector.Dense.Class.Internal( IOVector, STVector,
    BaseVector(..), ReadVector, WriteVector, doConjVector,
    withVectorPtr, stride, isConj, coerceVector )
import Data.Vector.Dense.Class.Copying( newCopyVector, unsafeCopyVector )
import Data.Vector.Dense.Class.Creating( newListVector )
import Data.Vector.Dense.Class.Operations( getConjVector, axpyVector )

import Data.Matrix.Herm
import Data.Matrix.Tri.Internal
import Data.Matrix.Dense.Class( BaseMatrix, ReadMatrix, WriteMatrix,
    isHermMatrix, arrayFromMatrix, matrixViewArray, colViews,
    coerceMatrix, newCopyMatrix, unsafeCopyMatrix, axpyMatrix )


class (MatrixShaped a e, Storable e) => BaseBanded_ a e where
    type VectorViewB a :: * -> * -> *
    bandedViewArray :: ForeignPtr e -> Ptr e -> Int -> Int -> Int -> Int -> Int -> Bool -> a mn e
    arrayFromBanded :: a mn e -> (ForeignPtr e, Ptr e, Int, Int, Int, Int, Int, Bool)

class (BaseBanded_ a e, BaseVector (VectorViewB a) e) => BaseBanded a e

class ( BaseBanded a e, BLAS2 e, UnsafeIOToM m, ReadTensor a (Int,Int) e m
      , MMatrix a e m, MMatrix (Herm a) e m, MMatrix (Tri a) e m
      , MSolve (Tri a) e m
      , ReadVector (VectorViewB a) e m) => 
    ReadBanded a e m where

class (ReadBanded a e m, WriteTensor a (Int,Int) e m,
           WriteVector (VectorViewB a) e m) => 
    WriteBanded a e m | m -> a where


------------------------- Basic Banded Properties ---------------------------

withBandedPtr :: (BaseBanded a e) => 
    a mn e -> (Ptr e -> IO b) -> IO b
withBandedPtr a f =
    let (fp,p,_,_,_,_,_,_) = arrayFromBanded a
    in do
        b <- f p
        touchForeignPtr fp
        return b

size1 :: (BaseBanded a e) => a mn e -> Int
size1 a = let (_,_,m,_,_,_,_,_) = arrayFromBanded a in m
{-# INLINE size1  #-}

size2 :: (BaseBanded a e) => a mn e -> Int
size2 a = let (_,_,_,n,_,_,_,_) = arrayFromBanded a in n
{-# INLINE size2 #-}

lowBW :: (BaseBanded a e) => a mn e -> Int
lowBW a = let (_,_,_,_,kl,_,_,_) = arrayFromBanded a in kl
{-# INLINE lowBW  #-}

upBW :: (BaseBanded a e) => a mn e -> Int
upBW a = let (_,_,_,_,_,ku,_,_) = arrayFromBanded a in ku
{-# INLINE upBW #-}

ldaOfBanded :: (BaseBanded a e) => a mn e -> Int
ldaOfBanded a = let (_,_,_,_,_,_,l,_) = arrayFromBanded a in l
{-# INLINE ldaOfBanded #-}

isHermBanded :: (BaseBanded a e) => a mn e -> Bool
isHermBanded a = let (_,_,_,_,_,_,_,h) = arrayFromBanded a in h
{-# INLINE isHermBanded #-}

matrixFromBanded :: (BaseBanded b e, BaseMatrix a e) => 
    b mn e -> ((Int,Int), (Int,Int), a mn' e, Bool)
matrixFromBanded b =
    let (f,p,m,n,kl,ku,ld,h) = arrayFromBanded b
        a = matrixViewArray f p (kl+1+ku) n ld False
    in ((m,n), (kl,ku), a, h)

bandedViewMatrix :: (BaseMatrix a e, BaseBanded b e) => 
    (Int,Int) -> (Int,Int) -> a mn e -> Bool -> Maybe (b mn' e)
bandedViewMatrix (m,n) (kl,ku) a h = 
    if isHermMatrix a 
        then Nothing
        else let (f,p,m',n',ld,_) = arrayFromMatrix a
             in case undefined of
                 _ | m' /= kl+1+ku -> 
                     error $ "bandedViewMatrix:"
                        ++ " number of rows must be equal to number of diagonals"
                 _ | n' /= n ->
                     error $ "bandedViewMatrix:"
                        ++ " numbers of columns must be equal"
                 _ ->
                     Just $ bandedViewArray f p m n kl ku ld h

bandwidth :: (BaseBanded a e) => a mn e -> (Int,Int)
bandwidth a =
    let (kl,ku) = (numLower a, numUpper a)
    in (negate kl, ku)
{-# INLINE bandwidth #-}

numLower :: (BaseBanded a e) =>  a mn e -> Int
numLower a | isHermBanded a = upBW a
           | otherwise      = lowBW a
{-# INLINE numLower #-}

numUpper :: (BaseBanded a e) =>  a mn e -> Int
numUpper a | isHermBanded a = lowBW a
           | otherwise      = upBW a
{-# INLINE numUpper #-}


-- | Cast the shape type of the matrix.
coerceBanded :: (BaseBanded a e) => a mn e -> a mn' e
coerceBanded = unsafeCoerce
{-# INLINE coerceBanded #-}


-------------------------- BaseTensor functions -----------------------------

shapeBanded :: (BaseBanded a e) => a mn e -> (Int,Int)
shapeBanded a | isHermBanded a = (size2 a, size1 a)
              | otherwise      = (size1 a, size2 a)
{-# INLINE shapeBanded #-}

boundsBanded :: (BaseBanded a e) => a mn e -> ((Int,Int), (Int,Int))
boundsBanded a = ((0,0), (m-1,n-1)) where (m,n) = shapeBanded a
{-# INLINE boundsBanded #-}


-------------------------- BaseMatrix functions -----------------------------

hermBanded :: (BaseBanded a e) => a (m,n) e -> a (n,m) e
hermBanded a = let (f,p,m,n,kl,ku,l,h) = arrayFromBanded a
               in bandedViewArray f p m n kl ku l (not h)
{-# INLINE hermBanded #-}


-------------------------- ReadTensor functions -----------------------------

getSizeBanded :: (ReadBanded a e m) => a mn e -> m Int
getSizeBanded = return . sizeBanded
{-# INLINE getSizeBanded #-}

getIndicesBanded :: (ReadBanded a e m) => a mn e -> m [(Int,Int)]
getIndicesBanded = return . indicesBanded
{-# INLINE getIndicesBanded #-}

getElemsBanded :: (ReadBanded a e m) => a mn e -> m [e]
getElemsBanded a = getAssocsBanded a >>= return . (map snd)

getAssocsBanded :: (ReadBanded a e m) => a mn e -> m [((Int,Int),e)]
getAssocsBanded a = do
    is <- getIndicesBanded a
    unsafeInterleaveM $ mapM (\i -> unsafeReadElem a i >>= \e -> return (i,e)) is
    
getIndicesBanded' :: (ReadBanded a e m) => a mn e -> m [(Int,Int)]
getIndicesBanded' = getIndicesBanded
{-# INLINE getIndicesBanded' #-}

getElemsBanded' :: (ReadBanded a e m) => a mn e -> m [e]
getElemsBanded' a = getAssocsBanded' a >>= return . (map snd)

getAssocsBanded' :: (ReadBanded a e m) => a mn e -> m [((Int,Int),e)]
getAssocsBanded' a = do
    is <- getIndicesBanded a
    mapM (\i -> unsafeReadElem a i >>= \e -> return (i,e)) is

unsafeReadElemBanded :: (ReadBanded a e m) => a mn e -> (Int,Int) -> m e
unsafeReadElemBanded a (i,j)
    | isHermBanded a = 
        unsafeReadElemBanded (hermBanded $ coerceBanded a) (j,i) 
        >>= return . conj
    | hasStorageBanded a (i,j) =
        unsafeIOToM $
            withBandedElemPtr a (i,j) peek
    | otherwise =
        return 0
{-# INLINE unsafeReadElemBanded #-}


------------------------- WriteTensor functions -----------------------------

-- | Create a new banded matrix of given shape and (lower,upper), bandwidths,
-- but do not initialize the elements.
newBanded_ :: (WriteBanded a e m) => (Int,Int) -> (Int,Int) -> m (a mn e)
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
        let m'  = kl + 1 + ku
            l   = m'
            h   = False
        in unsafeIOToM $ do    
            fp <- mallocForeignPtrArray (m' * n)
            let p = unsafeForeignPtrToPtr fp
            return $ bandedViewArray fp p m n kl ku l h
    where
      err s = fail $ "newBanded_ " ++ show (m,n) ++ " " ++ show (kl,ku) ++ ": " ++ s
                  
-- | Create a zero banded matrix of the specified shape and bandwidths.
newZeroBanded :: (WriteBanded a e m) => (Int,Int) -> (Int,Int) -> m (a mn e)
newZeroBanded mn bw = do
    a <- newBanded_ mn bw
    setZeroBanded a
    return a

-- | Create a constant banded matrix of the specified shape and bandwidths.
newConstantBanded :: (WriteBanded a e m) => (Int,Int) -> (Int,Int) -> e -> m (a mn e)
newConstantBanded mn bw e = do
    a <- newBanded_ mn bw
    setConstantBanded e a
    return a

setZeroBanded :: (WriteBanded a e m) => a mn e -> m ()    
setZeroBanded = setConstantBanded 0

setConstantBanded :: (WriteBanded a e m) => e -> a mn e -> m ()
setConstantBanded e a
    | isHermBanded a = setConstantBanded (conj e) a'
    | otherwise = do
        is <- getIndicesBanded a
        mapM_ (\i -> unsafeWriteElemBanded a i e) is
  where
    a' = (hermBanded . coerceBanded) a

unsafeWriteElemBanded :: (WriteBanded a e m) => 
    a mn e -> (Int,Int) -> e -> m ()
unsafeWriteElemBanded a (i,j) e
    | isHermBanded a  = unsafeWriteElemBanded a' (j,i) $ conj e
    | otherwise = unsafeIOToM $
                      withBandedElemPtr a (i,j) (`poke` e)
  where
    a' = (hermBanded . coerceBanded) a

modifyWithBanded :: (WriteBanded a e m) => (e -> e) -> a mn e -> m ()
modifyWithBanded f a = do
    ies <- getAssocsBanded a
    mapM_ (\(ij,e) -> unsafeWriteElemBanded a ij (f e)) ies

canModifyElemBanded :: (WriteBanded a e m) => a mn e -> (Int,Int) -> m Bool
canModifyElemBanded a ij = return $ hasStorageBanded a ij
{-# INLINE canModifyElemBanded #-}


------------------------------ Vector views ---------------------------------

unsafeRowViewBanded :: (BaseBanded a e) => 
    a mn e -> Int -> (Int, VectorViewB a k e, Int)
unsafeRowViewBanded a i =
    if h then
        case unsafeColViewBanded a' i of (nb, v, na) -> (nb, conj v, na)        
    else
        let nb  = max (i - kl)         0
            na  = max (n - 1 - i - ku) 0
            r   = min (ku + i)         (kl + ku)
            c   = max (i - kl)         0 
            p'  = p `advancePtr` (r + c * ld)
            inc = ld - 1
            len = n - (nb + na)
        in if len >= 0 
            then (nb, vectorViewArray f p' len inc False, na)
            else (n , vectorViewArray f p' 0   inc False,  0)
  where
    (f,p,_,n,kl,ku,ld,h) = arrayFromBanded a
    a' = (hermBanded . coerceBanded) a

unsafeColViewBanded :: (BaseBanded a e) => 
    a mn e -> Int -> (Int, VectorViewB a k e, Int)
unsafeColViewBanded a j =
    if h then
        case unsafeRowViewBanded a' j of (nb, v, na) -> (nb, conj v, na)
    else
        let nb  = max (j - ku)         0
            na  = max (m - 1 - j - kl) 0
            r   = max (ku - j) 0 
            c   = j 
            p'  = p `advancePtr` (r + c * ld)
            inc = 1
            len = m - (nb + na)
        in if len >= 0
            then (nb, vectorViewArray f p' len inc False, na)
            else (m , vectorViewArray f p' 0   inc False,  0)
  where
    (f,p,m,_,kl,ku,ld,h) = arrayFromBanded a
    a' = (hermBanded . coerceBanded) a

unsafeGetRowBanded :: (ReadBanded a e m, WriteVector y e m) => 
    a (k,l) e -> Int -> m (y l e)
unsafeGetRowBanded a i = 
    let (nb,x,na) = unsafeRowViewBanded a i
        n = numCols a
    in do
        es <- getElems x
        newListVector n $ (replicate nb 0) ++ es ++ (replicate na 0)

unsafeGetColBanded :: (ReadBanded a e m, WriteVector y e m) => 
    a (k,l) e -> Int -> m (y k e)
unsafeGetColBanded a j = unsafeGetRowBanded (hermBanded a) j >>= return . conj


-------------------------- Matrix multiplication ----------------------------

-- | @gbmv alpha a x beta y@ replaces @y := alpha a * x + beta y@
gbmv :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) => 
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
gbmv alpha a x beta y
    | numRows a == 0 || numCols a == 0 =
        scaleBy beta y
    | isConj x = do
        x' <- getConjVector (conj x)
        gbmv alpha a x' beta y
    | isConj y = do
        doConjVector y
        gbmv alpha a x beta (conj y)
        doConjVector y
    | otherwise =
        let transA = blasTransOf a
            (m,n)  = case (isHermBanded a) of
                         False -> shape a
                         True  -> (flipShape . shape) a
            (kl,ku) = case (isHermBanded a) of
                          False -> (numLower a, numUpper a)
                          True  -> (numUpper a, numLower a)
            ldA    = ldaOfBanded a
            incX   = stride x
            incY   = stride y
        in unsafeIOToM $
               withBandedPtr a $ \pA ->
               withVectorPtr x $ \pX ->
               withVectorPtr y $ \pY -> do
                   BLAS.gbmv transA m n kl ku alpha pA ldA pX incX beta pY incY

-- | @gbmm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gbmm :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
gbmm alpha a b beta c =
    sequence_ $
        zipWith (\x y -> gbmv alpha a x beta y) (colViews b) (colViews c)

hbmv :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) => 
    e -> Herm a (k,k) e -> x k e -> e -> y k e -> m ()
hbmv alpha h x beta y
    | numRows h == 0 =
        return ()
    | isConj y = do
        doConj y
        hbmv alpha h x beta (conj y)
        doConj y
    | isConj x = do
        x' <- newCopyVector x
        doConj x'
        hbmv alpha h (conj x') beta y
    | otherwise =
        let (u,a) = hermToBase h
            n     = numCols a
            k     = case u of 
                        Upper -> numUpper a
                        Lower -> numLower a      
            u'    = case (isHermBanded a) of
                        True  -> flipUpLo u
                        False -> u
            uploA = cblasUpLo u'
            ldA   = ldaOfBanded a
            incX  = stride x
            incY  = stride y
            withPtrA 
                  = case u' of Upper -> withBandedPtr a
                               Lower -> withBandedElemPtr a (0,0)
        in unsafeIOToM $
               withPtrA $ \pA ->
               withVectorPtr x $ \pX ->
               withVectorPtr y $ \pY -> do
                   BLAS.hbmv uploA n k alpha pA ldA pX incX beta pY incY

hbmm :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> Herm a (k,k) e -> b (k,l) e -> e -> c (k,l) e -> m ()
hbmm alpha h b beta c =
    zipWithM_ (\x y -> hbmv alpha h x beta y) (colViews b) (colViews c)

hbmv' :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) => 
    e -> Herm a (r,s) e -> x s e -> e -> y r e -> m ()
hbmv' alpha a x beta y = 
    hbmv alpha (coerceHerm a) x beta (coerceVector y)

hbmm' :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> Herm a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
hbmm' alpha a b beta c = 
    hbmm alpha (coerceHerm a) b beta (coerceMatrix c)

tbmv :: (ReadBanded a e m, WriteVector y e m) => 
    e -> Tri a (k,k) e -> y n e -> m ()
tbmv alpha t x | isConj x = do
    doConj x
    tbmv alpha t (conj x)
    doConj x

tbmv alpha t x =
    let (u,d,a) = triToBase t
        (transA,u') 
                  = if isHermBanded a 
                        then (conjTrans, flipUpLo u) else (noTrans, u)
        uploA     = cblasUpLo u'
        diagA     = cblasDiag d
        n         = numCols a
        k         = case u of Upper -> numUpper a 
                              Lower -> numLower a
        ldA       = ldaOfBanded a
        incX      = stride x
        withPtrA  = case u' of 
                        Upper -> withBandedPtr a
                        Lower -> withBandedElemPtr a (0,0)
    in do
        scaleBy alpha x
        unsafeIOToM $
            withPtrA $ \pA ->
            withVectorPtr x $ \pX -> do
                BLAS.tbmv uploA transA diagA n k pA ldA pX incX

tbmm :: (ReadBanded a e m, WriteMatrix b e m) =>
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
tbmm 1     t b = mapM_ (\x -> tbmv 1 t x) (colViews b)
tbmm alpha t b = scaleBy alpha b >> tbmm 1 t b

tbmv' :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) => 
    e -> Tri a (r,s) e -> x s e -> e -> y r e -> m ()
tbmv' alpha a x beta y 
    | beta /= 0 = do
        x' <- newCopyVector x
        tbmv alpha (coerceTri a) x'
        scaleBy beta y
        axpyVector 1 x' (coerceVector y)
    | otherwise = do
        unsafeCopyVector (coerceVector y) x
        tbmv alpha (coerceTri a) (coerceVector y)

tbmm' :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) => 
    e -> Tri a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
tbmm' alpha a b beta c
    | beta /= 0 = do
        b' <- newCopyMatrix b
        tbmm alpha (coerceTri a) b'
        scaleBy beta c
        axpyMatrix 1 b' (coerceMatrix c)
    | otherwise = do
        unsafeCopyMatrix (coerceMatrix c) b
        tbmm alpha (coerceTri a) (coerceMatrix c)

tbsv :: (ReadBanded a e m, WriteVector y e m) => 
    e -> Tri a (k,k) e -> y n e -> m ()
tbsv alpha t x | isConj x = do
    doConj x
    tbsv alpha t (conj x)
    doConj x
    
tbsv alpha t x = 
    let (u,d,a) = triToBase t
        (transA,u') = if isHermBanded a then (conjTrans, flipUpLo u) else (noTrans, u)
        uploA     = cblasUpLo u'
        diagA     = cblasDiag d
        n         = numCols a
        k         = case u of Upper -> numUpper a 
                              Lower -> numLower a        
        ldA       = ldaOfBanded a
        incX      = stride x
        withPtrA  = case u' of 
                        Upper -> withBandedPtr a
                        Lower -> withBandedElemPtr a (0,0)
    in do
        scaleBy alpha x
        unsafeIOToM $
            withPtrA $ \pA ->
            withVectorPtr x $ \pX -> do
                BLAS.tbsv uploA transA diagA n k pA ldA pX incX

tbsm :: (ReadBanded a e m, WriteMatrix b e m) => 
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
tbsm 1     t b = mapM_ (\x -> tbsv 1 t x) (colViews b)
tbsm alpha t b = scaleBy alpha b >> tbsm 1 t b

unsafeDoSSolveTriBanded :: (ReadBanded a e m,
    ReadVector y e m, WriteVector x e m) =>
        e -> Tri a (k,l) e -> y k e -> x l e -> m ()
unsafeDoSSolveTriBanded alpha a y x = do
    unsafeCopyVector (coerceVector x) y
    tbsv alpha (coerceTri a) (coerceVector x)

unsafeDoSSolveMatTriBanded :: (ReadBanded a e m,
    ReadMatrix c e m, WriteMatrix b e m) =>
        e -> Tri a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
unsafeDoSSolveMatTriBanded alpha a c b = do
    unsafeCopyMatrix (coerceMatrix b) c
    tbsm alpha (coerceTri a) b


--------------------------- Utility functions -------------------------------

withBandedElemPtr :: (BaseBanded a e) => 
    a mn e -> (Int,Int) -> (Ptr e -> IO b) -> IO b
withBandedElemPtr a (i,j) f
    | isHermBanded a  = withBandedElemPtr (hermBanded $ coerceBanded a) (j,i) f
    | otherwise = withBandedPtr a $ \ptr ->
                      f $ ptr `advancePtr` (indexOfBanded a (i,j))

indexOfBanded :: (BaseBanded a e) => a mn e -> (Int,Int) -> Int
indexOfBanded a (i,j) =
    let (_,_,_,_,_,ku,ld,h) = arrayFromBanded a
        (i',j')           = if h then (j,i) else (i,j)
    in ku + (i' - j') + j' * ld

hasStorageBanded :: (BaseBanded a e) => a mn e -> (Int,Int) -> Bool
hasStorageBanded a (i,j) =
    let (_,_,m,_,kl,ku,_,h) = arrayFromBanded a
        (i',j')             = if h then (j,i) else (i,j)
    in inRange (max 0 (j'-ku), min (m-1) (j'+kl)) i'

sizeBanded :: (BaseBanded a e) => a mn e -> Int
sizeBanded a =
    let (_,_,m,n,kl,ku,_,_) = arrayFromBanded a
    in foldl' (+) 0 $ map (diagLen (m,n)) [(-kl)..ku]

indicesBanded :: (BaseBanded a e) => a mn e -> [(Int,Int)]
indicesBanded a =
    let is = if isHermBanded a 
                 then [ (i,j) | i <- range (0,m-1), j <- range (0,n-1) ]
                 else [ (i,j) | j <- range (0,n-1), i <- range (0,m-1) ]
    in filter (hasStorageBanded a) is
  where (m,n) = shapeBanded a

blasTransOf :: (BaseBanded a e) => a mn e -> CBLASTrans
blasTransOf a = 
    case (isHermBanded a) of
          False -> noTrans
          True  -> conjTrans

flipShape :: (Int,Int) -> (Int,Int)
flipShape (m,n) = (n,m)


------------------------------------ Instances ------------------------------

-- | The Banded matrix data type.
data IOBanded mn e =
    BM {-# UNPACK #-} !(ForeignPtr e) -- storage
       {-# UNPACK #-} !(Ptr e)        -- base pointer
       {-# UNPACK #-} !Int            -- numer of rows
       {-# UNPACK #-} !Int            -- number of columns
       {-# UNPACK #-} !Int            -- lower bandwidth
       {-# UNPACK #-} !Int            -- upper bandwidth
       {-# UNPACK #-} !Int            -- lda of storage
       {-# UNPACK #-} !Bool           -- isHerm flag

newtype STBanded s mn e = ST (IOBanded mn e)

unsafeIOBandedToSTBanded :: IOBanded mn e -> STBanded s mn e
unsafeIOBandedToSTBanded = ST

unsafeSTBandedToIOBanded :: STBanded s mn e -> IOBanded mn e
unsafeSTBandedToIOBanded (ST x) = x

instance (Storable e) => BaseBanded_ IOBanded e where
    type VectorViewB IOBanded = IOVector
    bandedViewArray f p m n kl ku ld h      = BM f p m n kl ku ld h
    arrayFromBanded (BM f p m n kl ku ld h) = (f,p,m,n,kl,ku,ld,h)

instance (Storable e) => BaseBanded_ (STBanded s) e where
    type VectorViewB (STBanded s) = (STVector s)
    bandedViewArray f p m n kl ku ld h           = ST (BM f p m n kl ku ld h)
    arrayFromBanded (ST (BM f p m n kl ku ld h)) = (f,p,m,n,kl,ku,ld,h)

instance (Storable e) => BaseBanded IOBanded e
instance (Storable e) => BaseBanded (STBanded s) e

instance (Storable e) => BaseTensor IOBanded (Int,Int) e where
    shape  = shapeBanded
    bounds = boundsBanded
    
instance (Storable e) => BaseTensor (STBanded s) (Int,Int) e where
    shape  = shapeBanded
    bounds = boundsBanded

instance (Storable e) => MatrixShaped IOBanded e where
    herm = hermBanded
    
instance (Storable e) => MatrixShaped (STBanded s) e where
    herm = hermBanded

instance (BLAS2 e) => ReadBanded IOBanded     e IO
instance (BLAS2 e) => ReadBanded (STBanded s) e (ST s)

instance (BLAS2 e) => ReadTensor IOBanded (Int,Int) e IO where
    getSize        = getSizeBanded
    getAssocs      = getAssocsBanded
    getIndices     = getIndicesBanded
    getElems       = getElemsBanded
    getAssocs'     = getAssocsBanded'
    getIndices'    = getIndicesBanded'
    getElems'      = getElemsBanded'
    unsafeReadElem = unsafeReadElemBanded
    
instance (BLAS2 e) => ReadTensor (STBanded s) (Int,Int) e (ST s) where
    getSize        = getSizeBanded
    getAssocs      = getAssocsBanded
    getIndices     = getIndicesBanded
    getElems       = getElemsBanded
    getAssocs'     = getAssocsBanded'
    getIndices'    = getIndicesBanded'
    getElems'      = getElemsBanded'
    unsafeReadElem = unsafeReadElemBanded

instance (BLAS2 e) => WriteBanded IOBanded     e IO where
instance (BLAS2 e) => WriteBanded (STBanded s) e (ST s) where

instance (BLAS2 e) => WriteTensor IOBanded (Int,Int) e IO where
    setConstant     = setConstantBanded
    setZero         = setZeroBanded
    modifyWith      = modifyWithBanded
    unsafeWriteElem = unsafeWriteElemBanded
    canModifyElem   = canModifyElemBanded

instance (BLAS2 e) => WriteTensor (STBanded s) (Int,Int) e (ST s) where
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

instance (BLAS2 e) => MMatrix (Herm (STBanded s)) e (ST s) where
    unsafeDoSApplyAdd    = hbmv'
    unsafeDoSApplyAddMat = hbmm'

instance (BLAS2 e) => MMatrix (Herm IOBanded) e IO where
    unsafeDoSApplyAdd    = hbmv'
    unsafeDoSApplyAddMat = hbmm'

instance (BLAS2 e) => MMatrix (Tri (STBanded s)) e (ST s) where
    unsafeDoSApply_      = tbmv
    unsafeDoSApplyMat_   = tbmm
    unsafeDoSApplyAdd    = tbmv'
    unsafeDoSApplyAddMat = tbmm'

instance (BLAS2 e) => MMatrix (Tri IOBanded) e IO where
    unsafeDoSApply_      = tbmv
    unsafeDoSApplyMat_   = tbmm
    unsafeDoSApplyAdd    = tbmv'
    unsafeDoSApplyAddMat = tbmm'

instance (BLAS2 e) => MSolve (Tri IOBanded) e IO where
    unsafeDoSSolve     = unsafeDoSSolveTriBanded
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriBanded
    unsafeDoSSolve_    = tbsv
    unsafeDoSSolveMat_ = tbsm

instance (BLAS2 e) => MSolve (Tri (STBanded s)) e (ST s) where
    unsafeDoSSolve     = unsafeDoSSolveTriBanded
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriBanded
    unsafeDoSSolve_    = tbsv
    unsafeDoSSolveMat_ = tbsm
