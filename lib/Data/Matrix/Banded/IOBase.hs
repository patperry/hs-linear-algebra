{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances,
        TypeFamilies, ScopedTypeVariables, GADTs #-}
{-# OPTIONS_HADDOCK hide #-}
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
import Text.Printf
import Unsafe.Coerce

import BLAS.Internal( clearArray, diagLen, diagStart )
import Data.Elem.BLAS.Base( Elem, conjugate )
import qualified Data.Elem.BLAS.Base   as BLAS
import Data.Elem.BLAS.Level2( BLAS2 )
import qualified Data.Elem.BLAS.Level2 as BLAS

import Data.Matrix.Class
import Data.Matrix.Class.MMatrixBase
import Data.Matrix.Class.MSolveBase

import Data.Matrix.Herm
import Data.Matrix.Tri

import Data.Tensor.Class
import Data.Tensor.Class.MTensor

import Data.Matrix.Dense.IOBase( IOMatrix(..) )
import Data.Matrix.Dense.Base( ReadMatrix, WriteMatrix, coerceMatrix,
    newCopyMatrix, colViews, scaleByMatrix, unsafeCopyMatrix, 
    unsafeAxpyMatrix, unsafeMatrixToIOMatrix )
import Data.Vector.Dense.IOBase( IOVector(..), withIOVector )
import Data.Vector.Dense.Base( ReadVector, WriteVector, dim, coerceVector,
    stride, conj, isConj, unsafeVectorToIOVector, newListVector,
    newCopyVector, newCopyVector', scaleByVector, doConjVector,
    unsafeCopyVector, unsafeAxpyVector, unsafeConvertIOVector )

-- | Banded matrix in the 'IO' monad.  The type arguments are as follows:
--
--     * @np@: a phantom type for the shape of the matrix.  Most functions
--       will demand that this be specified as a pair.  When writing a function
--       signature, you should always prefer @IOBanded (n,p) e@ to
--       @IOBanded np e@.
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
data IOBanded np e =
       Elem e =>
       IOBanded { transEnumIOBanded :: {-# UNPACK #-} !TransEnum
                , numRowsIOBanded   :: {-# UNPACK #-} !Int
                , numColsIOBanded   :: {-# UNPACK #-} !Int
                , numLowerIOBanded  :: {-# UNPACK #-} !Int
                , numUpperIOBanded  :: {-# UNPACK #-} !Int
                , fptrIOBanded      :: {-# UNPACK #-} !(ForeignPtr e)
                , ptrIOBanded       :: {-# UNPACK #-} !(Ptr e)
                , ldaIOBanded       :: {-# UNPACK #-} !Int
                }

hermIOBanded :: IOBanded np e -> IOBanded pn e
hermIOBanded (IOBanded h m n kl ku f p l) = 
    (IOBanded (flipTrans h) n m ku kl f p l)
{-# INLINE hermIOBanded #-}

isHermIOBanded :: IOBanded np e -> Bool
isHermIOBanded = (ConjTrans ==) . transEnumIOBanded
{-# INLINE isHermIOBanded #-}

indexIOBanded :: IOBanded np e -> (Int,Int) -> Int
indexIOBanded (IOBanded h _ _ kl ku _ _ ld) (i,j) =
    let (ku',i',j') = if h == ConjTrans then (kl,j,i) else (ku,i,j)
    in ku' + (i' - j') + j' * ld
{-# INLINE indexIOBanded #-}

hasStorageIOBanded :: IOBanded np e -> (Int,Int) -> Bool
hasStorageIOBanded (IOBanded h m n kl ku _ _ _) (i,j) =
    let (m',kl',ku',i',j') = if h == ConjTrans
                                 then (n,ku,kl,j,i) 
                                 else (m,kl,ku,i,j)
    in inRange (max 0 (j'-ku'), min (m'-1) (j'+kl')) i'
{-# INLINE hasStorageIOBanded #-}

-- | Execute an 'IO' action with a pointer to the first element in the
-- banded matrix.
withIOBanded :: IOBanded (n,p) e -> (Ptr e -> IO a) -> IO a
withIOBanded (IOBanded _ _ _ _ _ f p _) g = do
    a <- g p
    touchForeignPtr f
    return a
{-# INLINE withIOBanded #-}

withIOBandedElem :: IOBanded np e -> (Int,Int) -> (Ptr e -> IO a) -> IO a
withIOBandedElem a@(IOBanded _ _ _ _ _ _ _ _) (i,j) f =
    withIOBanded (coerceIOBanded a) $ \p ->
        f $ p `advancePtr` (indexIOBanded a (i,j))
{-# INLINE withIOBandedElem #-}

maybeMatrixStorageFromIOBanded :: IOBanded (n,p) e -> Maybe (IOMatrix (k,p) e)
maybeMatrixStorageFromIOBanded (IOBanded h _ n kl ku f p l)
    | h == ConjTrans = Nothing
    | otherwise      = Just $ IOMatrix NoTrans (kl+1+ku) n f p l
{-# INLINE maybeMatrixStorageFromIOBanded #-}

maybeIOBandedFromMatrixStorage :: (Int,Int)
                               -> (Int,Int)
                               -> IOMatrix np e
                               -> Maybe (IOBanded np' e)
maybeIOBandedFromMatrixStorage (m,n) (kl,ku) (IOMatrix h m' n' f p l)
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
    | m' /= kl+1+ku =
        err $ "number of rows in the underlying matrix"
              ++ " must be equal to number of diagonals."
    | n' /= n =
        err $ "numbers of columns must be equal"
            ++ " to the number of columns in the underlying matrix."
    | h == ConjTrans =
        Nothing    
    | otherwise =
        Just $ IOBanded NoTrans m n kl ku f p l
  where
    err s = 
      error $ "maybeBandedFromMatrixStorage " ++ show (m,n) ++ " " ++ show (kl,ku) 
            ++ " <matrix of shape " ++ show (m',n') ++ ">: " ++ show s

viewVectorAsIOBanded :: (Int,Int) -> IOVector k e -> IOBanded (n,p) e
viewVectorAsIOBanded (m,n) (IOVector c k f p l)
    | k /= m && k /= n =
        error $ "viewVectorAsBanded " ++ show (m,n) ++ " "
              ++ " <vector of dim " ++ show k ++ ">: "
              ++ "vector must have length equal to one of the specified"
              ++ " diemsion sizes"
    | otherwise =
        let h = if c == Conj then ConjTrans else NoTrans
        in IOBanded h m n 0 0 f p l
{-# INLINE viewVectorAsIOBanded #-}

maybeViewIOBandedAsVector :: IOBanded (n,p) e -> Maybe (IOVector k e)
maybeViewIOBandedAsVector (IOBanded h m n kl ku f p l)
    | kl == 0 && ku == 0 =
        let c = if h == ConjTrans then Conj else NoConj
        in Just $ IOVector c (min m n) f p l
    | otherwise =
        Nothing
{-# INLINE maybeViewIOBandedAsVector #-}

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
sizeIOBanded (IOBanded _ m n kl ku _ _ _) =
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

getAssocsIOBanded :: IOBanded np e -> IO [((Int,Int),e)]
getAssocsIOBanded a =
    let go []     = return []
        go (i:is) = unsafeInterleaveIO $ do
                        e   <- unsafeReadElemIOBanded a i
                        ies <- go is
                        return $ (i,e):ies
    in go (indicesIOBanded a)
{-# INLINE getAssocsIOBanded #-}

getAssocsIOBanded' :: IOBanded np e -> IO [((Int,Int),e)]
getAssocsIOBanded' a =
    sequence [ liftM ((,) i) (unsafeReadElemIOBanded a i)
             | i <- indicesIOBanded a ]
{-# INLINE getAssocsIOBanded' #-}

getElemsIOBanded :: IOBanded np e -> IO [e]
getElemsIOBanded a = liftM (snd . unzip) $ getAssocsIOBanded a
{-# INLINE getElemsIOBanded #-}

getElemsIOBanded' :: IOBanded np e -> IO [e]
getElemsIOBanded' a = liftM (snd . unzip) $ getAssocsIOBanded a
{-# INLINE getElemsIOBanded' #-}

unsafeReadElemIOBanded :: IOBanded np e -> (Int,Int) -> IO e
unsafeReadElemIOBanded a@(IOBanded _ _ _ _ _ _ _ _) (i,j)
    | hasStorageIOBanded a (i,j) =
        let f = if isHermIOBanded a then conjugate else id
        in liftM f $ withIOBandedElem a (i,j) peek
    | otherwise =
        return 0
{-# INLINE unsafeReadElemIOBanded #-}

unsafeWriteElemIOBanded :: IOBanded np e -> (Int,Int) -> e -> IO ()
unsafeWriteElemIOBanded a@(IOBanded _ _ _ _ _ _ _ _) ij e =
    let e' = if isHermIOBanded a then conjugate e else e
    in withIOBandedElem a ij (`poke` e')
{-# INLINE unsafeWriteElemIOBanded #-}

modifyWithIOBanded :: (e -> e) -> IOBanded np e -> IO ()
modifyWithIOBanded f a = do
    ies <- getAssocsIOBanded a
    forM_ ies $ \(i,e) -> do
         unsafeWriteElemIOBanded a i (f e)

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
            h   = NoTrans
        in do
            fp <- mallocForeignPtrArray (m' * n)
            let p = unsafeForeignPtrToPtr fp
            return $ IOBanded h m n kl ku fp p l
    where
      err s = fail $ "newBanded_ " ++ show (m,n) ++ " " ++ show (kl,ku) ++ ": " ++ s

newIOBanded :: (Elem e) 
            => (Int,Int)
            -> (Int,Int)
            -> [((Int,Int), e)]
            -> IO (IOBanded (n,p) e)
newIOBanded (m,n) (kl,ku) ies = do
    a <- newIOBanded_ (m,n) (kl,ku)
    withIOBanded a $ \p ->
        clearArray p ((kl+1+ku)*n)
    forM_ ies $ \((i,j),e) -> do
        when (   i < 0 || i >= m
              || j < 0 || j >= n 
              || (not $ hasStorageIOBanded a (i,j))) $ error $
            printf "newBanded (%d,%d) (%d,%d) [ ..., ((%d,%d),_), ... ]: invalid index"
                   m n kl ku i j
        unsafeWriteElemIOBanded a (i,j) e
    return a

unsafeNewIOBanded :: (Elem e) 
                  => (Int,Int) 
                  -> (Int,Int) 
                  -> [((Int,Int), e)] 
                  -> IO (IOBanded (n,p) e)
unsafeNewIOBanded (m,n) (kl,ku) ies = do
    a <- newIOBanded_ (m,n) (kl,ku)
    withIOBanded a $ \p ->
        clearArray p ((kl+1+ku)*n)
    forM_ ies $ \(i,e) ->
        unsafeWriteElemIOBanded a i e
    return a

newListsIOBanded :: (Elem e) 
                 => (Int,Int)
                 -> (Int,Int)
                 -> [[e]]
                 -> IO (IOBanded (n,p) e)
newListsIOBanded (m,n) (kl,ku) xs = do
    a <- newIOBanded_ (m,n) (kl,ku)
    zipWithM_ (writeDiagElems a) [(negate kl)..ku] xs
    return a
  where
    writeDiagElems :: (Elem e) 
                   => IOBanded (n,p) e -> Int -> [e] -> IO ()
    writeDiagElems a i es =
        let d   = unsafeDiagViewIOBanded a i
            nb  = max 0 (-i)
            es' = drop nb es
        in zipWithM_ (unsafeWriteElem d) [0..(dim d - 1)] es'

newZeroIOBanded :: (Elem e) => (Int,Int) -> (Int,Int) -> IO (IOBanded np e)
newZeroIOBanded mn bw = do
    a <- newIOBanded_ mn bw
    setZeroIOBanded a
    return a

setZeroIOBanded :: IOBanded np e -> IO ()
setZeroIOBanded a@(IOBanded _ _ _ _ _ _ _ _) = setConstantIOBanded 0 a
{-# INLINE setZeroIOBanded #-}

newConstantIOBanded :: (Elem e)
                    => (Int,Int) -> (Int,Int) -> e -> IO (IOBanded np e)
newConstantIOBanded mn bw e = do
    a <- newIOBanded_ mn bw
    setConstantIOBanded e a
    return a

setConstantIOBanded :: e -> IOBanded np e -> IO ()
setConstantIOBanded e a@(IOBanded _ _ _ _ _ _ _ _)
    | isHermIOBanded a = setConstantIOBanded (conjugate e) (hermIOBanded a)
    | otherwise = do
        is <- getIndicesIOBanded a
        mapM_ (\i -> unsafeWriteElemIOBanded a i e) is
{-# INLINE setConstantIOBanded #-}

newCopyIOBanded :: IOBanded np e -> IO (IOBanded np e)
newCopyIOBanded a@(IOBanded _ _ _ _ _ _ _ _)
    | isHermIOBanded a =
        liftM hermIOBanded $ newCopyIOBanded (hermIOBanded a)
    | otherwise = do
        a' <- newIOBanded_ (shapeIOBanded a) (numLowerIOBanded a, numUpperIOBanded a)
        unsafeCopyIOBanded a' a
        return a'

unsafeCopyIOBanded :: IOBanded mn e -> IOBanded mn e -> IO ()
unsafeCopyIOBanded dst@(IOBanded _ _ _ _ _ _ _ _) src
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
        BLAS.copy NoConj NoConj (m*n) pSrc 1 pDst 1

    copyCols pDst pSrc nleft
        | nleft == 0 = return ()
        | otherwise = do
            BLAS.copy NoConj NoConj m pSrc 1 pDst 1
            copyCols (pDst `advancePtr` ldDst) (pSrc `advancePtr` ldSrc) 
                     (nleft-1)

    diagViews a = 
        case bandwidthsIOBanded a of
            (kl,ku) -> map (unsafeDiagViewIOBanded a) $ range (-kl,ku)


unsafeRowViewIOBanded :: IOBanded np e
                      -> Int
                      -> (Int, IOVector k e, Int)
unsafeRowViewIOBanded a@(IOBanded h _ n kl ku f p ld) i =
    if h == ConjTrans then
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
        in (nb, IOVector NoConj len f p' inc, na)

unsafeColViewIOBanded :: IOBanded np e
                      -> Int
                      -> (Int, IOVector k e, Int)
unsafeColViewIOBanded a@(IOBanded h m _ kl ku f p ld) j =
    if h == ConjTrans then
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
        in (nb, IOVector NoConj len f p' inc, na)

unsafeDiagViewIOBanded :: IOBanded np e -> Int -> IOVector k e
unsafeDiagViewIOBanded a@(IOBanded h m n _ _ fp p ld) d
    | h == ConjTrans =
         conj $ unsafeDiagViewIOBanded (hermIOBanded a) (negate d)
    | otherwise =
        let off = indexIOBanded a (diagStart d)
            p'  = p `advancePtr` off
            len = diagLen (m,n) d
            inc = ld
        in (IOVector NoConj len fp p' inc)


unsafeGetRowIOBanded :: IOBanded np e -> Int -> IO (IOVector p e)
unsafeGetRowIOBanded a@(IOBanded _ _ _ _ _ _ _ _) i =
    let (nb,x,na) = unsafeRowViewIOBanded a i
        n = numCols a
    in do
        es <- getElems x
        newListVector n $ (replicate nb 0) ++ es ++ (replicate na 0)

unsafeGetColIOBanded :: IOBanded np e -> Int -> IO (IOVector n e)
unsafeGetColIOBanded a j =
    liftM conj $ unsafeGetRowIOBanded (hermIOBanded a) j

-- | @gbmv alpha a x beta y@ replaces @y := alpha a * x + beta y@
gbmv :: (ReadVector x IO, WriteVector y IO, BLAS2 e)
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
        let transA = transEnumIOBanded a
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
gbmm :: (ReadMatrix b IO, WriteMatrix c IO, BLAS2 e)
     => e -> IOBanded (n,p) e -> b (p,q) e -> e -> c (n,q) e -> IO ()
gbmm alpha a b beta c =
    sequence_ $
        zipWith (\x y -> gbmv alpha a x beta y) (colViews b) (colViews c)

unsafeGetColHermIOBanded :: (WriteVector y IO, Elem e)
                         => Herm IOBanded (n,p) e -> Int -> IO (y n e)
unsafeGetColHermIOBanded =
    error "TODO: unsafeGetColHermIOBanded is not implemented"

unsafeGetRowHermIOBanded :: (WriteVector y IO, Elem e)
                         => Herm IOBanded (n,p) e -> Int -> IO (y p e)
unsafeGetRowHermIOBanded =
    error "TODO: unsafeGetRowHermIOBanded is not implemented"

hbmv :: (ReadVector x IO, BLAS2 e) => 
    e -> Herm IOBanded (n,p) e -> x p e -> e -> IOVector n e -> IO ()
hbmv alpha h (x :: x p e) beta y
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

hbmm :: (ReadMatrix b IO, BLAS2 e) => 
    e -> Herm IOBanded (n,p) e -> b (p,q) e -> e -> IOMatrix (n,q) e -> IO ()
hbmm alpha h b beta c =
    zipWithM_ (\x y -> hbmv alpha h x beta y) (colViews b) (colViews c)

unsafeGetColTriIOBanded :: (WriteVector y IO, Elem e)
                        => Tri IOBanded (n,p) e -> Int -> IO (y n e)
unsafeGetColTriIOBanded =
    error "TODO: unsafeGetColTriIOBanded is not implemented"

unsafeGetRowTriIOBanded :: (WriteVector y IO, Elem e)
                        => Tri IOBanded (n,p) e -> Int -> IO (y p e)
unsafeGetRowTriIOBanded =
    error "TODO: unsafeGetRowTriIOBanded is not implemented"

tbmv :: (WriteVector y IO, BLAS2 e) => 
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

tbmm :: (WriteMatrix b IO, BLAS2 e) =>
    e -> Tri IOBanded (n,n) e -> b (n,p) e -> IO ()
tbmm 1     t b = mapM_ (\x -> tbmv 1 t x) (colViews b)
tbmm alpha t b = scaleByMatrix alpha b >> tbmm 1 t b

tbmv' :: (ReadVector x IO, WriteVector y IO, BLAS2 e) => 
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

tbmm' :: (ReadMatrix b IO, WriteMatrix c IO, BLAS2 e) => 
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

tbsv :: (WriteVector y IO, BLAS2 e) => 
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

tbsm :: (WriteMatrix b IO, BLAS2 e) => 
    e -> Tri IOBanded (k,k) e -> b (k,l) e -> IO ()
tbsm 1     t b = mapM_ (\x -> tbsv 1 t x) (colViews b)
tbsm alpha t b = scaleByMatrix alpha b >> tbsm 1 t b

tbsv' :: (ReadVector y IO, WriteVector x IO, BLAS2 e)
      => e -> Tri IOBanded (k,l) e -> y k e -> x l e -> IO ()
tbsv' alpha a y x = do
    unsafeCopyVector (coerceVector x) y
    tbsv alpha (coerceTri a) (coerceVector x)

tbsm' :: (ReadMatrix c IO, WriteMatrix b IO, BLAS2 e) 
      => e -> Tri IOBanded (r,s) e -> c (r,t) e -> b (s,t) e -> IO ()
tbsm' alpha a c b = do
    unsafeCopyMatrix (coerceMatrix b) c
    tbsm alpha (coerceTri a) b

instance HasVectorView IOBanded where
    type VectorView IOBanded = IOVector

instance HasMatrixStorage IOBanded where
    type MatrixStorage IOBanded = IOMatrix

instance Shaped IOBanded (Int,Int) where
    shape  = shapeIOBanded
    {-# INLINE shape #-}
    bounds = boundsIOBanded
    {-# INLINE bounds #-}

instance MatrixShaped IOBanded where
    herm = hermIOBanded
    {-# INLINE herm #-}

instance ReadTensor IOBanded (Int,Int) IO where
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

instance WriteTensor IOBanded (Int,Int) IO where
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

instance MMatrix IOBanded IO where
    unsafeDoSApplyAddVector    = gbmv
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = gbmm
    {-# INLINE unsafeDoSApplyAddMatrix #-}
    unsafeGetRow a i = unsafeConvertIOVector $ unsafeGetRowIOBanded a i
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol a j = unsafeConvertIOVector $ unsafeGetColIOBanded a j
    {-# INLINE unsafeGetCol #-}
    getRows = getRowsIO
    {-# INLINE getRows #-}
    getCols = getColsIO
    {-# INLINE getCols #-}

instance MMatrix (Herm IOBanded) IO where
    unsafeDoSApplyAddVector alpha a x beta y = 
        hbmv alpha a x beta (unsafeVectorToIOVector y)
    {-# INLINE unsafeDoSApplyAddVector #-}    
    unsafeDoSApplyAddMatrix alpha a b beta c = 
        hbmm alpha a b beta (unsafeMatrixToIOMatrix c)
    {-# INLINE unsafeDoSApplyAddMatrix #-}    
    getRows = getRowsIO
    {-# INLINE getRows #-}
    getCols = getColsIO
    {-# INLINE getCols #-}
    unsafeGetRow = unsafeGetRowHermIOBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColHermIOBanded
    {-# INLINE unsafeGetCol #-}

instance MMatrix (Tri IOBanded) IO where
    unsafeDoSApplyVector_       = tbmv
    {-# INLINE unsafeDoSApplyVector_  #-}        
    unsafeDoSApplyMatrix_   = tbmm
    {-# INLINE unsafeDoSApplyMatrix_ #-}    
    unsafeDoSApplyAddVector    = tbmv'
    {-# INLINE unsafeDoSApplyAddVector #-}    
    unsafeDoSApplyAddMatrix = tbmm'
    {-# INLINE unsafeDoSApplyAddMatrix #-}    
    getRows = getRowsIO
    {-# INLINE getRows #-}
    getCols = getColsIO
    {-# INLINE getCols #-}
    unsafeGetRow = unsafeGetRowTriIOBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColTriIOBanded
    {-# INLINE unsafeGetCol #-}

instance MSolve (Tri IOBanded) IO where
    unsafeDoSSolveVector_    = tbsv
    {-# INLINE unsafeDoSSolveVector_ #-}    
    unsafeDoSSolveMatrix_ = tbsm
    {-# INLINE unsafeDoSSolveMatrix_ #-}    
    unsafeDoSSolveVector     = tbsv'
    {-# INLINE unsafeDoSSolveVector #-}
    unsafeDoSSolveMatrix  = tbsm'
    {-# INLINE unsafeDoSSolveMatrix #-}    
