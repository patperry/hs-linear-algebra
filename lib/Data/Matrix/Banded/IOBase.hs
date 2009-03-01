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
import Data.Matrix.Dense.Base( ReadMatrix, WriteMatrix,
    newCopyMatrix, colViews, scaleByMatrix, unsafeCopyMatrix, 
    unsafeAxpyMatrix, unsafeMatrixToIOMatrix )
import Data.Vector.Dense.IOBase( IOVector(..), withIOVector )
import Data.Vector.Dense.Base( ReadVector, WriteVector, dim,
    stride, conj, conjEnum, unsafeVectorToIOVector, newListVector,
    newCopyVector, scaleByVector,unsafeCopyVector, unsafeAxpyVector,
    unsafeConvertIOVector, unsafePerformIOWithVector )

-- | Banded matrix in the 'IO' monad.  The type arguments are as follows:
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
data IOBanded e =
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

hermIOBanded :: IOBanded e -> IOBanded e
hermIOBanded (IOBanded h m n kl ku f p l) = 
    (IOBanded (flipTrans h) n m ku kl f p l)
{-# INLINE hermIOBanded #-}

isHermIOBanded :: IOBanded e -> Bool
isHermIOBanded = (ConjTrans ==) . transEnumIOBanded
{-# INLINE isHermIOBanded #-}

indexIOBanded :: IOBanded e -> (Int,Int) -> Int
indexIOBanded (IOBanded h _ _ kl ku _ _ ld) (i,j) =
    let (ku',i',j') = if h == ConjTrans then (kl,j,i) else (ku,i,j)
    in ku' + (i' - j') + j' * ld
{-# INLINE indexIOBanded #-}

hasStorageIOBanded :: IOBanded e -> (Int,Int) -> Bool
hasStorageIOBanded (IOBanded h m n kl ku _ _ _) (i,j) =
    let (m',kl',ku',i',j') = if h == ConjTrans
                                 then (n,ku,kl,j,i) 
                                 else (m,kl,ku,i,j)
    in inRange (max 0 (j'-ku'), min (m'-1) (j'+kl')) i'
{-# INLINE hasStorageIOBanded #-}

-- | Execute an 'IO' action with a pointer to the first element in the
-- banded matrix.
withIOBanded :: IOBanded e -> (Ptr e -> IO a) -> IO a
withIOBanded (IOBanded _ _ _ _ _ f p _) g = do
    a <- g p
    touchForeignPtr f
    return a
{-# INLINE withIOBanded #-}

withIOBandedElem :: IOBanded e -> (Int,Int) -> (Ptr e -> IO a) -> IO a
withIOBandedElem a@(IOBanded _ _ _ _ _ _ _ _) (i,j) f =
    withIOBanded a $ \p ->
        f $ p `advancePtr` (indexIOBanded a (i,j))
{-# INLINE withIOBandedElem #-}

maybeMatrixStorageFromIOBanded :: IOBanded e -> Maybe (IOMatrix e)
maybeMatrixStorageFromIOBanded (IOBanded h _ n kl ku f p l)
    | h == ConjTrans = Nothing
    | otherwise      = Just $ IOMatrix NoTrans (kl+1+ku) n f p l
{-# INLINE maybeMatrixStorageFromIOBanded #-}

maybeIOBandedFromMatrixStorage :: (Int,Int)
                               -> (Int,Int)
                               -> IOMatrix e
                               -> Maybe (IOBanded e)
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

viewVectorAsIOBanded :: (Int,Int) -> IOVector e -> IOBanded e
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

maybeViewIOBandedAsVector :: IOBanded e -> Maybe (IOVector e)
maybeViewIOBandedAsVector (IOBanded h m n kl ku f p l)
    | kl == 0 && ku == 0 =
        let c = if h == ConjTrans then Conj else NoConj
        in Just $ IOVector c (min m n) f p l
    | otherwise =
        Nothing
{-# INLINE maybeViewIOBandedAsVector #-}

bandwidthsIOBanded :: IOBanded e -> (Int,Int)
bandwidthsIOBanded a =
    (numLowerIOBanded a, numUpperIOBanded a)
{-# INLINE bandwidthsIOBanded #-}

shapeIOBanded :: IOBanded e -> (Int,Int)
shapeIOBanded a = (numRowsIOBanded a, numColsIOBanded a)
{-# INLINE shapeIOBanded #-}

boundsIOBanded :: IOBanded e -> ((Int,Int), (Int,Int))
boundsIOBanded a = ((0,0), (m-1,n-1)) where (m,n) = shapeIOBanded a
{-# INLINE boundsIOBanded #-}

sizeIOBanded :: IOBanded e -> Int
sizeIOBanded (IOBanded _ m n kl ku _ _ _) =
    foldl' (+) 0 $ map (diagLen (m,n)) [(-kl)..ku]

indicesIOBanded :: IOBanded e -> [(Int,Int)]
indicesIOBanded a =
    let is = if isHermIOBanded a
                 then [ (i,j) | i <- range (0,m-1), j <- range (0,n-1) ]
                 else [ (i,j) | j <- range (0,n-1), i <- range (0,m-1) ]
    in filter (hasStorageIOBanded a) is
  where (m,n) = shapeIOBanded a

getSizeIOBanded :: IOBanded e -> IO Int
getSizeIOBanded = return . sizeIOBanded
{-# INLINE getSizeIOBanded #-}

getIndicesIOBanded :: IOBanded e -> IO [(Int,Int)]
getIndicesIOBanded = return . indicesIOBanded
{-# INLINE getIndicesIOBanded #-}

getIndicesIOBanded' :: IOBanded e -> IO [(Int,Int)]
getIndicesIOBanded' = getIndicesIOBanded
{-# INLINE getIndicesIOBanded' #-}

getAssocsIOBanded :: IOBanded e -> IO [((Int,Int),e)]
getAssocsIOBanded a =
    let go []     = return []
        go (i:is) = unsafeInterleaveIO $ do
                        e   <- unsafeReadElemIOBanded a i
                        ies <- go is
                        return $ (i,e):ies
    in go (indicesIOBanded a)
{-# INLINE getAssocsIOBanded #-}

getAssocsIOBanded' :: IOBanded e -> IO [((Int,Int),e)]
getAssocsIOBanded' a =
    sequence [ liftM ((,) i) (unsafeReadElemIOBanded a i)
             | i <- indicesIOBanded a ]
{-# INLINE getAssocsIOBanded' #-}

getElemsIOBanded :: IOBanded e -> IO [e]
getElemsIOBanded a = liftM (snd . unzip) $ getAssocsIOBanded a
{-# INLINE getElemsIOBanded #-}

getElemsIOBanded' :: IOBanded e -> IO [e]
getElemsIOBanded' a = liftM (snd . unzip) $ getAssocsIOBanded a
{-# INLINE getElemsIOBanded' #-}

unsafeReadElemIOBanded :: IOBanded e -> (Int,Int) -> IO e
unsafeReadElemIOBanded a@(IOBanded _ _ _ _ _ _ _ _) (i,j)
    | hasStorageIOBanded a (i,j) =
        let f = if isHermIOBanded a then conjugate else id
        in liftM f $ withIOBandedElem a (i,j) peek
    | otherwise =
        return 0
{-# INLINE unsafeReadElemIOBanded #-}

unsafeWriteElemIOBanded :: IOBanded e -> (Int,Int) -> e -> IO ()
unsafeWriteElemIOBanded a@(IOBanded _ _ _ _ _ _ _ _) ij e =
    let e' = if isHermIOBanded a then conjugate e else e
    in withIOBandedElem a ij (`poke` e')
{-# INLINE unsafeWriteElemIOBanded #-}

modifyWithIOBanded :: (e -> e) -> IOBanded e -> IO ()
modifyWithIOBanded f a = do
    ies <- getAssocsIOBanded a
    forM_ ies $ \(i,e) -> do
         unsafeWriteElemIOBanded a i (f e)

canModifyElemIOBanded :: IOBanded e -> (Int,Int) -> IO Bool
canModifyElemIOBanded a ij = return $ hasStorageIOBanded a ij
{-# INLINE canModifyElemIOBanded #-}

-- | Create a new banded matrix of given shape and (lower,upper), bandwidths,
-- but do not initialize the elements.
newIOBanded_ :: (Elem e) => (Int,Int) -> (Int,Int) -> IO (IOBanded e)
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
            -> IO (IOBanded e)
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
                  -> IO (IOBanded e)
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
                 -> IO (IOBanded e)
newListsIOBanded (m,n) (kl,ku) xs = do
    a <- newIOBanded_ (m,n) (kl,ku)
    zipWithM_ (writeDiagElems a) [(negate kl)..ku] xs
    return a
  where
    writeDiagElems :: (Elem e) 
                   => IOBanded e -> Int -> [e] -> IO ()
    writeDiagElems a i es =
        let d   = unsafeDiagViewIOBanded a i
            nb  = max 0 (-i)
            es' = drop nb es
        in zipWithM_ (unsafeWriteElem d) [0..(dim d - 1)] es'

newZeroIOBanded :: (Elem e) => (Int,Int) -> (Int,Int) -> IO (IOBanded e)
newZeroIOBanded mn bw = do
    a <- newIOBanded_ mn bw
    setZeroIOBanded a
    return a

setZeroIOBanded :: IOBanded e -> IO ()
setZeroIOBanded a@(IOBanded _ _ _ _ _ _ _ _) = setConstantIOBanded 0 a
{-# INLINE setZeroIOBanded #-}

newConstantIOBanded :: (Elem e)
                    => (Int,Int) -> (Int,Int) -> e -> IO (IOBanded e)
newConstantIOBanded mn bw e = do
    a <- newIOBanded_ mn bw
    setConstantIOBanded e a
    return a

setConstantIOBanded :: e -> IOBanded e -> IO ()
setConstantIOBanded e a@(IOBanded _ _ _ _ _ _ _ _)
    | isHermIOBanded a = setConstantIOBanded (conjugate e) (hermIOBanded a)
    | otherwise = do
        is <- getIndicesIOBanded a
        mapM_ (\i -> unsafeWriteElemIOBanded a i e) is
{-# INLINE setConstantIOBanded #-}

newCopyIOBanded :: IOBanded e -> IO (IOBanded e)
newCopyIOBanded a@(IOBanded _ _ _ _ _ _ _ _)
    | isHermIOBanded a =
        liftM hermIOBanded $ newCopyIOBanded (hermIOBanded a)
    | otherwise = do
        a' <- newIOBanded_ (shapeIOBanded a) (numLowerIOBanded a, numUpperIOBanded a)
        unsafeCopyIOBanded a' a
        return a'

unsafeCopyIOBanded :: IOBanded e -> IOBanded e -> IO ()
unsafeCopyIOBanded dst@(IOBanded _ _ _ _ _ _ _ _) src
    | isHermIOBanded dst = 
        unsafeCopyIOBanded (hermIOBanded dst) 
                           (hermIOBanded src)
                         
    | (not . isHermIOBanded) src =
        withIOBanded dst $ \pDst ->
        withIOBanded src $ \pSrc ->
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


unsafeRowViewIOBanded :: IOBanded e
                      -> Int
                      -> (Int, IOVector e, Int)
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

unsafeColViewIOBanded :: IOBanded e
                      -> Int
                      -> (Int, IOVector e, Int)
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

unsafeDiagViewIOBanded :: IOBanded e -> Int -> IOVector e
unsafeDiagViewIOBanded a@(IOBanded h m n _ _ fp p ld) d
    | h == ConjTrans =
         conj $ unsafeDiagViewIOBanded (hermIOBanded a) (negate d)
    | otherwise =
        let off = indexIOBanded a (diagStart d)
            p'  = p `advancePtr` off
            len = diagLen (m,n) d
            inc = ld
        in (IOVector NoConj len fp p' inc)


unsafeGetRowIOBanded :: IOBanded e -> Int -> IO (IOVector e)
unsafeGetRowIOBanded a@(IOBanded _ _ _ _ _ _ _ _) i =
    let (nb,x,na) = unsafeRowViewIOBanded a i
        n = numCols a
    in do
        es <- getElems x
        newListVector n $ (replicate nb 0) ++ es ++ (replicate na 0)

unsafeGetColIOBanded :: IOBanded e -> Int -> IO (IOVector e)
unsafeGetColIOBanded a j =
    liftM conj $ unsafeGetRowIOBanded (hermIOBanded a) j

-- | @gbmv alpha a x beta y@ replaces @y := alpha a * x + beta y@
gbmv :: (ReadVector x IO, WriteVector y IO, BLAS2 e)
     => e -> IOBanded e -> x e -> e -> y e -> IO ()
gbmv alpha a (x :: x e) beta y
    | n == 0 =
        scaleByVector beta y
    | otherwise =
        let transA = transEnumIOBanded a
            (m',n',kl',ku') = if isHermIOBanded a
                                   then (n,m,ku,kl) else (m,n,kl,ku)
            ldA    = ldaIOBanded a
            incX   = stride x
            incY   = stride y
        in unsafePerformIOWithVector y $ \y' ->
           withIOBanded a $ \pA ->
           withIOVector (unsafeVectorToIOVector x) $ \pX ->
           withIOVector y' $ \pY ->
                BLAS.gbmv transA (conjEnum x) (conjEnum y) 
                          m' n' kl' ku' alpha pA ldA pX incX beta pY incY
  where
    (m,n)   = shape a
    (kl,ku) = bandwidthsIOBanded a


-- | @gbmm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gbmm :: (ReadMatrix b IO, WriteMatrix c IO, BLAS2 e)
     => e -> IOBanded e -> b e -> e -> c e -> IO ()
gbmm alpha a b beta c =
    sequence_ $
        zipWith (\x y -> gbmv alpha a x beta y) (colViews b) (colViews c)

unsafeGetColHermIOBanded :: (WriteVector y IO, Elem e)
                         => Herm IOBanded e -> Int -> IO (y e)
unsafeGetColHermIOBanded =
    error "TODO: unsafeGetColHermIOBanded is not implemented"

unsafeGetRowHermIOBanded :: (WriteVector y IO, Elem e)
                         => Herm IOBanded e -> Int -> IO (y e)
unsafeGetRowHermIOBanded =
    error "TODO: unsafeGetRowHermIOBanded is not implemented"

hbmv :: (ReadVector x IO, BLAS2 e) => 
    e -> Herm IOBanded e -> x e -> e -> IOVector e -> IO ()
hbmv alpha h x beta y =
    let (u,a) = hermToBase h
        n     = numCols a
        k     = case u of 
                    Upper -> numUpperIOBanded a
                    Lower -> numLowerIOBanded a      
        uplo  = case (isHermIOBanded a) of
                    True  -> flipUpLo u
                    False -> u
        ldA   = ldaIOBanded a
        incX  = stride x
        incY  = stride y
        withPtrA 
              = case uplo of Upper -> withIOBanded a
                             Lower -> withIOBandedElem a (0,0)
    in unsafePerformIOWithVector y $ \y' ->
       withPtrA $ \pA ->
       withIOVector (unsafeVectorToIOVector x) $ \pX ->
       withIOVector y' $ \pY ->
            BLAS.hbmv uplo (conjEnum x) (conjEnum y) 
                      n k alpha pA ldA pX incX beta pY incY

hbmm :: (ReadMatrix b IO, BLAS2 e) => 
    e -> Herm IOBanded e -> b e -> e -> IOMatrix e -> IO ()
hbmm alpha h b beta c =
    zipWithM_ (\x y -> hbmv alpha h x beta y) (colViews b) (colViews c)

unsafeGetColTriIOBanded :: (WriteVector y IO, Elem e)
                        => Tri IOBanded e -> Int -> IO (y e)
unsafeGetColTriIOBanded =
    error "TODO: unsafeGetColTriIOBanded is not implemented"

unsafeGetRowTriIOBanded :: (WriteVector y IO, Elem e)
                        => Tri IOBanded e -> Int -> IO (y e)
unsafeGetRowTriIOBanded =
    error "TODO: unsafeGetRowTriIOBanded is not implemented"

tbmv :: (WriteVector y IO, BLAS2 e) => 
    e -> Tri IOBanded e -> y e -> IO ()
tbmv alpha t x =
    let (u,d,a) = triToBase t
        (transA,uploA) 
                  = if isHermIOBanded a then (ConjTrans, flipUpLo u) 
                                        else (NoTrans  , u)
        diagA     = d
        n         = numCols a
        k         = case u of Upper -> numUpperIOBanded a 
                              Lower -> numLowerIOBanded a
        ldA       = ldaIOBanded a
        incX      = stride x
        withPtrA  = case uploA of 
                        Upper -> withIOBanded a
                        Lower -> withIOBandedElem a (0,0)
    in unsafePerformIOWithVector x $ \x' ->
       withPtrA $ \pA ->
       withIOVector x' $ \pX -> do
           scaleByVector alpha x
           BLAS.tbmv uploA transA diagA (conjEnum x) n k pA ldA pX incX

tbmm :: (WriteMatrix b IO, BLAS2 e) =>
    e -> Tri IOBanded e -> b e -> IO ()
tbmm 1     t b = mapM_ (\x -> tbmv 1 t x) (colViews b)
tbmm alpha t b = scaleByMatrix alpha b >> tbmm 1 t b

tbmv' :: (ReadVector x IO, WriteVector y IO, BLAS2 e) => 
    e -> Tri IOBanded e -> x e -> e -> y e -> IO ()
tbmv' alpha a (x :: x e) beta (y  :: y e)
    | beta /= 0 = do
        (x' :: y e) <- newCopyVector x
        tbmv alpha a x'
        scaleByVector beta y
        unsafeAxpyVector 1 x' y
    | otherwise = do
        unsafeCopyVector y x
        tbmv alpha a y

tbmm' :: (ReadMatrix b IO, WriteMatrix c IO, BLAS2 e) => 
    e -> Tri IOBanded e -> b e -> e -> c e -> IO ()
tbmm' alpha a (b :: b e) beta (c :: c e)
    | beta /= 0 = do
        (b' :: c e) <- newCopyMatrix b
        tbmm alpha a b'
        scaleByMatrix beta c
        unsafeAxpyMatrix 1 b' c
    | otherwise = do
        unsafeCopyMatrix c b
        tbmm alpha a c

tbsv :: (WriteVector y IO, BLAS2 e) => 
    e -> Tri IOBanded e -> y e -> IO ()
tbsv alpha t x =
    let (u,d,a) = triToBase t
        (transA,uploA) 
                  = if isHermIOBanded a then (ConjTrans, flipUpLo u) 
                                        else (NoTrans  , u)
        diagA     = d
        n         = numCols a
        k         = case u of Upper -> numUpperIOBanded a 
                              Lower -> numLowerIOBanded a
        ldA       = ldaIOBanded a
        incX      = stride x
        withPtrA  = case uploA of 
                        Upper -> withIOBanded a
                        Lower -> withIOBandedElem a (0,0)
    in unsafePerformIOWithVector x $ \x' ->
       withPtrA $ \pA ->
       withIOVector x' $ \pX -> do
           scaleByVector alpha x
           BLAS.tbsv uploA transA diagA (conjEnum x) n k pA ldA pX incX

tbsm :: (WriteMatrix b IO, BLAS2 e) => 
    e -> Tri IOBanded e -> b e -> IO ()
tbsm 1     t b = mapM_ (\x -> tbsv 1 t x) (colViews b)
tbsm alpha t b = scaleByMatrix alpha b >> tbsm 1 t b

tbsv' :: (ReadVector y IO, WriteVector x IO, BLAS2 e)
      => e -> Tri IOBanded e -> y e -> x e -> IO ()
tbsv' alpha a y x = do
    unsafeCopyVector x y
    tbsv alpha a x

tbsm' :: (ReadMatrix c IO, WriteMatrix b IO, BLAS2 e) 
      => e -> Tri IOBanded e -> c e -> b e -> IO ()
tbsm' alpha a c b = do
    unsafeCopyMatrix b c
    tbsm alpha a b

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

instance HasHerm IOBanded where
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

instance MMatrix (Herm IOBanded) IO where
    unsafeDoSApplyAddVector alpha a x beta y = 
        hbmv alpha a x beta (unsafeVectorToIOVector y)
    {-# INLINE unsafeDoSApplyAddVector #-}    
    unsafeDoSApplyAddMatrix alpha a b beta c = 
        hbmm alpha a b beta (unsafeMatrixToIOMatrix c)
    {-# INLINE unsafeDoSApplyAddMatrix #-}    
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
