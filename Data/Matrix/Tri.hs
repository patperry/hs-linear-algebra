{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Tri
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Tri (
    Tri(..),
    UpLo(..), Diag(..),

    fromBase,
    toBase,
    mapTri,

    lower,
    lowerFat,
    lowerTall,
    
    lowerU,
    lowerUFat,
    lowerUTall,
    
    upper,
    upperFat,
    upperTall,
    
    upperU,
    upperUFat,
    upperUTall,

    coerceTri,
    
    module BLAS.Matrix.Apply,
    module BLAS.Matrix.Solve,
    ) where

import BLAS.Matrix.Apply
import BLAS.Matrix.Solve


import Control.Monad( when )
import Control.Monad.ST( ST )
import Unsafe.Coerce

import BLAS.C( BLAS3 )
import BLAS.Internal ( checkSquare, checkFat, checkTall )
import BLAS.UnsafeIOToM
import BLAS.UnsafeInterleaveM
import BLAS.Matrix hiding ( diag )
--import BLAS.Tensor
import BLAS.Types ( UpLo(..), Diag(..), Trans(..), flipTrans, flipUpLo )
import BLAS.C.Types ( cblasDiag, cblasUpLo, cblasTrans, colMajor, 
    noTrans, conjTrans, leftSide, rightSide )
import qualified BLAS.C as BLAS
import qualified BLAS.C.Types as BLAS

import Data.Matrix.Dense.Class hiding ( BaseMatrix )
import Data.Matrix.Dense( Matrix )
import Data.Matrix.Dense.IO( IOMatrix )
import Data.Matrix.Dense.ST( STMatrix, runSTMatrix )
import qualified Data.Matrix.Dense.Class as Dense
import Data.Vector.Dense.Class
import Data.Vector.Dense.ST( runSTVector )

data Tri a mn e = Tri UpLo Diag (a mn e)

-- | Coerce the shape type.
coerceTri :: Tri a mn e -> Tri a mn' e
coerceTri = unsafeCoerce

mapTri :: (a (m,n) e -> b (m,n) e) -> Tri a (m,n) e -> Tri b (m,n) e
mapTri f (Tri u d a) = Tri u d $ f a

fromBase :: UpLo -> Diag -> a (m,n) e -> Tri a (m,n) e
fromBase = Tri
        
toBase :: Tri a (m,n) e -> (UpLo, Diag, a (m,n) e)
toBase (Tri u d a) = (u,d,a)


lower :: (BaseMatrix a e) => a (n,n) e -> Tri a (n,n) e
lower a = checkSquare (shape a) $ Tri Lower NonUnit a

lowerFat :: (BaseMatrix a e) => a (m,n) e -> Tri a (m,m) e
lowerFat a = checkFat (shape a) $ Tri Lower NonUnit (unsafeCoerce a)

lowerTall :: (BaseMatrix a e) => a (m,n) e -> Tri a (m,n) e
lowerTall a = checkTall (shape a) $ Tri Lower NonUnit a


lowerU :: (BaseMatrix a e) => a (n,n) e -> Tri a (n,n) e
lowerU a = checkSquare (shape a) $ Tri Lower Unit a

lowerUFat :: (BaseMatrix a e) => a (m,n) e -> Tri a (m,m) e
lowerUFat a = checkFat (shape a) $ Tri Lower Unit (unsafeCoerce a)

lowerUTall :: (BaseMatrix a e) => a (m,n) e -> Tri a (m,n) e
lowerUTall a = checkTall (shape a) $ Tri Lower Unit a


upper :: (BaseMatrix a e) => a (n,n) e -> Tri a (n,n) e
upper a = checkSquare (shape a) $ Tri Upper NonUnit a

upperFat :: (BaseMatrix a e) => a (m,n) e -> Tri a (m,n) e
upperFat a = checkFat (shape a) $ Tri Upper NonUnit a

upperTall :: (BaseMatrix a e) => a (m,n) e -> Tri a (n,n) e
upperTall a = checkTall (shape a) $ Tri Upper NonUnit (unsafeCoerce a)


upperU :: (BaseMatrix a e) => a (n,n) e -> Tri a (n,n) e
upperU a = checkSquare (shape a) $ Tri Upper Unit a

upperUFat :: (BaseMatrix a e) => a (m,n) e -> Tri a (m,n) e
upperUFat a = checkFat (shape a) $ Tri Upper Unit a

upperUTall :: (BaseMatrix a e) => a (m,n) e -> Tri a (n,n) e
upperUTall a = checkTall (shape a) $ Tri Upper Unit (unsafeCoerce a)

      
instance BaseMatrix a e => BaseTensor (Tri a) (Int,Int) e where
    shape (Tri Lower _ a) = (numRows a, min (numRows a) (numCols a))
    shape (Tri Upper _ a) = (min (numRows a) (numCols a), numCols a)
    
    bounds a = ((0,0),(m-1,n-1)) where (m,n) = shape a
    
instance BaseMatrix a e => BaseMatrix (Tri a) e where
    herm (Tri u d a) = Tri (flipUpLo u) d (herm a)


instance (Show (a (m,n) e), BaseMatrix a e) => Show (Tri a (m,n) e) where
    show (Tri u d a) =
        constructor ++ suffix ++ " (" ++ show a ++ ")"
        where
          constructor = case (u,d) of
              (Lower, NonUnit) -> "lower"
              (Lower, Unit   ) -> "lowerU"
              (Upper, NonUnit) -> "upper"
              (Upper, Unit   ) -> "upperU"

          suffix = case undefined of
                       _ | isSquare a -> ""
                       _ | isFat a    -> "Fat"
                       _              -> "Tall"

------------------------ Tri Matrix Apply Functions -------------------------

instance (BLAS3 e) => IApply (Tri Matrix) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    

instance (BLAS3 e) => ReadApply (Tri IOMatrix) e IO where
    unsafeDoSApplyAdd    = unsafeDoSApplyAddTriMatrix
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatTriMatrix
    unsafeDoSApply_      = trmv
    unsafeDoSApplyMat_   = trmm

instance (BLAS3 e) => ReadApply (Tri (STMatrix s)) e (ST s) where
    unsafeDoSApplyAdd    = unsafeDoSApplyAddTriMatrix
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatTriMatrix
    unsafeDoSApply_      = trmv
    unsafeDoSApplyMat_   = trmm

instance (BLAS3 e, UnsafeIOToM m, UnsafeInterleaveM m) => ReadApply (Tri Matrix) e m where
    unsafeDoSApplyAdd    = unsafeDoSApplyAddTriMatrix
    unsafeDoSApplyAddMat = unsafeDoSApplyAddMatTriMatrix
    unsafeDoSApply_      = trmv
    unsafeDoSApplyMat_   = trmm



unsafeDoSApplyAddTriMatrix :: (BLAS3 e, ReadMatrix a z e m, ReadApply a e m, 
    ReadVector x e m, WriteVector y e m) =>
        e -> Tri a (k,l) e -> x l e -> e -> y k e -> m ()
unsafeDoSApplyAddTriMatrix alpha t x beta y =
    if beta == 0
        then unsafeDoSApplyTriMatrix alpha t x y
        else do
            y' <- newCopyVector y
            unsafeDoSApplyTriMatrix alpha t x y'
            scaleBy beta y
            unsafeAxpyVector 1 y' y

unsafeDoSApplyAddMatTriMatrix :: (BLAS3 e, ReadMatrix a z e m, ReadApply a e m, 
    ReadMatrix b x e m, WriteMatrix c y e m) =>
        e -> Tri a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
unsafeDoSApplyAddMatTriMatrix alpha t b beta c =
    if beta == 0
        then unsafeDoSApplyMatTriMatrix alpha t b c
        else do
            c' <- newCopyMatrix c
            unsafeDoSApplyMatTriMatrix alpha t b c'
            scaleBy beta c
            unsafeAxpyMatrix 1 c' c

unsafeDoSApplyTriMatrix :: (BLAS3 e, ReadMatrix a z e m, ReadApply a e m, 
    ReadVector x e m, WriteVector y e m) =>
        e -> Tri a (k,l) e -> x l e -> y k e -> m ()
unsafeDoSApplyTriMatrix alpha t x y =
    case (u, toLower d a, toUpper d a) of
        (Lower,Left t',_) -> do
            unsafeCopyVector y (coerceVector x)
            trmv alpha t' y
            
        (Lower,Right (t',r),_) -> do
            let y1 = unsafeSubvector y 0            (numRows t')
                y2 = unsafeSubvector y (numRows t') (numRows r)
            unsafeCopyVector y1 x
            trmv alpha t' y1
            unsafeDoSApplyAdd alpha r x 0 y2
            
        (Upper,_,Left t') -> do
            unsafeCopyVector (coerceVector y) x
            trmv alpha t' (coerceVector y)

        (Upper,_,Right (t',r)) ->
            let x1 = unsafeSubvector x 0            (numCols t')
                x2 = unsafeSubvector x (numCols t') (numCols r)
            in do
                unsafeCopyVector y x1
                trmv alpha t' y
                unsafeDoSApplyAdd alpha r x2 1 y
  where
    (u,d,a) = toBase t

unsafeDoSApplyMatTriMatrix :: (BLAS3 e, ReadMatrix a z e m, ReadApply a e m, 
    ReadMatrix b x e m, WriteMatrix c y e m) =>
        e -> Tri a (r,s) e -> b (s,t) e -> c (r,t) e -> m ()
unsafeDoSApplyMatTriMatrix alpha t b c =
    case (u, toLower d a, toUpper d a) of
        (Lower,Left t',_) -> do
            unsafeCopyMatrix c (coerceMatrix b)
            trmm alpha t' c
            
        (Lower,Right (t',r),_) -> do
            let c1 = unsafeSubmatrix c (0,0)          (numRows t',numCols c)
                c2 = unsafeSubmatrix c (numRows t',0) (numRows r ,numCols c)
            unsafeCopyMatrix c1 b
            trmm alpha t' c1
            unsafeDoSApplyAddMat alpha r b 0 c2
            
        (Upper,_,Left t') -> do
            unsafeCopyMatrix (coerceMatrix c) b
            trmm alpha t' (coerceMatrix c)

        (Upper,_,Right (t',r)) ->
            let b1 = unsafeSubmatrix b (0,0)          (numCols t',numCols b)
                b2 = unsafeSubmatrix b (numCols t',0) (numCols r ,numCols b)
            in do
                unsafeCopyMatrix c b1
                trmm alpha t' c
                unsafeDoSApplyAddMat alpha r b2 1 c
  where
    (u,d,a) = toBase t


toLower :: (Dense.BaseMatrix a x e) => Diag -> a (m,n) e 
        -> Either (Tri a (m,m) e) 
                  (Tri a (n,n) e, a (d,n) e)
toLower diag a =
    if m <= n
        then Left $  fromBase Lower diag (unsafeSubmatrix a (0,0) (m,m))
        else let t = fromBase Lower diag (unsafeSubmatrix a (0,0) (n,n))
                 r = unsafeSubmatrix a (n,0) (d,n)
             in Right (t,r)
  where
    (m,n) = shape a
    d     = m - n
    
toUpper :: (Dense.BaseMatrix a x e) => Diag -> a (m,n) e
        -> Either (Tri a (n,n) e)
                  (Tri a (m,m) e, a (m,d) e)
toUpper diag a =
    if n <= m
        then Left $  fromBase Upper diag (unsafeSubmatrix a (0,0) (n,n))
        else let t = fromBase Upper diag (unsafeSubmatrix a (0,0) (m,m))
                 r = unsafeSubmatrix a (0,m) (m,d)
             in Right (t,r)
  where
    (m,n) = shape a
    d     = n - m

trmv :: (ReadMatrix a x e m, WriteVector y e m, BLAS3 e) => 
    e -> Tri a (k,k) e -> y n e -> m ()
trmv alpha t x 
    | dim x == 0 = 
        return ()
        
    | isConj x =
        let (u,d,a) = toBase t
            order   = colMajor
            side    = rightSide
            (h,u')  = if isHerm a then (NoTrans, flipUpLo u) else (ConjTrans, u)
            uploA   = cblasUpLo u'
            transA  = cblasTrans h
            diagA   = cblasDiag d
            m       = 1
            n       = dim x
            alpha'  = conj alpha
            ldA     = lda a
            ldB     = stride x
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withVectorPtr x $ \pB ->
                   BLAS.trmm order side uploA transA diagA m n alpha' pA ldA pB ldB

    | otherwise =
        let (u,d,a)   = toBase t
            order     = colMajor
            (transA,u') = if isHerm a then (conjTrans, flipUpLo u) else (noTrans, u)
            uploA     = cblasUpLo u'
            diagA     = cblasDiag d
            n         = dim x
            ldA       = lda a
            incX      = stride x
        in do
            when (alpha /= 1) $ scaleBy alpha x
            unsafeIOToM $
                withMatrixPtr a $ \pA ->
                withVectorPtr x $ \pX -> do
                   BLAS.trmv order uploA transA diagA n pA ldA pX incX


trmm :: (ReadMatrix a x e m, WriteMatrix b y e m, BLAS3 e) => 
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
trmm _ _ b
    | numRows b == 0 || numCols b == 0 = return ()
trmm alpha t b =
    let (u,d,a)   = toBase t
        order     = colMajor
        (h,u')    = if isHerm a then (ConjTrans, flipUpLo u) else (NoTrans, u)
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if isHerm b
                        then (rightSide, flipTrans h, n, m, conj alpha)
                        else (leftSide , h          , m, n, alpha       )
        uploA     = cblasUpLo u'
        transA    = cblasTrans h'
        diagA     = cblasDiag d
        ldA       = lda a
        ldB       = lda b
    in unsafeIOToM $
           withMatrixPtr a $ \pA ->
           withMatrixPtr b $ \pB ->
               BLAS.trmm order side uploA transA diagA m' n' alpha' pA ldA pB ldB


------------------------ Tri Matrix Solve Functions -------------------------

instance (BLAS3 e) => ISolve (Tri Matrix) e where
    unsafeSSolve    alpha a y = runSTVector $ unsafeGetSSolve    alpha a y
    unsafeSSolveMat alpha a c = runSTMatrix $ unsafeGetSSolveMat alpha a c

instance (BLAS3 e) => ReadSolve (Tri IOMatrix) e IO where
    unsafeDoSSolve     = unsafeDoSSolveTriMatrix
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriMatrix
    unsafeDoSSolve_    = trsv
    unsafeDoSSolveMat_ = trsm

instance (BLAS3 e) => ReadSolve (Tri (STMatrix s)) e (ST s) where
    unsafeDoSSolve     = unsafeDoSSolveTriMatrix
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriMatrix
    unsafeDoSSolve_    = trsv
    unsafeDoSSolveMat_ = trsm

instance (BLAS3 e, UnsafeIOToM m, UnsafeInterleaveM m) => ReadSolve (Tri Matrix) e m where
    unsafeDoSSolve     = unsafeDoSSolveTriMatrix
    unsafeDoSSolveMat  = unsafeDoSSolveMatTriMatrix
    unsafeDoSSolve_    = trsv
    unsafeDoSSolveMat_ = trsm



unsafeDoSSolveTriMatrix :: (BLAS3 e, ReadMatrix a z e m, ReadApply a e m, 
    ReadVector y e m, WriteVector x e m) =>
        e -> Tri a (k,l) e -> y k e -> x l e -> m ()
unsafeDoSSolveTriMatrix alpha t y x =
    case (u, toLower d a, toUpper d a) of
        (Lower,Left t',_) -> do
            unsafeCopyVector x (coerceVector y)
            trsv alpha t' (coerceVector x)
            
        (Lower,Right (t',_),_) -> do
            let y1 = unsafeSubvector y 0            (numRows t')
            unsafeCopyVector x y1
            trsv alpha t' x
            
        (Upper,_,Left t') -> do
            unsafeCopyVector x (coerceVector y)
            trsv alpha t' x

        (Upper,_,Right (t',r)) ->
            let x1 = unsafeSubvector x 0            (numCols t')
                x2 = unsafeSubvector x (numCols t') (numCols r)
            in do
                unsafeCopyVector x1 y
                trsv alpha t' x1
                setZero x2
  where
    (u,d,a) = toBase t


unsafeDoSSolveMatTriMatrix :: (BLAS3 e, ReadMatrix a z e m, ReadApply a e m, 
    ReadMatrix c y e m, WriteMatrix b x e m) =>
        e -> Tri a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
unsafeDoSSolveMatTriMatrix alpha t c b =
    case (u, toLower d a, toUpper d a) of
        (Lower,Left t',_) -> do
            unsafeCopyMatrix b (coerceMatrix c)
            trsm alpha t' (coerceMatrix b)
            
        (Lower,Right (t',_),_) -> do
            let c1 = unsafeSubmatrix c (0,0)          (numRows t',numCols c)
            unsafeCopyMatrix b c1
            trsm alpha t' b
            
        (Upper,_,Left t') -> do
            unsafeCopyMatrix (coerceMatrix b) c
            trsm alpha t' (coerceMatrix b)

        (Upper,_,Right (t',r)) ->
            let b1 = unsafeSubmatrix b (0,0)          (numCols t',numCols b)
                b2 = unsafeSubmatrix b (numCols t',0) (numCols r ,numCols b)
            in do
                unsafeCopyMatrix b1 c
                trsm alpha t' b1
                setZero b2
  where
    (u,d,a) = toBase t


trsv :: (ReadMatrix a x e m, WriteVector y e m, BLAS3 e) => 
    e -> Tri a (k,k) e -> y n e -> m ()
trsv alpha t x
    | dim x == 0 = return ()

    | isConj x =
        let (u,d,a) = toBase t
            order   = colMajor
            side    = rightSide
            (h,u')  = if isHerm a then (NoTrans, flipUpLo u) else (ConjTrans, u)
            uploA   = cblasUpLo u'
            transA  = cblasTrans h
            diagA   = cblasDiag d
            m       = 1
            n       = dim x
            alpha'  = conj alpha
            ldA     = lda a
            ldB     = stride x
        in unsafeIOToM $
               withMatrixPtr a $ \pA ->
               withVectorPtr x $ \pB ->
                   BLAS.trsm order side uploA transA diagA m n alpha' pA ldA pB ldB

    | otherwise =
        let (u,d,a) = toBase t
            order     = colMajor
            (transA,u') = if isHerm a then (conjTrans, flipUpLo u) else (noTrans, u)
            uploA     = cblasUpLo u'
            diagA     = cblasDiag d
            n         = dim x
            ldA       = lda a
            incX      = stride x
        in do
            when (alpha /= 1) $ scaleBy alpha x
            unsafeIOToM $
                withMatrixPtr a $ \pA ->
                withVectorPtr x $ \pX ->
                    BLAS.trsv order uploA transA diagA n pA ldA pX incX

trsm :: (ReadMatrix a x e m, WriteMatrix b y e m, BLAS3 e) => 
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
trsm _ _ b
    | numRows b == 0 || numCols b == 0 = return ()
trsm alpha t b =
    let (u,d,a)   = toBase t
        order     = colMajor
        (h,u')    = if isHerm a then (ConjTrans, flipUpLo u) else (NoTrans, u)
        (m,n)     = shape b
        (side,h',m',n',alpha')
                  = if isHerm b
                        then (rightSide, flipTrans h, n, m, conj alpha)
                        else (leftSide , h          , m, n, alpha     )
        uploA     = cblasUpLo u'
        transA    = cblasTrans h'
        diagA     = cblasDiag d
        ldA       = lda a
        ldB       = lda b
    in unsafeIOToM $     
           withMatrixPtr a $ \pA ->
           withMatrixPtr b $ \pB -> do
               BLAS.trsm order side uploA transA diagA m' n' alpha' pA ldA pB ldB
