{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Operations
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Operations (
    module BLAS.Matrix.Immutable,
    module BLAS.Matrix.ReadOnly,
    
    -- * Matrix Arithmetic
    -- ** Pure
    scale,
    invScale,
    
    -- ** Impure
    getScaled,
    getInvScaled,
    
    -- * In-place operations
    doConj,
    scaleBy,
    invScaleBy,

    -- * BLAS operations
    gbmv,
    gbmm,
    
    -- * Unsafe operations
    unsafeGbmv,
    unsafeGbmm
    
    ) where

import System.IO.Unsafe
import Unsafe.Coerce

import Data.Matrix.Banded.Internal
import Data.Matrix.Dense.Internal ( DMatrix, IOMatrix)
import Data.Vector.Dense.Internal hiding ( unsafeWithElemPtr, unsafeThaw, 
    unsafeFreeze )
import qualified Data.Vector.Dense.Operations as V
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Matrix.Dense.Internal as M
import qualified Data.Matrix.Dense.Operations as M

import BLAS.Access
import BLAS.C ( CBLASTrans, colMajor, noTrans, conjTrans )
import qualified BLAS.C as BLAS
import BLAS.Elem ( BLAS1, BLAS2  )
import qualified BLAS.Elem as E
import BLAS.Internal ( checkMatVecMultAdd, checkMatMatMultAdd )
import BLAS.Matrix.Immutable
import BLAS.Matrix.ReadOnly

infixl 7 `scale`, `invScale`


-- | Form a new matrix by multiplying every element by a value.
getScaled :: (BLAS1 e) => e -> BMatrix t (m,n) e -> IO (BMatrix r (m,n) e)
getScaled k = unaryOp (scaleBy k)

-- | Form a new matrix by dividing every element by a value.
getInvScaled :: (BLAS1 e) => e -> BMatrix t (m,n) e -> IO (BMatrix r (m,n) e)
getInvScaled k = unaryOp (invScaleBy k)

-- | Conjugate every element in a matrix.
doConj  :: (BLAS1 e) => IOBanded (m,n) e -> IO ()
doConj a = let (_,_,a',_) = toRawMatrix a
           in M.doConj a'

-- | Scale every element in a matrix by the given value.
scaleBy :: (BLAS1 e) => e -> IOBanded (m,n) e -> IO ()
scaleBy k a = 
    let (_,_,a',h) = toRawMatrix a
        k' = if h then E.conj k else k
    in M.scaleBy k' a'
    
-- | Divide every element by the given value.
invScaleBy :: (BLAS1 e) => e -> IOBanded (m,n) e -> IO ()
invScaleBy k a = 
    let (_,_,a',h) = toRawMatrix a
        k' = if h then E.conj k else k
    in M.invScaleBy k' a'

blasTransOf :: BMatrix t (m,n) e -> CBLASTrans
blasTransOf a = 
    case (isHerm a) of
          False -> noTrans
          True  -> conjTrans

flipShape :: (Int,Int) -> (Int,Int)
flipShape (m,n) = (n,m)

-- | @gemv alpha a x beta y@ replaces @y := alpha a * x + beta y@
gbmv :: (BLAS2 e) => e -> BMatrix s (m,n) e -> DVector t n e -> e -> IOVector m e -> IO ()
gbmv alpha a x beta y =
    checkMatVecMultAdd (shape a) (dim x) (dim y) $ 
        unsafeGbmv alpha a x beta y

unsafeGbmv :: (BLAS2 e) => e -> BMatrix s (m,n) e -> DVector t n e -> e -> IOVector m e -> IO ()
unsafeGbmv alpha a x beta y
    | numRows a == 0 || numCols a == 0 =
        return ()
    | isConj x = do
        x' <- V.getConj (conj x)
        gbmv alpha a x' beta y
    | isConj y = do
        V.doConj y
        gbmv alpha a x beta (conj y)
        V.doConj y
    | otherwise =
        let order  = colMajor
            transA = blasTransOf a
            (m,n)  = case (isHerm a) of
                         False -> shape a
                         True  -> (flipShape . shape) a
            (kl,ku) = case (isHerm a) of
                          False -> (numLower a, numUpper a)
                          True  -> (numUpper a, numLower a)
            ldA    = ldaOf a
            incX   = V.strideOf x
            incY   = V.strideOf y
        in unsafeWithBasePtr a $ \pA ->
               V.unsafeWithElemPtr x 0 $ \pX ->
                    V.unsafeWithElemPtr y 0 $ \pY -> do
                        BLAS.gbmv order transA m n kl ku alpha pA ldA pX incX beta pY incY

-- | @gemm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gbmm :: (BLAS2 e) => e -> BMatrix s (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> IO ()
gbmm alpha a b beta c =
    checkMatMatMultAdd (shape a) (shape b) (shape c) $ 
        unsafeGbmm alpha a b beta c

unsafeGbmm :: (BLAS2 e) => e -> BMatrix s (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> IO ()
unsafeGbmm alpha a b beta c =
    sequence_ $
        zipWith (\x y -> unsafeGbmv alpha a x beta y) (M.cols b) (M.cols c)

unaryOp :: (BLAS1 e) => (IOBanded (m,n) e -> IO ()) 
    -> BMatrix t (m,n) e -> IO (BMatrix r (m,n) e)
unaryOp f a = do
    a' <- newCopy a
    f (unsafeThaw a')
    return (unsafeCoerce a')

-- | Create a new matrix by scaling another matrix by the given value.
scale :: (BLAS1 e) => e -> Banded (m,n) e -> Banded (m,n) e
scale k a = unsafePerformIO $ getScaled k a
{-# NOINLINE scale #-}

-- | Form a new matrix by dividing every element by a value.
invScale :: (BLAS1 e) => e -> Banded (m,n) e -> Banded (m,n) e
invScale k a = unsafePerformIO $ getInvScaled k a
{-# NOINLINE invScale #-}

instance (BLAS2 e) => RMatrix (BMatrix s) e where
    unsafeGetSApply alpha a x = do
        y <- V.newZero (numRows a)
        unsafeGbmv alpha a x 0 y
        return (unsafeCoerce y)

    unsafeGetSApplyMat alpha a b = do
        c <- newZero (numRows a, numCols b)
        unsafeGbmm alpha a b 0 c
        return (unsafeCoerce c)

instance (BLAS2 e) => IMatrix (BMatrix Imm) e where
    unsafeSApply alpha a x = unsafePerformIO $ unsafeGetSApply alpha a x
    {-# NOINLINE unsafeSApply #-}

    unsafeSApplyMat alpha a b = unsafePerformIO $ unsafeGetSApplyMat alpha a b
    {-# NOINLINE unsafeSApplyMat #-}

