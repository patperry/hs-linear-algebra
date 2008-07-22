{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Operations
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Operations (
    -- * Matrix multiplication
    -- ** Pure
    apply,
    applyMat,
    sapply,
    sapplyMat,
    
    -- ** Impure
    getApply,
    getApplyMat,
    getSApply,
    getSApplyMat,

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
    
    ) where

-- import Data.Maybe ( fromJust )
import System.IO.Unsafe
import Unsafe.Coerce

import BLAS.Internal ( checkMatVecMult, checkMatMatMult )
import Data.Matrix.Banded.Internal
import Data.Matrix.Dense.Internal ( DMatrix(..), Matrix, IOMatrix )
import Data.Vector.Dense.Internal hiding ( unsafeWithElemPtr, unsafeThaw, 
    unsafeFreeze )
-- import qualified Data.Vector.Dense.Operations as V
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Matrix.Dense.Operations as M

-- import BLAS.C ( CBLASTrans, colMajor, noTrans, conjTrans )
-- import qualified BLAS.C as BLAS
import BLAS.Elem ( BLAS1, BLAS2  )
import qualified BLAS.Elem as E

infixl 7 `apply`, `applyMat`, `scale`, `invScale`

-- | Multiply a matrix by a vector.
getApply :: (BLAS2 e) => BMatrix s (m,n) e -> DVector t n e -> IO (DVector r m e)
getApply = getSApply 1

-- | Multiply a scaled matrix by a vector.
getSApply :: (BLAS2 e) => e -> BMatrix s (m,n) e -> DVector t n e -> IO (DVector r m e)
getSApply alpha a x = checkMatVecMult (shape a) (V.dim x) >> unsafeGetSApply alpha a x

unsafeGetSApply :: (BLAS2 e) => e -> BMatrix s (m,n) e -> DVector t n e -> IO (DVector r m e)
unsafeGetSApply alpha a x = do
    y <- V.newZero (numRows a)
    gbmv alpha a x 0 y
    return (unsafeCoerce y)

-- | Multiply a matrix by a matrix.
getApplyMat :: (BLAS2 e) => BMatrix s (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
getApplyMat = getSApplyMat 1

-- | Multiply a scaled matrix by a matrix.
getSApplyMat :: (BLAS2 e) => e -> BMatrix s (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
getSApplyMat alpha a b = checkMatMatMult (shape a) (shape b) >> unsafeGetSApplyMat alpha a b

unsafeGetSApplyMat :: (BLAS2 e) => e -> BMatrix s (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
unsafeGetSApplyMat alpha a b = do
    c <- newZero (numRows a, numCols b)
    gbmm alpha a b 0 c
    return (unsafeCoerce c)

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

{-
blasTransOf :: BMatrix t (m,n) e -> CBLASTrans
blasTransOf a = 
    case (isHerm a) of
          False -> noTrans
          True  -> conjTrans

flipShape :: (Int,Int) -> (Int,Int)
flipShape (m,n) = (n,m)
-}

-- | @gemv alpha a x beta y@ replaces @y := alpha a * x + beta y@
gbmv :: (BLAS2 e) => e -> BMatrix s (m,n) e -> DVector t n e -> e -> IOVector m e -> IO ()
gbmv = undefined -- alpha a x beta y =
{-
    | numRows a == 0 || numCols a == 0 =
        return ()
    | isConj x =
        let b = fromJust $ maybeFromCol x
        in case maybeFromCol y of
               Just c  -> gemm alpha a b beta c
               Nothing -> do
                   V.doConj y
                   gemm alpha a b beta (fromJust $ maybeFromCol $ V.conj y)
                   V.doConj y
    | isConj y = 
        let c = fromJust $ maybeFromCol y
        in case maybeFromCol x of
               Just b  -> gemm alpha a b beta c
               Nothing -> do
                   x' <- V.newCopy (V.unsafeThaw x)
                   V.doConj x'
                   gemm alpha a (fromJust $ maybeFromCol $ V.conj x') beta c
    | otherwise =
        let order  = colMajor
            transA = blasTransOf a
            (m,n)  = case (isHerm a) of
                         False -> shape a
                         True  -> (flipShape . shape) a
            ldA    = ldaOf a
            incX   = V.strideOf x
            incY   = V.strideOf y
        in unsafeWithElemPtr a (0,0) $ \pA ->
               V.unsafeWithElemPtr x 0 $ \pX ->
                    V.unsafeWithElemPtr y 0 $ \pY -> do
                        BLAS.gemv order transA m n alpha pA ldA pX incX beta pY incY
-}

-- | @gemm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gbmm :: (BLAS2 e) => e -> BMatrix s (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> IO ()
gbmm = undefined -- alpha a b beta c
{-
    | numRows a == 0 || numCols a == 0 || numCols b == 0 = return ()
    | isHerm c = gemm (E.conj alpha) (herm b) (herm a) (E.conj beta) (herm c)
    | otherwise =
        let order  = colMajor
            transA = blasTransOf a
            transB = blasTransOf b
            (m,n)  = shape c
            k      = numCols a
            ldA    = ldaOf a
            ldB    = ldaOf b
            ldC    = ldaOf c
        in unsafeWithElemPtr a (0,0) $ \pA ->
               unsafeWithElemPtr b (0,0) $ \pB ->
                   unsafeWithElemPtr c (0,0) $ \pC ->
                       BLAS.gemm order transA transB m n k alpha pA ldA pB ldB beta pC ldC
-}

unaryOp :: (BLAS1 e) => (IOBanded (m,n) e -> IO ()) 
    -> BMatrix t (m,n) e -> IO (BMatrix r (m,n) e)
unaryOp f a = do
    a' <- newCopy a
    f (unsafeThaw a')
    return (unsafeCoerce a')

{-
binaryOp :: (BLAS1 e) => String -> (IOMatrix (m,n) e -> DMatrix t (m,n) e -> IO ()) 
    -> DMatrix s (m,n) e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
binaryOp name f a b = 
    checkMatMatOp name (shape a) (shape b) >> do
        a' <- newCopy a
        f (unsafeThaw a') b
        return (unsafeCoerce a')
-}        
        
-- | Multiply a matrix by a vector.
apply :: (BLAS2 e) => Banded (m,n) e -> Vector n e -> Vector m e
apply = sapply 1

-- | Multiply a scaled matrix by a vector.
sapply :: (BLAS2 e) => e -> Banded (m,n) e -> Vector n e -> Vector m e
sapply alpha a x = unsafePerformIO $ getSApply alpha a x
{-# NOINLINE sapply #-}

-- | Multiply a scaled matrix by a matrix.
sapplyMat :: (BLAS2 e) => e -> Banded (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
sapplyMat alpha a b = unsafePerformIO $ getSApplyMat alpha a b
{-# NOINLINE sapplyMat #-}
    
-- | Multiply a matrix by a matrix.
applyMat :: (BLAS2 e) => Banded (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
applyMat = sapplyMat 1

-- | Create a new matrix by scaling another matrix by the given value.
scale :: (BLAS1 e) => e -> Banded (m,n) e -> Banded (m,n) e
scale k a = unsafePerformIO $ getScaled k a
{-# NOINLINE scale #-}

-- | Form a new matrix by dividing every element by a value.
invScale :: (BLAS1 e) => e -> Banded (m,n) e -> Banded (m,n) e
invScale k a = unsafePerformIO $ getInvScaled k a
{-# NOINLINE invScale #-}

{-# RULES
"scale.apply/sapply"       forall k a x. apply (scale k a) x = sapply k a x
"scale.applyMat/sapplyMat" forall k a b. applyMat (scale k a) b = sapplyMat k a b
  #-}
  