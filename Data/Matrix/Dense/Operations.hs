{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Operations
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Operations (
    -- * Copy and swap
    copyMatrix,
    swapMatrices,

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
    shift,
    scale,
    invScale,
    add,
    plus,
    minus,
    times,
    divide,
    
    -- ** Impure
    getShifted,
    getScaled,
    getInvScaled,
    getSum,
    getDiff,
    getProduct,
    getRatio,
    
    -- * In-place operations
    doConj,
    shiftBy,
    scaleBy,
    invScaleBy,
    axpy,
    (+=),
    (-=),
    (*=),
    (//=),
    
    -- * BLAS operations
    gemv,
    gemm,
    
    ) where

import Data.Maybe ( fromJust )
import System.IO.Unsafe
import Unsafe.Coerce

import BLAS.Internal ( checkMatMatOp, checkMatVecMult, checkMatMatMult )
import Data.Matrix.Dense.Internal
import Data.Vector.Dense.Internal hiding ( unsafeWithElemPtr, unsafeThaw, 
    unsafeFreeze )
import qualified Data.Vector.Dense.Operations as V
import qualified Data.Vector.Dense.Internal as V

import BLAS.C ( CBLASTrans, colMajor, noTrans, conjTrans )
import qualified BLAS.C as BLAS
import BLAS.Elem ( BLAS1, BLAS2, BLAS3 )
import qualified BLAS.Elem as E

infixl 7 `apply`, `applyMat`, `scale`, `invScale`
infixl 6 `shift`
infixl 1 +=, -=, *=, //=

-- | @copy dst src@ copies the elements from the second argument to the first.
copyMatrix :: (BLAS1 e) => IOMatrix (m,n) e -> DMatrix t (m,n) e -> IO ()
copyMatrix a b = checkMatMatOp "copyMatrix" (shape a) (shape b) >> unsafeCopyMatrix a b

unsafeCopyMatrix :: (BLAS1 e) => IOMatrix (m,n) e -> DMatrix t (m,n) e -> IO ()
unsafeCopyMatrix = liftV2 (V.unsafeCopyVector)

-- | @swap a b@ exchanges the elements stored in two matrices.
swapMatrices :: (BLAS1 e) => IOMatrix (m,n) e -> IOMatrix (m,n) e -> IO ()
swapMatrices a b = checkMatMatOp "swapMatrices" (shape a) (shape b) >> unsafeSwapMatrices a b

unsafeSwapMatrices :: (BLAS1 e) => IOMatrix (m,n) e -> IOMatrix (m,n) e -> IO ()
unsafeSwapMatrices = liftV2 (V.unsafeSwapVectors)

-- | Multiply a matrix by a vector.
getApply :: (BLAS3 e) => DMatrix s (m,n) e -> DVector t n e -> IO (DVector r m e)
getApply = getSApply 1

-- | Multiply a scaled matrix by a vector.
getSApply :: (BLAS3 e) => e -> DMatrix s (m,n) e -> DVector t n e -> IO (DVector r m e)
getSApply alpha a x = checkMatVecMult (shape a) (V.dim x) >> unsafeGetSApply alpha a x

unsafeGetSApply :: (BLAS3 e) => e -> DMatrix s (m,n) e -> DVector t n e -> IO (DVector r m e)
unsafeGetSApply alpha a x = do
    y <- V.newZero (numRows a)
    gemv alpha a x 0 y
    return (unsafeCoerce y)

-- | Multiply a matrix by a matrix.
getApplyMat :: (BLAS3 e) => DMatrix s (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
getApplyMat = getSApplyMat 1

-- | Multiply a scaled matrix by a matrix.
getSApplyMat :: (BLAS3 e) => e -> DMatrix s (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
getSApplyMat alpha a b = checkMatMatMult (shape a) (shape b) >> unsafeGetSApplyMat alpha a b

unsafeGetSApplyMat :: (BLAS3 e) => e -> DMatrix s (m,k) e -> DMatrix t (k,n) e -> IO (DMatrix r (m,n) e)
unsafeGetSApplyMat alpha a b = do
    c <- newZero (numRows a, numCols b)
    gemm alpha a b 0 c
    return (unsafeCoerce c)

-- | Form a new matrix by adding a value to every element in a matrix.
getShifted :: (BLAS1 e) => e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
getShifted k = unaryOp (shiftBy k)

-- | Form a new matrix by multiplying every element by a value.
getScaled :: (BLAS1 e) => e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
getScaled k = unaryOp (scaleBy k)

-- | Form a new matrix by dividing every element by a value.
getInvScaled :: (BLAS1 e) => e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
getInvScaled k = unaryOp (invScaleBy k)

-- | Create a new matrix by taking the elementwise sum of two matrices.
getSum :: (BLAS1 e) => e -> DMatrix s (m,n) e -> e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
getSum alpha a beta b = checkMatMatOp "getSum" (shape a) (shape b) 
    >> unsafeGetSum alpha a beta b

unsafeGetSum :: (BLAS1 e) => e -> DMatrix s (m,n) e -> e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
unsafeGetSum alpha a@(H _) beta b = do
    s <- unsafeGetSum (E.conj alpha) (herm a) (E.conj beta) (herm b)
    return (herm s)
unsafeGetSum alpha a@(DM _ _ _ _ _) beta b = do
    s <- getScaled alpha a
    axpy beta b (unsafeThaw s)
    return (unsafeCoerce s)

-- | Create a new matrix by taking the elementwise difference of two matrices.
getDiff :: (BLAS1 e) => DMatrix s (m,n) e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
getDiff a b = checkMatMatOp "minus" (shape a) (shape b) >> unsafeGetSum 1 a (-1) b

-- | Create a new matrix by taking the elementwise product of two matrices.
getProduct :: (BLAS2 e) => DMatrix s (m,n) e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
getProduct = binaryOp "times" (*=)

-- | Create a new matrix by taking the elementwise ratio of two matrices.
getRatio :: (BLAS2 e) => DMatrix s (m,n) e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
getRatio = binaryOp "getRatio" (//=)

-- | Conjugate every element in a matrix.
doConj  :: (BLAS1 e) => IOMatrix (m,n) e -> IO (IOMatrix (m,n) e)
doConj (H a) = do
    a' <- doConj a
    return (H a')
doConj a@(DM _ _ _ _ _) = do
    liftV (\x -> V.doConj x >> return ()) a
    return a

-- | Scale every element in a matrix by the given value.
scaleBy :: (BLAS1 e) => e -> IOMatrix (m,n) e -> IO ()
scaleBy k = liftV (\x -> V.scaleBy k x >> return ())
    
-- | Scale every element by the given value.
shiftBy :: (BLAS1 e) => e -> IOMatrix (m,n) e -> IO ()
shiftBy k = liftV (V.shiftBy k)

-- | Divide every element by the given value.
invScaleBy :: (BLAS1 e) => e -> IOMatrix (m,n) e -> IO ()
invScaleBy k = liftV (V.invScaleBy k)

axpy :: (BLAS1 e) => e -> DMatrix t (m,n) e -> IOMatrix (m,n) e -> IO ()
axpy alpha = liftV2 (V.axpy alpha)

-- | In-place elementwise add.
(+=) :: (BLAS1 e) => IOMatrix (m,n) e -> DMatrix t (m,n) e -> IO ()
(+=) = liftV2 (V.+=)

-- | In-place elementwise subtract.
(-=) :: (BLAS1 e) => IOMatrix (m,n) e -> DMatrix t (m,n) e -> IO ()
(-=) = liftV2 (V.-=)

-- | In-place elementwise product.
(*=) :: (BLAS2 e) => IOMatrix (m,n) e -> DMatrix t (m,n) e -> IO ()
(*=) = liftV2 (V.*=)

-- | In-place elementwise divide.
(//=) :: (BLAS2 e) => IOMatrix (m,n) e -> DMatrix t (m,n) e -> IO ()
(//=) = liftV2 (V.//=)

blasTransOf :: DMatrix t (m,n) e -> CBLASTrans
blasTransOf a = 
    case (isHerm a) of
          False -> noTrans
          True  -> conjTrans

flipShape :: (Int,Int) -> (Int,Int)
flipShape (m,n) = (n,m)

-- | @gemv alpha a x beta y@ replaces @y := alpha a * x + beta y@
gemv :: (BLAS3 e) => e -> DMatrix s (m,n) e -> DVector t n e -> e -> IOVector m e -> IO ()
gemv alpha a x beta y
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


-- | @gemm alpha a b beta c@ replaces @c := alpha a * b + beta c@.
gemm :: (BLAS3 e) => e -> DMatrix s (m,k) e -> DMatrix t (k,n) e -> e -> IOMatrix (m,n) e -> IO ()
gemm alpha a b beta c
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


unaryOp :: (BLAS1 e) => (IOMatrix (m,n) e -> IO ()) 
    -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
unaryOp f a = do
    a' <- newCopy a
    f (unsafeThaw a')
    return (unsafeCoerce a')

binaryOp :: (BLAS1 e) => String -> (IOMatrix (m,n) e -> DMatrix t (m,n) e -> IO ()) 
    -> DMatrix s (m,n) e -> DMatrix t (m,n) e -> IO (DMatrix r (m,n) e)
binaryOp name f a b = 
    checkMatMatOp name (shape a) (shape b) >> do
        a' <- newCopy a
        f (unsafeThaw a') b
        return (unsafeCoerce a')
        
        
-- | Multiply a matrix by a vector.
apply :: (BLAS3 e) => Matrix (m,n) e -> Vector n e -> Vector m e
apply = sapply 1

-- | Multiply a scaled matrix by a vector.
sapply :: (BLAS3 e) => e -> Matrix (m,n) e -> Vector n e -> Vector m e
sapply alpha a x = unsafePerformIO $ getSApply alpha a x
{-# NOINLINE sapply #-}

-- | Multiply a scaled matrix by a matrix.
sapplyMat :: (BLAS3 e) => e -> Matrix (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
sapplyMat alpha a b = unsafePerformIO $ getSApplyMat alpha a b
{-# NOINLINE sapplyMat #-}
    
-- | Multiply a matrix by a matrix.
applyMat :: (BLAS3 e) => Matrix (m,k) e -> Matrix (k,n) e -> Matrix (m,n) e
applyMat = sapplyMat 1

-- | Create a new matrix by scaling another matrix by the given value.
scale :: (BLAS1 e) => e -> Matrix (m,n) e -> Matrix (m,n) e
scale k a = unsafePerformIO $ getScaled k a
{-# NOINLINE scale #-}

-- | Form a new matrix by adding a value to every element in a matrix.
shift :: (BLAS1 e) => e -> Matrix (m,n) e -> Matrix (m,n) e
shift k a = unsafePerformIO $ getShifted k a
{-# NOINLINE shift #-}

-- | Form a new matrix by dividing every element by a value.
invScale :: (BLAS1 e) => e -> Matrix (m,n) e -> Matrix (m,n) e
invScale k a = unsafePerformIO $ getInvScaled k a
{-# NOINLINE invScale #-}

-- | Create a new matrix by taking the elementwise sum of two matrices.
add :: (BLAS1 e) => e -> Matrix (m,n) e -> e -> Matrix (m,n) e -> Matrix (m,n) e
add alpha a beta b = unsafePerformIO $ getSum alpha a beta b
{-# NOINLINE add #-}

-- | Create a new matrix by taking the elementwise sum of two matrices.
plus :: (BLAS1 e) => Matrix (m,n) e -> Matrix (m,n) e -> Matrix (m,n) e
plus a b = add 1 a 1 b

-- | Create a new matrix by taking the elementwise difference of two matrices.
minus :: (BLAS1 e) => Matrix (m,n) e -> Matrix (m,n) e -> Matrix (m,n) e
minus a b = unsafePerformIO $ getDiff a b
{-# NOINLINE minus #-}

-- | Create a new matrix by taking the elementwise product of two matrices.
times :: (BLAS2 e) => Matrix (m,n) e -> Matrix (m,n) e -> Matrix (m,n) e
times a b = unsafePerformIO $ getProduct a b
{-# NOINLINE times #-}

-- | Create a new matrix by taking the elementwise ratio of two matrices.
divide :: (BLAS2 e) => Matrix (m,n) e -> Matrix (m,n) e -> Matrix (m,n) e
divide a b = unsafePerformIO $ getRatio a b
{-# NOINLINE divide #-}        

{-# RULES
"scale.plus/add"   forall k l x y. plus (scale k x) (scale l y) = add k x l y
"scale1.plus/add"  forall k x y.   plus (scale k x) y = add k x 1 y
"scale2.plus/add"  forall k x y.   plus x (scale k y) = add 1 x k y

"scale.minus/add"  forall k l x y. minus (scale k x) (scale l y) = add k x (-l) y
"scale1.minus/add" forall k x y.   minus (scale k x) y = add k x (-1) y
"scale2.minus/add" forall k x y.   minus x (scale k y) = add 1 x (-k) y

"scale.apply/sapply"       forall k a x. apply (scale k a) x = sapply k a x
"scale.applyMat/sapplyMat" forall k a b. applyMat (scale k a) b = sapplyMat k a b
  #-}
  