{-# LANGUAGE Rank2Types, FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Vector.Dense.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.ST (
    -- * The @STVector@ data type
    STVector,
    runSTVector,

    module Data.Vector.Dense.IO,
    ) where

import BLAS.Elem ( Elem, BLAS1 )

import Control.Monad.ST

import Data.Vector.Dense ( Vector, UnsafeFreezeVector(..), UnsafeThawVector(..) )
import Data.Vector.Dense.Class.Base
import Data.Vector.Dense.IO hiding ( IOVector )
import Data.Vector.Dense.Internal ( newCopyVector, unsafeCopyVector, 
    unsafeAxpyVector, unsafeMulVector, unsafeDivVector, unsafeDoVectorOp2 )
import qualified Data.Vector.Dense.IO as IO

newtype STVector s n e = ST (IO.IOVector n e)

liftSTVector :: (IO.IOVector n e -> a) -> STVector s n e -> a
liftSTVector f (ST x) = f x

unsafeLiftSTVector :: (IO.IOVector n e -> IO a) -> STVector s n e -> ST s a
unsafeLiftSTVector f = unsafeIOToST . liftSTVector f

instance UnsafeFreezeVector (STVector s) where
    unsafeFreezeVector = liftSTVector unsafeFreezeVector
instance UnsafeThawVector (STVector s) where
    unsafeThawVector = ST . unsafeThawVector

runSTVector :: (forall s . ST s (STVector s n e)) -> Vector n e
runSTVector x = runST $ x >>= return . unsafeFreezeVector

instance BaseTensor (STVector s) Int e where
    shape  = liftSTVector shape
    bounds = liftSTVector bounds

instance (Elem e) => ReadTensor (STVector s) Int e (ST s) where
    getSize            = unsafeLiftSTVector getSize
    
    getAssocs          = unsafeLiftSTVector getAssocs
    getIndices         = unsafeLiftSTVector getIndices
    getElems           = unsafeLiftSTVector getElems
    
    getAssocs'         = unsafeLiftSTVector getAssocs'
    getIndices'        = unsafeLiftSTVector getIndices'
    getElems'          = unsafeLiftSTVector getElems'

    unsafeReadElem x i = unsafeLiftSTVector (flip unsafeReadElem i) x

instance (Elem e) => WriteTensor (STVector s) Int e (ST s) where
    newZero n                  = unsafeIOToST $ newZero n >>= return . ST
    newConstant e n            = unsafeIOToST $ newConstant e n >>= return . ST
    unsafeSwap (ST x) (ST y)   = unsafeIOToST $ unsafeSwap x y
    setConstant k              = unsafeLiftSTVector (setConstant k)
    setZero                    = unsafeLiftSTVector setZero
    modifyWith f               = unsafeLiftSTVector (modifyWith f)
    unsafeWriteElem (ST x) i e = unsafeIOToST $ unsafeWriteElem x i e
    canModifyElem x i          = unsafeLiftSTVector (flip canModifyElem i) x
    
instance (Elem e) => BaseVector (STVector s) e where
    stride                    = liftSTVector stride
    isConj                    = liftSTVector isConj
    conjVector (ST x)         = ST $ conjVector x
    unsafeSubvectorWithStride s (ST x) o n = ST $ unsafeSubvectorWithStride s x o n
    vectorViewArray f o n s c = ST $ vectorViewArray f o n s c
    withVectorPtr             = liftSTVector withVectorPtr

instance (Elem e) => ReadVector (STVector s) e (ST s) where
    
instance (Elem e) => WriteVector (STVector s) e (ST s) where
    newVector_ n = unsafeIOToST $ newVector_ n >>= return . ST

instance (BLAS1 e) => ReadNumeric (STVector s) Int e (ST s) where
    
instance (BLAS1 e) => WriteNumeric (STVector s) Int e (ST s) where
    doConj                     = unsafeLiftSTVector doConj
    scaleBy k                  = unsafeLiftSTVector (scaleBy k)
    shiftBy k                  = unsafeLiftSTVector (shiftBy k)

instance (BLAS1 e, ReadVector x e (ST s)) => CopyTensor x (STVector s) Int e (ST s) where
    newCopyTensor    = newCopyVector    
    unsafeCopyTensor = unsafeCopyVector
    
instance (BLAS1 e, ReadVector x e (ST s)) => Numeric2 x (STVector s) Int e (ST s) where
    unsafeAxpy = unsafeAxpyVector
    unsafeMul  = unsafeMulVector
    unsafeDiv  = unsafeDivVector

instance (BLAS1 e, ReadVector x e (ST s), ReadVector y e (ST s)) => Numeric3 x y (STVector s) Int e (ST s) where
    unsafeDoAdd = unsafeDoVectorOp2 $ flip $ unsafeAxpyVector 1
    unsafeDoSub = unsafeDoVectorOp2 $ flip $ unsafeAxpyVector (-1)
    unsafeDoMul = unsafeDoVectorOp2 $ unsafeMulVector
    unsafeDoDiv = unsafeDoVectorOp2 $ unsafeDivVector
