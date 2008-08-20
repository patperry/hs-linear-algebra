{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses, UndecidableInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.ST (
    -- * The @STMatrix@ data type
    STMatrix,
    -- runSTMatrix,

    module Data.Matrix.Dense.IO,

    unsafeIOMatrixToSTMatrix,
    unsafeSTMatrixToIOMatrix,
    ) where

import BLAS.Elem ( Elem, BLAS1 )

import Control.Monad.ST

-- import Data.Vector.Dense ( Vector, UnsafeFreezeVector(..), UnsafeThawVector(..) )
import Data.Vector.Dense.ST
import Data.Matrix.Dense.IO hiding ( IOMatrix )
import qualified Data.Matrix.Dense.IO as IO
import qualified BLAS.Matrix.Base as BLAS

newtype STMatrix s mn e = ST (IO.IOMatrix mn e)

unsafeIOMatrixToSTMatrix :: IO.IOMatrix mn e -> STMatrix s mn e
unsafeIOMatrixToSTMatrix = ST

unsafeSTMatrixToIOMatrix :: STMatrix s mn e -> IO.IOMatrix mn e
unsafeSTMatrixToIOMatrix (ST a) = a


liftSTMatrix :: (IO.IOMatrix mn e -> a) -> STMatrix s mn e -> a
liftSTMatrix f (ST a) = f a

unsafeLiftSTMatrix :: (IO.IOMatrix n e -> IO a) -> STMatrix s n e -> ST s a
unsafeLiftSTMatrix f = unsafeIOToST . liftSTMatrix f

instance BaseTensor (STMatrix s) (Int,Int) e where
    shape  = liftSTMatrix shape
    bounds = liftSTMatrix bounds

instance (Elem e) => ReadTensor (STMatrix s) (Int,Int) e (ST s) where
    getSize            = unsafeLiftSTMatrix getSize
    getAssocs          = unsafeLiftSTMatrix getAssocs
    getIndices         = unsafeLiftSTMatrix getIndices
    getElems           = unsafeLiftSTMatrix getElems
    getAssocs'         = unsafeLiftSTMatrix getAssocs'
    getIndices'        = unsafeLiftSTMatrix getIndices'
    getElems'          = unsafeLiftSTMatrix getElems'
    unsafeReadElem x i = unsafeLiftSTMatrix (flip unsafeReadElem i) x    

instance (Elem e) => WriteTensor (STMatrix s) (Int,Int) e (ST s) where
    newZero mn                 = unsafeIOToST $ newZero mn >>= return . ST
    newConstant e mn           = unsafeIOToST $ newConstant e mn >>= return . ST
    unsafeSwap (ST a) (ST b)   = unsafeIOToST $ unsafeSwap a b
    setConstant k              = unsafeLiftSTMatrix (setConstant k)
    setZero                    = unsafeLiftSTMatrix setZero
    modifyWith f               = unsafeLiftSTMatrix (modifyWith f)
    unsafeWriteElem (ST a) i e = unsafeIOToST $ unsafeWriteElem a i e
    canModifyElem a i          = unsafeLiftSTMatrix (flip canModifyElem i) a

instance (Elem e) => BLAS.BaseMatrix (STMatrix s) e where
    herm (ST a) = ST (herm a)
    
instance (Elem e) => BaseMatrix (STMatrix s) e where
    lda                          = liftSTMatrix lda
    isHerm                       = liftSTMatrix isHerm
    unsafeSubmatrix (ST a) ij mn = ST $ unsafeSubmatrix a ij mn
    withMatrixPtr                = liftSTMatrix withMatrixPtr
    matrixViewArray f o mn l h   = ST $ matrixViewArray f o mn l h
    arrayFromMatrix (ST a)       = arrayFromMatrix a

instance (Elem e) => ReadMatrix (STMatrix s) (STVector s) e (ST s) where
    
instance (Elem e) => WriteMatrix (STMatrix s) (STVector s) e (ST s) where
    newMatrix_ mn = unsafeIOToST $ newMatrix_ mn >>= return . ST

instance (Elem e) => RowColView (STMatrix s) (STVector s) e where
    unsafeRowView (ST a) i = unsafeIOVectorToSTVector $ unsafeRowView a i
    unsafeColView (ST a) j = unsafeIOVectorToSTVector $ unsafeColView a j

instance (BLAS1 e) => RowColRead (STMatrix s) (STVector s) e (ST s) where
    unsafeGetRow a i = newCopyVector (unsafeRowView a i)
    unsafeGetCol a j = newCopyVector (unsafeColView a j)

instance (Elem e) => DiagView (STMatrix s) (STVector s) e where
    unsafeDiagView (ST a) i = unsafeIOVectorToSTVector $ unsafeDiagView a i

instance (BLAS1 e) => DiagRead (STMatrix s) (STVector s) e (ST s) where
    unsafeGetDiag a i = newCopyVector (unsafeDiagView a i)

instance (BLAS1 e) => ReadNumeric (STMatrix s) (Int,Int) e (ST s) where

instance (BLAS1 e) => WriteNumeric (STMatrix s) (Int,Int) e (ST s) where
    doConj    = liftMatrix doConj
    scaleBy k = liftMatrix (scaleBy k)
    shiftBy k = liftMatrix (shiftBy k)

instance (BLAS1 e, ReadMatrix a x e (ST s)) => CopyTensor a (STMatrix s) (Int,Int) e (ST s) where
    newCopyTensor    = newCopyMatrix
    unsafeCopyTensor = unsafeCopyMatrix

instance (BLAS1 e, ReadMatrix a x e (ST s)) => Numeric2 a (STMatrix s) (Int,Int) e (ST s) where
    unsafeAxpy k = liftMatrix2 (unsafeAxpy k)
    unsafeMul    = liftMatrix2 unsafeMul
    unsafeDiv    = liftMatrix2 unsafeDiv

instance (BLAS1 e, ReadMatrix a x e (ST s), ReadMatrix b y e (ST s)) => Numeric3 a b (STMatrix s) (Int,Int) e (ST s) where
    unsafeDoAdd = unsafeDoMatrixOp2 $ flip $ unsafeAxpy 1
    unsafeDoSub = unsafeDoMatrixOp2 $ flip $ unsafeAxpy (-1)
    unsafeDoMul = unsafeDoMatrixOp2 $ unsafeMul
    unsafeDoDiv = unsafeDoMatrixOp2 $ unsafeDiv
