{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor.Dense.Read
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Numeric.Write (
    WriteNumeric(..),

    -- * Unary operations that return a new result.
    getConj,
    getScaled,
    getShifted,
    
    module BLAS.Numeric.Read,
    ) where

import BLAS.Elem
import BLAS.Tensor
import BLAS.Numeric.Read

class (ReadNumeric x i e m, WriteTensor x i e m) => WriteNumeric x i e m where
    -- | Replace every element with its complex conjugate.
    doConj :: x n e -> m ()
    doConj = modifyWith conj

    -- | Scale every element in the vector by the given value.
    scaleBy :: e -> x n e -> m ()
    scaleBy 1 = const $ return ()
    scaleBy k = modifyWith (k*)

    -- | Add a value to every element in a vector.
    shiftBy :: e -> x n e -> m ()
    shiftBy 0 = const $ return ()
    shiftBy k = modifyWith (k+)


-- | Create a new object by conjugating another object.
getConj :: (CopyTensor x y i e m, WriteNumeric y i e m) =>
    x n e -> m (y n e)
getConj = getUnaryOp doConj
{-# INLINE getConj #-}

-- | Create a new object by scaling another object.
getScaled :: (CopyTensor x y i e m, WriteNumeric y i e m) =>
    e -> x n e -> m (y n e)
getScaled k = getUnaryOp (scaleBy k)
{-# INLINE getScaled #-}

-- | Create a new object by shifting another object.
getShifted :: (CopyTensor x y i e m, WriteNumeric y i e m) =>
    e -> x n e -> m (y n e)
getShifted k = getUnaryOp (shiftBy k)
{-# INLINE getShifted #-}


getUnaryOp :: (CopyTensor x y i e m) =>
    (y n e -> m ()) -> x n e -> m (y n e)
getUnaryOp f x = do
    y <- newCopy x
    f y
    return y
