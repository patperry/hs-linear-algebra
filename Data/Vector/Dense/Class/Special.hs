-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Special
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Special (
    -- * Special vectors
    newZeroVector,
    newConstantVector,
    newBasisVector,
    setZeroVector,
    setConstantVector,
    setBasisVector,
        
    ) where

import Foreign

import BLAS.Elem
import BLAS.Tensor

import Data.Vector.Dense.Class.Internal
import Data.Vector.Dense.Class.Creating


-- | @newBasisVector n i@ creates a vector of length @n@ that is all zero 
-- except for at position @i@, where it equal to one.
newBasisVector :: (WriteVector x m, Elem e) => Int -> Int -> m (x n e)
newBasisVector n i = do
    x <- newVector_ n
    setBasisVector i x
    return x


-- | @setBasis x i@ sets the @i@th coordinate of @x@ to @1@, and all other
-- coordinates to @0@.  If the vector has been scaled, it is possible that
-- @readVector x i@ will not return exactly @1@.  See 'setElem'.
setBasisVector :: (WriteVector x m, Elem e) => Int -> x n e -> m ()
setBasisVector i x
    | i < 0 || i >= dim x =
        fail $ "tried to set a vector of dimension `" ++ show (dim x) ++ "'"
               ++ " to basis vector `" ++ show i ++ "'"
    | otherwise = do
        setZeroVector x
        unsafeWriteElem x i 1 
