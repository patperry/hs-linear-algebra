{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Elem.VFractional
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector operations.
--

module Numeric.LinearAlgebra.Elem.VFractional (
    VFractional(..)
    ) where
     
import Foreign( Ptr, Storable, peek, poke, advancePtr )
import Data.Complex( Complex(..) )

import Numeric.LinearAlgebra.Elem.VNum
import Numeric.LinearAlgebra.Elem.Double  
-- import Numeric.LinearAlgebra.Elem.Zomplex
        
-- | Types with vectorized 'Fractional' operations.
class (VNum a, Fractional a) => VFractional a where
    vInv :: Int -> Ptr a -> Ptr a -> IO ()
    vInv = vop recip
   
    vDiv :: Int -> Ptr a -> Ptr a -> Ptr a -> IO ()
    vDiv = vop2 (/)

vop :: (Storable a, Storable b)
    => (a -> b) -> Int -> Ptr a -> Ptr b -> IO ()
vop f n src dst | n <= 0    = return ()
                | otherwise = do
    a <- peek src
    poke dst $ f a
    vop f (n-1) (src `advancePtr` 1) (dst `advancePtr` 1)
{-# INLINE vop #-}

vop2 :: (Storable a1, Storable a2, Storable b)
     => (a1 -> a2 -> b) -> Int -> Ptr a1 -> Ptr a2 -> Ptr b -> IO ()
vop2 f n src1 src2 dst | n <= 0    = return ()
                       | otherwise = do
    a1 <- peek src1
    a2 <- peek src2
    poke dst $ f a1 a2
    vop2 f (n-1) (src1 `advancePtr` 1) (src2 `advancePtr` 1)
         (dst `advancePtr` 1)
{-# INLINE vop2 #-}


instance VFractional Double where
    vInv = vdInv
    {-# INLINE vInv #-}
    vDiv = vdDiv
    {-# INLINE vDiv #-}

instance VFractional (Complex Double) where
    {- These functions disagree with Haskell
    vInv = vzInv
    {-# INLINE vInv #-}
    vDiv = vzDiv
    {-# INLINE vDiv #-}
    -}
