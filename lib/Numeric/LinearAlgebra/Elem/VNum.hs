{-# LANGUAGE FlexibleInstances, BangPatterns #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Elem.VNum
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector operations.
--

module Numeric.LinearAlgebra.Elem.VNum (
    VNum(..)
    ) where
     
import Foreign( Storable, Ptr, with, peek, poke, advancePtr )
import Foreign.Storable.Complex()
import Data.Complex( Complex(..) )

import Numeric.LinearAlgebra.Elem.Double  
import Numeric.LinearAlgebra.Elem.Zomplex
        
-- | Types with vectorized 'Num' operations.
class (Storable a, Num a) => VNum a where
    vShift :: Int -> a -> Ptr a -> Ptr a -> IO ()
    vShift n k = vop (k+) n
    
    vAdd :: Int -> Ptr a -> Ptr a -> Ptr a -> IO ()
    vAdd = vop2 (+)
    
    vSub :: Int -> Ptr a -> Ptr a -> Ptr a -> IO ()
    vSub = vop2 (-)
    
    vAxpby :: Int -> a -> Ptr a -> a -> Ptr a -> Ptr a -> IO ()
    vAxpby n a px b = vop2 (\e f -> a * e + b * f) n px

    vScale :: Int -> a -> Ptr a -> Ptr a -> IO ()
    vScale n k = vop (k*) n
    
    vMul :: Int -> Ptr a -> Ptr a -> Ptr a -> IO ()
    vMul = vop2 (*)

    vConj :: Int -> Ptr a -> Ptr a -> IO ()

    vNeg :: Int -> Ptr a -> Ptr a -> IO ()
    vNeg = vop negate
    
    vAbs :: Int -> Ptr a -> Ptr a -> IO ()
    vAbs = vop abs
    
    vSgn :: Int -> Ptr a -> Ptr a -> IO ()
    vSgn = vop signum

    vSum :: Int -> Ptr a -> IO a
    vSum = vSumAcc 0
      where vSumAcc !acc n p | n <= 0    = return acc
                             | otherwise = do
                                 e <- peek p
                                 vSumAcc (acc+e) (n-1) (p `advancePtr` 1)


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


instance VNum Double where
    vShift = vdShift    
    {-# INLINE vShift #-}
    vAdd = vdAdd
    {-# INLINE vAdd #-}
    vSub = vdSub
    {-# INLINE vSub #-}
    vScale = vdScale
    {-# INLINE vScale #-}
    vMul = vdMul
    {-# INLINE vMul #-}
    vConj n src dst = dcopy n src 1 dst 1
    {-# INLINE vConj #-}
    vNeg = vdNeg
    {-# INLINE vNeg #-}
    vAbs = vdAbs
    {-# INLINE vAbs #-}
    vSgn = vdSgn
    {-# INLINE vSgn #-}
    {- This funciton disagrees with Haskell
    vAxpby = vdAxpby
    {-# INLINE vAxpby #-}
    -}

instance VNum (Complex Double) where
    vShift n alpha pX pZ = with alpha $ \pAlpha ->
        vzShift n pAlpha pX pZ
    {-# INLINE vShift #-}
    vAdd = vzAdd
    {-# INLINE vAdd #-}
    vSub = vzSub
    {-# INLINE vSub #-}
    {-
    vAxpby n alpha pX beta pY pZ =
        with alpha $ \pAlpha ->
        with beta $ \pBeta ->
            vzAxpby n pAlpha pX pBeta pY pZ
    {-# INLINE vAxpby #-}
    -}
    vScale n alpha pX pZ = with alpha $ \pAlpha ->
        vzScale n pAlpha pX pZ
    {-# INLINE vScale #-}
    vMul = vzMul
    {-# INLINE vMul #-}
    vConj = vzConj
    {-# INLINE vConj #-}
    vNeg = vzNeg
    {-# INLINE vNeg #-}
    vAbs = vzAbs
    {-# INLINE vAbs #-}
    vSgn = vzSgn
    {-# INLINE vSgn #-}
