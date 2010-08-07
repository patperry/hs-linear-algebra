{-# LANGUAGE FlexibleInstances, BangPatterns #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.VMath.VNum
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector Num operations.
--

module Foreign.VMath.VNum
    where
     
import Foreign( Storable, Ptr, with, peek, advancePtr )
import Foreign.Storable.Complex()
import Data.Complex( Complex(..) )

import Foreign.BLAS( copy )
import Foreign.VMath.Double
import Foreign.VMath.Zomplex
        
-- | Types with vectorized 'Num' operations.
class (Storable a, Num a) => VNum a where
    vSum :: Int -> Ptr a -> IO a
    vShift :: Int -> a -> Ptr a -> Ptr a -> IO ()        
    vAdd :: Int -> Ptr a -> Ptr a -> Ptr a -> IO ()
    vSub :: Int -> Ptr a -> Ptr a -> Ptr a -> IO ()
    vAxpby :: Int -> a -> Ptr a -> a -> Ptr a -> Ptr a -> IO ()

    vScale :: Int -> a -> Ptr a -> Ptr a -> IO ()    
    vMul :: Int -> Ptr a -> Ptr a -> Ptr a -> IO ()

    vConj :: Int -> Ptr a -> Ptr a -> IO ()    
    vNeg :: Int -> Ptr a -> Ptr a -> IO ()
    vAbs :: Int -> Ptr a -> Ptr a -> IO ()
    vSgn :: Int -> Ptr a -> Ptr a -> IO ()

    vSum = vSumAcc 0
      where vSumAcc !acc n p | n <= 0    = return acc
                             | otherwise = do
                                 e <- peek p
                                 vSumAcc (acc+e) (n-1) (p `advancePtr` 1)

instance VNum Double where
    vShift = vdShift    
    {-# INLINE vShift #-}
    vAdd = vdAdd
    {-# INLINE vAdd #-}
    vSub = vdSub
    {-# INLINE vSub #-}
    vAxpby = vdAxpby
    {-# INLINE vAxpby #-}
    vScale = vdScale
    {-# INLINE vScale #-}
    vMul = vdMul
    {-# INLINE vMul #-}
    vConj n src dst = copy n src 1 dst 1
    {-# INLINE vConj #-}
    vNeg = vdNeg
    {-# INLINE vNeg #-}
    vAbs = vdAbs
    {-# INLINE vAbs #-}
    vSgn = vdSgn
    {-# INLINE vSgn #-}

instance VNum (Complex Double) where
    vShift n alpha pX pZ = with alpha $ \pAlpha ->
        vzShift n pAlpha pX pZ
    {-# INLINE vShift #-}
    vAdd = vzAdd
    {-# INLINE vAdd #-}
    vSub = vzSub
    {-# INLINE vSub #-}
    vAxpby n alpha pX beta pY pZ =
        with alpha $ \pAlpha ->
        with beta $ \pBeta ->
            vzAxpby n pAlpha pX pBeta pY pZ
    {-# INLINE vAxpby #-}
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
