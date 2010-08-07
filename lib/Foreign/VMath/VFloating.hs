{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.VMath.VFloating
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector Floating operations.
--

module Foreign.VMath.VFloating (
    VFloating(..)
    ) where
     
import Foreign( Ptr, Storable, peek, poke, advancePtr )
import Data.Complex( Complex(..) )

import Foreign.VMath.VFractional
import Foreign.VMath.Double  
import Foreign.VMath.Zomplex

-- | Types with vectorized 'Floating' operations.
class (VFractional a, Floating a) => VFloating a where
    vExp :: Int -> Ptr a -> Ptr a -> IO ()
    vSqrt :: Int -> Ptr a -> Ptr a -> IO ()
    vLog :: Int -> Ptr a -> Ptr a -> IO ()        
    vPow :: Int -> Ptr a -> Ptr a -> Ptr a -> IO ()            
    vSin :: Int -> Ptr a -> Ptr a -> IO ()                
    vCos :: Int -> Ptr a -> Ptr a -> IO ()                    
    vTan :: Int -> Ptr a -> Ptr a -> IO ()                        
    vASin :: Int -> Ptr a -> Ptr a -> IO ()                
    vACos :: Int -> Ptr a -> Ptr a -> IO ()                    
    vATan :: Int -> Ptr a -> Ptr a -> IO ()                        
    vSinh :: Int -> Ptr a -> Ptr a -> IO ()                
    vCosh :: Int -> Ptr a -> Ptr a -> IO ()                    
    vTanh :: Int -> Ptr a -> Ptr a -> IO ()                        
    vASinh :: Int -> Ptr a -> Ptr a -> IO ()                
    vACosh :: Int -> Ptr a -> Ptr a -> IO ()                    
    vATanh :: Int -> Ptr a -> Ptr a -> IO ()

    vExp = vop exp
    vSqrt = vop sqrt
    vLog = vop log
    vPow = vop2 (**)
    vSin = vop sin
    vCos = vop cos
    vTan = vop tan
    vASin = vop asin
    vACos = vop acos
    vATan = vop atan
    vSinh = vop sinh
    vCosh = vop cosh
    vTanh = vop tanh
    vASinh = vop asinh
    vACosh = vop acosh
    vATanh = vop atanh
    
    
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
    
instance VFloating Double where
    vExp = vdExp
    {-# INLINE vExp #-}
    vSqrt = vdSqrt
    {-# INLINE vSqrt #-}
    vLog = vdLog
    {-# INLINE vLog #-}
    vPow = vdPow
    {-# INLINE vPow #-}
    vSin = vdSin
    {-# INLINE vSin #-}
    vCos = vdCos
    {-# INLINE vCos #-}
    vTan = vdTan
    {-# INLINE vTan #-}
    vASin = vdASin
    {-# INLINE vASin #-}
    vACos = vdACos
    {-# INLINE vACos #-}
    vATan = vdATan
    {-# INLINE vATan #-}
    vSinh = vdSinh
    {-# INLINE vSinh #-}
    vCosh = vdCosh
    {-# INLINE vCosh #-}
    vTanh = vdTanh
    {-# INLINE vTanh #-}
    vASinh = vdASinh
    {-# INLINE vASinh #-}
    vACosh = vdACosh
    {-# INLINE vACosh #-}
    vATanh = vdATanh
    {-# INLINE vATanh #-}


instance VFloating (Complex Double) where
    vSqrt = vzSqrt
    {-# INLINE vSqrt #-}
    vLog = vzLog
    {-# INLINE vLog #-}
    
    {- These functions have branch cuts in the wrong places 
    vExp = vzExp
    {-# INLINE vExp #-}
    vPow = vzPow
    {-# INLINE vPow #-}
    vSin = vzSin
    {-# INLINE vSin #-}
    vCos = vzCos
    {-# INLINE vCos #-}
    vTan = vzTan
    {-# INLINE vTan #-}
    vASin = vzASin
    {-# INLINE vASin #-}
    vACos = vzACos
    {-# INLINE vACos #-}
    vATan = vzATan
    {-# INLINE vATan #-}
    vSinh = vzSinh
    {-# INLINE vSinh #-}
    vCosh = vzCosh
    {-# INLINE vCosh #-}
    vTanh = vzTanh
    {-# INLINE vTanh #-}
    vASinh = vzASinh
    {-# INLINE vASinh #-}
    vACosh = vzACosh
    {-# INLINE vACosh #-}
    vATanh = vzATanh
    {-# INLINE vATanh #-}
    -}
