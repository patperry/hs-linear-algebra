{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Elem.VFloating
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector operations.
--

module BLAS.Elem.VFloating
    where
     
import Foreign( Ptr )
import Data.Complex( Complex(..) )

import BLAS.Elem.VFractional
import BLAS.Elem.Double  
import BLAS.Elem.Zomplex

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
    vExp = vzExp
    {-# INLINE vExp #-}
    vSqrt = vzSqrt
    {-# INLINE vSqrt #-}
    vLog = vzLog
    {-# INLINE vLog #-}
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
