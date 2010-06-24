{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Elem.VFractional
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector operations.
--

module BLAS.Elem.VFractional
    where
     
import Foreign( Ptr )
import Data.Complex( Complex(..) )

import BLAS.Elem.VNum
import BLAS.Elem.Double  
import BLAS.Elem.Zomplex
        
-- | Types with vectorized 'Fractional' operations.
class (VNum a, Fractional a) => VFractional a where
    vInv :: Int -> Ptr a -> Ptr a -> IO ()    
    vDiv :: Int -> Ptr a -> Ptr a -> Ptr a -> IO ()

instance VFractional Double where
    vInv = vdInv
    {-# INLINE vInv #-}
    vDiv = vdDiv
    {-# INLINE vDiv #-}

instance VFractional (Complex Double) where
    vInv = vzInv
    {-# INLINE vInv #-}
    vDiv = vzDiv
    {-# INLINE vDiv #-}
