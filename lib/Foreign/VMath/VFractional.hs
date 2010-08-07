{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.VMath.VFractional
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Vector Fractional operations.
--

module Foreign.VMath.VFractional
    where
     
import Foreign( Ptr )
import Data.Complex( Complex(..) )

import Foreign.VMath.VNum
import Foreign.VMath.Double  
import Foreign.VMath.Zomplex
        
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
