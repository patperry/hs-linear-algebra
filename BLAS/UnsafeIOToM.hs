-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.UnsafeIOToM
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.UnsafeIOToM (
    UnsafeIOToM(..),
    ) where

import Control.Monad.ST ( ST, unsafeIOToST, unsafeInterleaveST )
import System.IO.Unsafe( unsafeInterleaveIO )


class (Monad m) => UnsafeIOToM m where
    unsafeIOToM :: IO a -> m a
    unsafeInterleaveM :: m a -> m a
    
instance UnsafeIOToM IO where
    unsafeIOToM = id
    {-# INLINE unsafeIOToM #-}
    
    unsafeInterleaveM = unsafeInterleaveIO
    {-# INLINE unsafeInterleaveM #-}

instance UnsafeIOToM (ST s) where
    unsafeIOToM = unsafeIOToST
    {-# INLINE unsafeIOToM #-}

    unsafeInterleaveM = unsafeInterleaveST
    {-# INLINE unsafeInterleaveM #-}
