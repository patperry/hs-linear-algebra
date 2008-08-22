-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.UnsafeInterleaveM
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.UnsafeInterleaveM (
    UnsafeInterleaveM(..),
    ) where

import Control.Monad.ST( ST, unsafeInterleaveST )
import System.IO.Unsafe( unsafeInterleaveIO )

class (Monad m) => UnsafeInterleaveM m where
    unsafeInterleaveM :: m a -> m a

instance UnsafeInterleaveM IO where
    unsafeInterleaveM = unsafeInterleaveIO
instance UnsafeInterleaveM (ST s) where
    unsafeInterleaveM = unsafeInterleaveST
