-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.IO
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.IO (
    module Data.Vector.Dense.Internal,
    
    ) where

import Data.Vector.Dense.Internal hiding ( DV, storageIOVector, offsetIOVector, 
    dimIOVector, strideIOVector, isConjIOVector, 
    unsafeAxpyVector, unsafeMulVector, unsafeDivVector, 
    unsafeDoVectorOp2 )
import BLAS.Numeric
