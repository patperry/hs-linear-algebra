{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Read
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Read (
    ReadVector,
    
    unsafeIOToM,
    module Data.Vector.Dense.Class.Base,
    
    ) where

import BLAS.Internal ( UnsafeIOToM(..) )
import BLAS.Tensor
import Data.Vector.Dense.Class.Base

class (ReadTensor x Int e m, BaseVector x e, UnsafeIOToM m) => ReadVector x e m where
