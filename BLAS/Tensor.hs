-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Tensor
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Tensor (
    module BLAS.Tensor.Base,
    module BLAS.Tensor.Immutable,
    module BLAS.Tensor.Read,
    module BLAS.Tensor.Copy,
    module BLAS.Tensor.Write,
    module BLAS.Tensor.Swap,
    ) where

import BLAS.Tensor.Base
import BLAS.Tensor.Immutable 
import BLAS.Tensor.Read
import BLAS.Tensor.Copy
import BLAS.Tensor.Swap
import BLAS.Tensor.Write
