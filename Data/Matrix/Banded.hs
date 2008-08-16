{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded (
    -- * Banded matrix type
    Banded,
    
    module BLAS.Matrix.Base,
    module BLAS.Tensor.Base,
    module BLAS.Tensor.Dense.Immutable,
    module BLAS.Tensor.Immutable,
    module BLAS.Tensor.Scalable,
    module Data.Matrix.Banded.Operations,

    -- * Creating banded matrices
    banded, 
    listsBanded,
    
    -- * Properties
    numLower,
    numUpper,
    bandwidth,
    
    -- * Rows and columns
    row,
    col,
    rows,
    cols,

    -- * Diagonals
    diag,
    toLists,

    -- * Casting matrices
    coerceBanded,
    
    -- * Converting between vectors and matrices
    
    -- * Unsafe operations
    unsafeBanded,
    unsafeRow,
    unsafeCol,
    unsafeDiag,
    
    ) where

import BLAS.Access
import BLAS.Elem ( BLAS1 )
import BLAS.Matrix.Base hiding ( Matrix )
import BLAS.Tensor.Base
import BLAS.Tensor.Dense.Immutable
import BLAS.Tensor.Immutable
import BLAS.Tensor.Scalable

import Data.Matrix.Banded.Internal
import qualified Data.Matrix.Banded.Internal as B
import Data.Matrix.Banded.Operations hiding ( RMatrix(..), getScaled, 
    getInvScaled, doConj, scaleBy, invScaleBy )
import Data.Vector.Dense hiding ( scale, invScale )

instance (BLAS1 e) => Scalable (BMatrix Imm (m,n)) e where
    (*>) = scale
