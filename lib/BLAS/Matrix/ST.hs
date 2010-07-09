-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.ST
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Mutable matrices in the ST monad.

module BLAS.Matrix.ST (
    -- * The @STMatrix@ data type
    STMatrix,
    -- runMatrix,
    
    -- * Read-only Dense matrix type class
    RMatrix( dimMatrix, unsafeWithMatrix ),
    
    -- * Creating new matrices
    newMatrix_,
    newMatrix,
    
    -- * Copying matrices
    newCopyMatrix,
    copyToMatrix,

    -- * Matrix views
    spliceMatrix,
    splitRowsMatrixAt,
    splitColsMatrixAt,    

    -- * Reading and writing matrix elements
    readMatrix,
    writeMatrix,
    indicesMatrix,
    getElemsMatrix,
    getElemsMatrix',
    getAssocsMatrix,
    getAssocsMatrix',
    setElemsMatrix,
    setAssocsMatrix,    

    -- * Conversions between mutable and immutable matrices
    {- freezeMatrix,
    unsafeFreezeMatrix,
    thawMatrix,
    unsafeThawMatrix, -}
    
    -- * Matrix views of arrays
    matrixViewArray,

    ) where

-- import BLAS.Matrix.Base
import BLAS.Matrix.STBase
