-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.LAPACK.Types
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Foreign.LAPACK.Types (
    LAEigJob,
    EigJob(..),
    withEigJob,
    
    LAEigRange,
    EigRange(..),
    withEigRange,
    
    ) where

import Foreign
import Foreign.C.String

import Foreign.BLAS.Types( LAInt )

newtype LAEigJob = LAEigJob CString

-- | A type of eigenvalue computation.
data EigJob = NoEigVectors -- ^ no eigenvectors (eigenvalues only)
            | EigVectors   -- ^ eigenvalues and eigenvectors
    deriving (Eq, Show)

withEigJob :: EigJob -> (LAEigJob -> IO a) -> IO a
withEigJob evjob f = flip withCString (f . LAEigJob) $ case evjob of
    NoEigVectors -> "N"
    EigVectors   -> "V"

newtype LAEigRange = LAEigRange CString

-- | A range of eigenvalues to compute.
data EigRange = AllEigs                       -- ^ All eigenvalues are computed
              | EigValueRange !Double !Double -- ^ Eigenvalues in the half-open
                                              --   interval @(vl,vu]@ are
                                              --   computed
              | EigIndexRange !Int    !Int    -- ^ Eigenvalues with indices in
                                              --   the closed interval
                                              --   @[il,iu]@ are computed.
                                              --   Index @0@ corresponds to
                                              --   the smallest eigenvalue.
    deriving (Eq, Show)
             
withEigRange :: EigRange
             -> (LAEigRange -> Ptr Double -> Ptr Double -> Ptr LAInt -> Ptr LAInt -> IO a)
             -> IO a
withEigRange eigrange f = case eigrange of
    AllEigs ->
        withCString "A" $ \pA ->
        alloca $ \pvl ->
        alloca $ \pvu ->
        alloca $ \pil ->
        alloca $ \piu ->
            f (LAEigRange pA) pvl pvu pil piu
    EigValueRange vl vu ->
        withCString "V" $ \pV ->
        with vl $ \pvl ->
        with vu $ \pvu ->
        alloca $ \pil ->
        alloca $ \piu ->
            f (LAEigRange pV) pvl pvu pil piu
    EigIndexRange il iu ->
        withCString "I" $ \pI ->
        alloca $ \pvl ->
        alloca $ \pvu ->
        with (toEnum $ il + 1) $ \pil ->
        with (toEnum $ iu + 1) $ \piu ->
            f (LAEigRange pI) pvl pvu pil piu
