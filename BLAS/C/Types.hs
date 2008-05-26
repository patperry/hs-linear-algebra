-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.C.Types
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.C.Types (
    CBLASOrder,
    CBLASTrans,
    CBLASUpLo,
    CBLASDiag,
    CBLASSide,

    cblasOrder,
    cblasTrans,
    cblasUpLo,
    cblasDiag,
    cblasSide,
    ) where
        
import BLAS.Types
        
newtype CBLASOrder = CBLASOrder Int deriving (Eq, Show)
newtype CBLASTrans = CBLASTrans Int deriving (Eq, Show)
newtype CBLASUpLo  = CBLASUpLo  Int deriving (Eq, Show)
newtype CBLASDiag  = CBLASDiag  Int deriving (Eq, Show)
newtype CBLASSide  = CBLASSide  Int deriving (Eq, Show)

cblasOrder :: Order -> CBLASOrder
cblasOrder RowMajor = CBLASOrder 101
cblasOrder ColMajor = CBLASOrder 102

cblasTrans :: Trans -> CBLASTrans
cblasTrans NoTrans   = CBLASTrans 111
--cblasTrans Trans     = CBLASTrans 112
cblasTrans ConjTrans = CBLASTrans 113

cblasUpLo :: UpLo -> CBLASUpLo
cblasUpLo Upper = CBLASUpLo 121
cblasUpLo Lower = CBLASUpLo 122

cblasDiag :: Diag -> CBLASDiag
cblasDiag NonUnit = CBLASDiag 131
cblasDiag Unit    = CBLASDiag 132

cblasSide :: Side -> CBLASSide
cblasSide LeftSide  = CBLASSide 141
cblasSide RightSide = CBLASSide 142
