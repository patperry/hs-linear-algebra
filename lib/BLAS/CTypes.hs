{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.CTypes
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module BLAS.CTypes (
    CBLASOrder,
    CBLASTrans,
    CBLASUpLo,
    CBLASDiag,
    CBLASSide,

    rowMajor,
    colMajor,
    
    noTrans,
    trans,
    conjTrans,
    
    upper,
    lower,
    
    nonUnit,
    unit,
    
    leftSide,
    rightSide,

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

rowMajor, colMajor :: CBLASOrder
rowMajor = CBLASOrder 101
colMajor = CBLASOrder 102

cblasOrder :: OrderEnum -> CBLASOrder
cblasOrder RowMajor = CBLASOrder 101
cblasOrder ColMajor = CBLASOrder 102

noTrans, trans, conjTrans :: CBLASTrans
noTrans   = CBLASTrans 111
trans     = CBLASTrans 112
conjTrans = CBLASTrans 113

cblasTrans :: TransEnum -> CBLASTrans
cblasTrans NoTrans   = CBLASTrans 111
cblasTrans ConjTrans = CBLASTrans 113

upper, lower :: CBLASUpLo
upper = CBLASUpLo 121
lower = CBLASUpLo 122

cblasUpLo :: UpLoEnum -> CBLASUpLo
cblasUpLo Upper = CBLASUpLo 121
cblasUpLo Lower = CBLASUpLo 122

nonUnit, unit :: CBLASDiag
nonUnit = CBLASDiag 131
unit    = CBLASDiag 132

cblasDiag :: DiagEnum -> CBLASDiag
cblasDiag NonUnit = CBLASDiag 131
cblasDiag Unit    = CBLASDiag 132

leftSide, rightSide :: CBLASSide
leftSide  = CBLASSide 141
rightSide = CBLASSide 142

cblasSide :: SideEnum -> CBLASSide
cblasSide LeftSide  = CBLASSide 141
cblasSide RightSide = CBLASSide 142
