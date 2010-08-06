{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.CTypes
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Numeric.LinearAlgebra.CTypes (
    CBLASOrder,
    CBLASTrans,
    CBLASUplo,
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

    cblasTrans,
    cblasUplo,
    cblasDiag,
    cblasSide,
    ) where
        
import Numeric.LinearAlgebra.Types
        
newtype CBLASOrder = CBLASOrder Int deriving (Eq, Show)
newtype CBLASTrans = CBLASTrans Int deriving (Eq, Show)
newtype CBLASUplo  = CBLASUplo  Int deriving (Eq, Show)
newtype CBLASDiag  = CBLASDiag  Int deriving (Eq, Show)
newtype CBLASSide  = CBLASSide  Int deriving (Eq, Show)

rowMajor, colMajor :: CBLASOrder
rowMajor = CBLASOrder 101
colMajor = CBLASOrder 102

noTrans, trans, conjTrans :: CBLASTrans
noTrans   = CBLASTrans 111
trans     = CBLASTrans 112
conjTrans = CBLASTrans 113

cblasTrans :: Trans -> CBLASTrans
cblasTrans NoTrans   = CBLASTrans 111
cblasTrans ConjTrans = CBLASTrans 113

upper, lower :: CBLASUplo
upper = CBLASUplo 121
lower = CBLASUplo 122

cblasUplo :: Uplo -> CBLASUplo
cblasUplo Upper = CBLASUplo 121
cblasUplo Lower = CBLASUplo 122

nonUnit, unit :: CBLASDiag
nonUnit = CBLASDiag 131
unit    = CBLASDiag 132

cblasDiag :: Diag -> CBLASDiag
cblasDiag NonUnitDiag = CBLASDiag 131
cblasDiag UnitDiag    = CBLASDiag 132

leftSide, rightSide :: CBLASSide
leftSide  = CBLASSide 141
rightSide = CBLASSide 142

cblasSide :: Side -> CBLASSide
cblasSide LeftSide  = CBLASSide 141
cblasSide RightSide = CBLASSide 142
