-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Types
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module BLAS.Types (
    OrderEnum(..),
    flipOrder,

    TransEnum(..),
    flipTrans,

    UpLoEnum(..),
    flipUpLo,

    DiagEnum(..),

    SideEnum(..),
    flipSide,

    ConjEnum(..),
    flipConj,
        
    ) where

-- | Matrix element storage order.
data OrderEnum = RowMajor | ColMajor deriving (Eq, Show)

-- | Exchange @RowMajor@ and @ColMajor@.
flipOrder :: OrderEnum -> OrderEnum
flipOrder RowMajor = ColMajor
flipOrder ColMajor = RowMajor

-- | Transpose type.
data TransEnum = NoTrans | ConjTrans deriving (Eq, Show)

-- | Exchange @NoTrans@ and @ConjTrans@.
flipTrans :: TransEnum -> TransEnum
flipTrans NoTrans = ConjTrans
flipTrans ConjTrans = NoTrans

-- | Lower or upper triangular storage.
data UpLoEnum = Upper | Lower deriving (Eq, Show)

-- | Exchange @Upper@ and @Lower@.
flipUpLo :: UpLoEnum -> UpLoEnum
flipUpLo Upper = Lower
flipUpLo Lower = Upper

-- | Diagonal storage.
data DiagEnum = Unit | NonUnit deriving (Eq, Show)

-- | Multiplication side.
data SideEnum = LeftSide | RightSide deriving (Eq, Show)

-- | Exchange @LeftSide@ and @RigthSide@.
flipSide :: SideEnum -> SideEnum
flipSide LeftSide  = RightSide
flipSide RightSide = LeftSide

-- | Vector conjugacy
data ConjEnum = NoConj | Conj deriving (Eq, Show)

-- | Exchange @NoConj@ and @Conj@.
flipConj :: ConjEnum -> ConjEnum
flipConj NoConj = Conj
flipConj Conj   = NoConj
