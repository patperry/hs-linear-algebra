-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Types
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Types (
    OrderEnum(..),
    TransEnum(..),
    UpLoEnum(..),
    DiagEnum(..),
    SideEnum(..),
    
    flipOrder,
    flipTrans,
    flipUpLo,
    flipSide,
    ) where

-- | Matrix element storage order.
data OrderEnum = RowMajor | ColMajor deriving (Eq, Show, Enum)

-- | Transpose type.
data TransEnum = NoTrans | ConjTrans deriving (Eq, Show, Enum)

-- | Lower or upper triangular storage.
data UpLoEnum = Upper | Lower deriving (Eq, Show, Enum)

-- | Diagonal storage.
data DiagEnum = Unit | NonUnit deriving (Eq, Show, Enum)

-- | Multiplication side.
data SideEnum = LeftSide | RightSide deriving (Eq, Show, Enum)

-- | Exchange @RowMajor@ and @ColMajor@.
flipOrder :: OrderEnum -> OrderEnum
flipOrder RowMajor = ColMajor
flipOrder ColMajor = RowMajor

-- | Exchange @NoTrans@ and @ConjTrans@.
flipTrans :: TransEnum -> TransEnum
flipTrans NoTrans = ConjTrans
flipTrans ConjTrans = NoTrans

-- | Exchange @Upper@ and @Lower@.
flipUpLo :: UpLoEnum -> UpLoEnum
flipUpLo Upper = Lower
flipUpLo Lower = Upper
        
-- | Exchange @LeftSide@ and @RigthSide@.
flipSide :: SideEnum -> SideEnum
flipSide LeftSide  = RightSide
flipSide RightSide = LeftSide
