-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Types
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module BLAS.Types (
    Order(..),
    Trans(..),
    UpLo(..),
    Diag(..),
    Side(..),
    
    flipOrder,
    flipTrans,
    flipUpLo,
    flipSide,
    ) where

-- | Matrix element storage order
data Order = RowMajor | ColMajor deriving (Eq, Show)

-- | Transpose type
data Trans = NoTrans | ConjTrans deriving (Eq,Show)

-- | Lower or upper triangular storage
data UpLo = Upper | Lower deriving (Eq,Show)

-- | Diagonal storage
data Diag = Unit | NonUnit deriving (Eq,Show)

-- | Multiplication side
data Side = LeftSide | RightSide deriving (Eq,Show)


flipOrder :: Order -> Order
flipOrder RowMajor = ColMajor
flipOrder ColMajor = RowMajor

flipTrans :: Trans -> Trans
flipTrans NoTrans = ConjTrans
flipTrans ConjTrans = NoTrans

flipUpLo :: UpLo -> UpLo
flipUpLo Upper = Lower
flipUpLo Lower = Upper
        
flipSide :: Side -> Side
flipSide LeftSide  = RightSide
flipSide RightSide = LeftSide
