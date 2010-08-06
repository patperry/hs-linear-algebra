{-# LANGUAGE TypeFamilies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Types
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Numeric.LinearAlgebra.Types (
    HasVectorView(..),

    Order(..),
    flipOrder,

    Trans(..),
    flipTrans,

    UpLo(..),
    flipUpLo,

    Diag(..),

    Side(..),
    flipSide,
        
    ) where

-- | Types that can be viewed as vectors.
class HasVectorView (a :: * -> *) where
    type VectorView a :: * -> *

-- | Matrix element storage order.
data Order = RowMajor | ColMajor deriving (Eq, Show)

-- | Exchange @RowMajor@ and @ColMajor@.
flipOrder :: Order -> Order
flipOrder RowMajor = ColMajor
flipOrder ColMajor = RowMajor

-- | Transpose type.
data Trans = NoTrans | ConjTrans deriving (Eq, Show)

-- | Exchange @NoTrans@ and @ConjTrans@.
flipTrans :: Trans -> Trans
flipTrans NoTrans = ConjTrans
flipTrans ConjTrans = NoTrans

-- | Lower or upper triangular storage.
data UpLo = Upper | Lower deriving (Eq, Show)

-- | Exchange @Upper@ and @Lower@.
flipUpLo :: UpLo -> UpLo
flipUpLo Upper = Lower
flipUpLo Lower = Upper

-- | Diagonal storage.
data Diag = Unit | NonUnit deriving (Eq, Show)

-- | Multiplication side.
data Side = LeftSide | RightSide deriving (Eq, Show)

-- | Exchange @LeftSide@ and @RigthSide@.
flipSide :: Side -> Side
flipSide LeftSide  = RightSide
flipSide RightSide = LeftSide
