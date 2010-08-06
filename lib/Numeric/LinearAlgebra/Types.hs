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

    Trans(..),
    isNoTrans,
    isConjTrans,
    flipTrans,

    UpLo(..),
    isLower,
    isUpper,
    flipUpLo,

    Diag(..),
    isUnitDiag,
    isNonUnitDiag,
    flipDiag,

    Side(..),
    isLeftSide,
    isRightSide,
    flipSide,
        
    ) where

-- | Types that can be viewed as vectors.
class HasVectorView (a :: * -> *) where
    type VectorView a :: * -> *

-- | Transpose type.
data Trans = NoTrans | ConjTrans deriving (Eq, Show)

-- | Indicates if the value is @NoTrans@.
isNoTrans :: Trans -> Bool
isNoTrans NoTrans = True
isNoTrans _       = False

-- | Indicates if the value is @ConjTrans@.
isConjTrans :: Trans -> Bool
isConjTrans ConjTrans = True
isConjTrans _         = False

-- | Exchange @NoTrans@ and @ConjTrans@.
flipTrans :: Trans -> Trans
flipTrans NoTrans = ConjTrans
flipTrans ConjTrans = NoTrans


-- | Lower or upper triangular storage.
data UpLo = Upper | Lower deriving (Eq, Show)

-- | Indicates if the value is @Lower@.
isLower :: UpLo -> Bool
isLower Lower = True
isLower _     = False

-- | Indicates if the value is @Upper@.
isUpper :: UpLo -> Bool
isUpper Upper = True
isUpper _     = False

-- | Exchange @Upper@ and @Lower@.
flipUpLo :: UpLo -> UpLo
flipUpLo Upper = Lower
flipUpLo Lower = Upper


-- | Diagonal storage.
data Diag = UnitDiag | NonUnitDiag deriving (Eq, Show)

-- | Indicates if the value is @UnitDiag@.
isUnitDiag :: Diag -> Bool
isUnitDiag UnitDiag = True
isUnitDiag _        = False

-- | Indicates if the value is @NonUnitDiag@.
isNonUnitDiag :: Diag -> Bool
isNonUnitDiag NonUnitDiag = True
isNonUnitDiag _           = False

-- | Exchange @UnitDiag@ and @NonUnitDiag@.
flipDiag :: Diag -> Diag
flipDiag UnitDiag    = NonUnitDiag
flipDiag NonUnitDiag = UnitDiag


-- | Multiplication side.
data Side = LeftSide | RightSide deriving (Eq, Show)

-- | Indicates if the value is @LeftSide@.
isLeftSide :: Side -> Bool
isLeftSide LeftSide = True
isLeftSide _        = False

-- | Indicates if the value is @RightSide@.
isRightSide :: Side -> Bool
isRightSide RightSide = True
isRightSide _         = False

-- | Exchange @LeftSide@ and @RigthSide@.
flipSide :: Side -> Side
flipSide LeftSide  = RightSide
flipSide RightSide = LeftSide
