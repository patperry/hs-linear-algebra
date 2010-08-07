-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.BLAS.Types
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Foreign.BLAS.Types (
    BLASTrans,
    Trans(..),
    withTrans,
    
    BLASUplo,
    Uplo(..),
    withUplo,
    
    BLASDiag,
    Diag(..),
    withDiag,
    
    BLASSide,
    Side(..),
    withSide,
    
    ) where

import Foreign.C.String

type BLASTrans = CString
data Trans = NoTrans | Trans | ConjTrans deriving (Eq, Show)

withTrans :: Trans -> (BLASTrans -> IO a) -> IO a
withTrans trans = withCString $ case trans of
    NoTrans   -> "N"
    Trans     -> "T"
    ConjTrans -> "C"

type BLASUplo = CString
data Uplo = Upper | Lower deriving (Eq, Show)

withUplo :: Uplo -> (BLASUplo -> IO a) -> IO a
withUplo uplo = withCString $ case uplo of
    Upper -> "U"
    Lower -> "L"

type BLASSide = CString
data Side = LeftSide | RightSide deriving (Eq, Show)

withSide :: Side -> (BLASSide -> IO a) -> IO a
withSide side = withCString $ case side of
    LeftSide  -> "L"
    RightSide -> "R"

type BLASDiag = CString
data Diag = NonUnit | Unit deriving (Eq, Show)

withDiag :: Diag -> (BLASDiag -> IO a) -> IO a
withDiag diag = withCString $ case diag of
    NonUnit -> "N"
    Unit    -> "U"
