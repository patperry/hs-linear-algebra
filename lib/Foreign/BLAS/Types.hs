{-# LANGUAGE GeneralizedNewtypeDeriving #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.BLAS.Types
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Foreign.BLAS.Types (
    BLASTrans(..),
    Trans(..),
    withTrans,
    
    BLASUplo(..),
    Uplo(..),
    withUplo,
    
    BLASDiag(..),
    Diag(..),
    withDiag,
    
    BLASSide(..),
    Side(..),
    withSide,
    
    LAInt(..),
    ) where

import Foreign
import Foreign.C.Types
import Foreign.C.String

newtype BLASTrans = BLASTrans CString
data Trans = NoTrans | Trans | ConjTrans deriving (Eq, Show)

withTrans :: Trans -> (BLASTrans -> IO a) -> IO a
withTrans trans f = flip withCString (f . BLASTrans) $ case trans of
    NoTrans   -> "N"
    Trans     -> "T"
    ConjTrans -> "C"

newtype BLASUplo = BLASUplo CString
data Uplo = Upper | Lower deriving (Eq, Show)

withUplo :: Uplo -> (BLASUplo -> IO a) -> IO a
withUplo uplo f = flip withCString (f . BLASUplo) $ case uplo of
    Upper -> "U"
    Lower -> "L"

newtype BLASSide = BLASSide CString
data Side = LeftSide | RightSide deriving (Eq, Show)

withSide :: Side -> (BLASSide -> IO a) -> IO a
withSide side f = flip withCString (f . BLASSide) $ case side of
    LeftSide  -> "L"
    RightSide -> "R"

newtype BLASDiag = BLASDiag CString
data Diag = NonUnit | Unit deriving (Eq, Show)

withDiag :: Diag -> (BLASDiag -> IO a) -> IO a
withDiag diag f = flip withCString (f . BLASDiag) $ case diag of
    NonUnit -> "N"
    Unit    -> "U"

newtype LAInt = LAInt CInt deriving (Eq, Show, Enum, Num, Ord, Storable)
