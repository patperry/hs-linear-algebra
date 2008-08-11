{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Herm
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Herm (
    Herm(..),
    UpLo(..),

    fromBase,
    toBase,
    mapHerm,

    hermL,
    hermU,

    coerceHerm,

    ) where

import Unsafe.Coerce

import BLAS.Matrix
import BLAS.Types ( UpLo(..) )

data Herm a nn e = Herm UpLo (a nn e)

coerceHerm :: Herm a mn e -> Herm a mn' e
coerceHerm = unsafeCoerce

mapHerm :: (a (n,n) e -> b (n,n) e) -> Herm a (n,n) e -> Herm b (n,n) e
mapHerm f (Herm u a) = Herm u $ f a

fromBase :: UpLo -> a (n,n) e -> Herm a (n,n) e
fromBase = Herm
        
toBase :: Herm a (n,n) e -> (UpLo, a (n,n) e)
toBase (Herm u a) = (u,a)

hermL :: a (n,n) e -> Herm a (n,n) e
hermL = Herm Lower

hermU :: a (n,n) e -> Herm a (n,n) e
hermU = Herm Upper
      
instance Matrix a => Matrix (Herm a) where
    numRows (Herm _ a) = numRows a
    numCols (Herm _ a) = numCols a
    herm = coerceHerm
    
instance Show (a mn e) => Show (Herm a mn e) where
    show (Herm u a) = constructor ++ " (" ++ show a ++ ")"
      where
        constructor = case u of
            Lower -> "hermL"
            Upper -> "hermU"
