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

    hermFromBase,
    hermToBase,
    mapHerm,

    hermL,
    hermU,

    coerceHerm,

    ) where

import Unsafe.Coerce

import Data.Matrix.Class
import Data.Tensor.Class
import BLAS.Types ( UpLo(..) )

data Herm a nn e = Herm UpLo (a nn e)

coerceHerm :: Herm a mn e -> Herm a mn' e
coerceHerm = unsafeCoerce

mapHerm :: (a nn e -> b nn e) -> Herm a nn e -> Herm b nn e
mapHerm f (Herm u a) = Herm u $ f a

hermFromBase :: UpLo -> a (n,n) e -> Herm a (n,n) e
hermFromBase = Herm
        
hermToBase :: Herm a (n,n) e -> (UpLo, a (n,n) e)
hermToBase (Herm u a) = (u,a)

hermL :: a (n,n) e -> Herm a (n,n) e
hermL = Herm Lower

hermU :: a (n,n) e -> Herm a (n,n) e
hermU = Herm Upper
      
instance MatrixShaped a e => Shaped (Herm a) (Int,Int) e where
    shape  (Herm _ a) = (n,n)             where n = min (numRows a) (numCols a)
    bounds (Herm _ a) = ((0,0),(n-1,n-1)) where n = min (numRows a) (numCols a)
      
instance MatrixShaped a e => MatrixShaped (Herm a) e where
    herm = coerceHerm
    
instance Show (a mn e) => Show (Herm a mn e) where
    show (Herm u a) = constructor ++ " (" ++ show a ++ ")"
      where
        constructor = case u of
            Lower -> "hermL"
            Upper -> "hermU"
