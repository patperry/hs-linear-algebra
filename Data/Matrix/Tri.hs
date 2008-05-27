{-# LANGUAGE MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Tri
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Tri (
    Tri(..),
    UpLo(..), Diag(..),

    fromBase,
    toBase,
    mapTri,

    lower,
    lowerU,
    upper,
    upperU,

    ) where

import qualified BLAS.Elem as E
import BLAS.Matrix
import BLAS.Tensor
import BLAS.Types ( UpLo(..), Diag(..), flipUpLo )

data Tri a nn e = Tri UpLo Diag e (a nn e)

mapTri :: (a (n,n) e -> b (n,n) e) -> Tri a (n,n) e -> Tri b (n,n) e
mapTri f (Tri u d n a) = Tri u d n $ f a

fromBase :: UpLo -> Diag -> e -> a (n,n) e -> Tri a (n,n) e
fromBase = Tri
        
toBase :: Tri a (n,n) e -> (UpLo, Diag, e, a (n,n) e)
toBase (Tri u d e a) = (u,d,e,a)

lower :: (Num e) => a (n,n) e -> Tri a (n,n) e
lower = Tri Lower NonUnit 1

lowerU :: (Num e) => a (n,n) e -> Tri a (n,n) e
lowerU = Tri Lower Unit 1

upper :: (Num e) => a (n,n) e -> Tri a (n,n) e
upper = Tri Upper NonUnit 1

upperU :: (Num e) => a (n,n) e -> Tri a (n,n) e
upperU = Tri Upper Unit 1
      
instance Matrix a => Matrix (Tri a) where
    numRows (Tri _ _ _ a) = numRows a
    numCols (Tri _ _ _ a) = numCols a
    herm    (Tri u d e a) = Tri (flipUpLo u) d (E.conj e) (herm a)
    
instance (Num e) => Scalable (Tri a nn) e where
    (*>) k (Tri u d e a) = Tri u d (k*e) a
