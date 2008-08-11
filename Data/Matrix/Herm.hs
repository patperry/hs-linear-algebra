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

import qualified BLAS.Elem as E
import BLAS.Matrix
import BLAS.Tensor
import BLAS.Types ( UpLo(..) )

data Herm a nn e = Herm UpLo e (a nn e)

coerceHerm :: Herm a mn e -> Herm a mn' e
coerceHerm = unsafeCoerce

mapHerm :: (a (n,n) e -> b (n,n) e) -> Herm a (n,n) e -> Herm b (n,n) e
mapHerm f (Herm u e a) = Herm u e $ f a

fromBase :: UpLo -> e -> a (n,n) e -> Herm a (n,n) e
fromBase = Herm
        
toBase :: Herm a (n,n) e -> (UpLo, e, a (n,n) e)
toBase (Herm u e a) = (u,e,a)

hermL :: (Num e) => a (n,n) e -> Herm a (n,n) e
hermL = Herm Lower 1

hermU :: (Num e) => a (n,n) e -> Herm a (n,n) e
hermU = Herm Upper 1
      
instance Matrix a => Matrix (Herm a) where
    numRows (Herm _ _ a) = numRows a
    numCols (Herm _ _ a) = numCols a
    herm    (Herm u e a) = Herm u (E.conj e) (unsafeCoerce a)
    
instance (Num e) => Scalable (Herm a nn) e where
    (*>) k (Herm u e a) = Herm u (k*e) a

instance (Show (a mn e), Show e, Num e) => Show (Herm a mn e) where
    show (Herm u k a) 
        | k /= 1 = "(" ++ show k ++ ") *> " ++ show (Herm u 1 a)
        | otherwise =
            constructor ++ " (" ++ show a ++ ")"
        where
          constructor = case u of
              Lower -> "hermL"
              Upper -> "hermU"
