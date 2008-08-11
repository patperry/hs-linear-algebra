{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances #-}
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
    lowerFat,
    lowerTall,
    
    lowerU,
    lowerUFat,
    lowerUTall,
    
    upper,
    upperFat,
    upperTall,
    
    upperU,
    upperUFat,
    upperUTall,

    coerceTri,

    ) where

import BLAS.Internal ( checkSquare, checkFat, checkTall )
import BLAS.Matrix
import BLAS.Tensor
import BLAS.Types ( UpLo(..), Diag(..), flipUpLo )

import Unsafe.Coerce

data Tri a mn e = Tri UpLo Diag (a mn e)

-- | Coerce the shape type.
coerceTri :: Tri a mn e -> Tri a mn' e
coerceTri = unsafeCoerce

mapTri :: (a (m,n) e -> b (m,n) e) -> Tri a (m,n) e -> Tri b (m,n) e
mapTri f (Tri u d a) = Tri u d $ f a

fromBase :: UpLo -> Diag -> a (m,n) e -> Tri a (m,n) e
fromBase = Tri
        
toBase :: Tri a (m,n) e -> (UpLo, Diag, a (m,n) e)
toBase (Tri u d a) = (u,d,a)


lower :: (Matrix a) => a (n,n) e -> Tri a (n,n) e
lower a = checkSquare (shape a) $ Tri Lower NonUnit a

lowerFat :: (Matrix a) => a (m,n) e -> Tri a (m,m) e
lowerFat a = checkFat (shape a) $ Tri Lower NonUnit (unsafeCoerce a)

lowerTall :: (Matrix a) => a (m,n) e -> Tri a (m,n) e
lowerTall a = checkTall (shape a) $ Tri Lower NonUnit a


lowerU :: (Matrix a) => a (n,n) e -> Tri a (n,n) e
lowerU a = checkSquare (shape a) $ Tri Lower Unit a

lowerUFat :: (Matrix a) => a (m,n) e -> Tri a (m,m) e
lowerUFat a = checkFat (shape a) $ Tri Lower Unit (unsafeCoerce a)

lowerUTall :: (Matrix a) => a (m,n) e -> Tri a (m,n) e
lowerUTall a = checkTall (shape a) $ Tri Lower Unit a


upper :: (Matrix a) => a (n,n) e -> Tri a (n,n) e
upper a = checkSquare (shape a) $ Tri Upper NonUnit a

upperFat :: (Matrix a) => a (m,n) e -> Tri a (m,n) e
upperFat a = checkFat (shape a) $ Tri Upper NonUnit a

upperTall :: (Matrix a) => a (m,n) e -> Tri a (n,n) e
upperTall a = checkTall (shape a) $ Tri Upper NonUnit (unsafeCoerce a)


upperU :: (Matrix a) => a (n,n) e -> Tri a (n,n) e
upperU a = checkSquare (shape a) $ Tri Upper Unit a

upperUFat :: (Matrix a) => a (m,n) e -> Tri a (m,n) e
upperUFat a = checkFat (shape a) $ Tri Upper Unit a

upperUTall :: (Matrix a) => a (m,n) e -> Tri a (n,n) e
upperUTall a = checkTall (shape a) $ Tri Upper Unit (unsafeCoerce a)

      
instance Matrix a => Matrix (Tri a) where
    numRows (Tri Lower _ a) = numRows a
    numRows (Tri Upper _ a) = min (numRows a) (numCols a)
    
    numCols (Tri Lower _ a) = min (numRows a) (numCols a)
    numCols (Tri Upper _ a) = numCols a
    
    herm (Tri u d a) = Tri (flipUpLo u) d (herm a)


instance (Show (a (m,n) e), Matrix a) => Show (Tri a (m,n) e) where
    show (Tri u d a) =
        constructor ++ suffix ++ " (" ++ show a ++ ")"
        where
          constructor = case (u,d) of
              (Lower, NonUnit) -> "lower"
              (Lower, Unit   ) -> "lowerU"
              (Upper, NonUnit) -> "upper"
              (Upper, Unit   ) -> "upperU"

          suffix = case undefined of
                       _ | isSquare a -> ""
                       _ | isFat a    -> "Fat"
                       _              -> "Tall"
