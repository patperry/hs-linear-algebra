{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Tri.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Tri.Internal (
    Tri(..),
    UpLo(..), Diag(..),

    triFromBase,
    triToBase,
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

import Unsafe.Coerce

import BLAS.Internal ( checkSquare, checkFat, checkTall )
import Data.Matrix.Shaped
import BLAS.Tensor.Base
import BLAS.Types ( UpLo(..), Diag(..), flipUpLo )

data Tri a mn e = Tri UpLo Diag (a mn e)

-- | Coerce the shape type.
coerceTri :: Tri a mn e -> Tri a mn' e
coerceTri = unsafeCoerce

mapTri :: (a (m,n) e -> b (m,n) e) -> Tri a (m,n) e -> Tri b (m,n) e
mapTri f (Tri u d a) = Tri u d $ f a

triFromBase :: UpLo -> Diag -> a (m,n) e -> Tri a (m,n) e
triFromBase = Tri
        
triToBase :: Tri a (m,n) e -> (UpLo, Diag, a (m,n) e)
triToBase (Tri u d a) = (u,d,a)


lower :: (MatrixShaped a e) => a (n,n) e -> Tri a (n,n) e
lower a = checkSquare (shape a) $ Tri Lower NonUnit a

lowerFat :: (MatrixShaped a e) => a (m,n) e -> Tri a (m,m) e
lowerFat a = checkFat (shape a) $ Tri Lower NonUnit (unsafeCoerce a)

lowerTall :: (MatrixShaped a e) => a (m,n) e -> Tri a (m,n) e
lowerTall a = checkTall (shape a) $ Tri Lower NonUnit a


lowerU :: (MatrixShaped a e) => a (n,n) e -> Tri a (n,n) e
lowerU a = checkSquare (shape a) $ Tri Lower Unit a

lowerUFat :: (MatrixShaped a e) => a (m,n) e -> Tri a (m,m) e
lowerUFat a = checkFat (shape a) $ Tri Lower Unit (unsafeCoerce a)

lowerUTall :: (MatrixShaped a e) => a (m,n) e -> Tri a (m,n) e
lowerUTall a = checkTall (shape a) $ Tri Lower Unit a


upper :: (MatrixShaped a e) => a (n,n) e -> Tri a (n,n) e
upper a = checkSquare (shape a) $ Tri Upper NonUnit a

upperFat :: (MatrixShaped a e) => a (m,n) e -> Tri a (m,n) e
upperFat a = checkFat (shape a) $ Tri Upper NonUnit a

upperTall :: (MatrixShaped a e) => a (m,n) e -> Tri a (n,n) e
upperTall a = checkTall (shape a) $ Tri Upper NonUnit (unsafeCoerce a)


upperU :: (MatrixShaped a e) => a (n,n) e -> Tri a (n,n) e
upperU a = checkSquare (shape a) $ Tri Upper Unit a

upperUFat :: (MatrixShaped a e) => a (m,n) e -> Tri a (m,n) e
upperUFat a = checkFat (shape a) $ Tri Upper Unit a

upperUTall :: (MatrixShaped a e) => a (m,n) e -> Tri a (n,n) e
upperUTall a = checkTall (shape a) $ Tri Upper Unit (unsafeCoerce a)

      
instance MatrixShaped a e => BaseTensor (Tri a) (Int,Int) e where
    shape (Tri Lower _ a) = (numRows a, min (numRows a) (numCols a))
    shape (Tri Upper _ a) = (min (numRows a) (numCols a), numCols a)
    
    bounds a = ((0,0),(m-1,n-1)) where (m,n) = shape a
    
instance MatrixShaped a e => MatrixShaped (Tri a) e where
    herm (Tri u d a) = Tri (flipUpLo u) d (herm a)


instance (Show (a (m,n) e), MatrixShaped a e) => Show (Tri a (m,n) e) where
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
